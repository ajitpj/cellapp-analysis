"""Position-specific background and illumination correction.

The pipeline's stock corrections come from two blank wells: a FluoroBrite well
supplies the background map, a DMEM well the intensity map, and both are applied
to every position on the plate. They do not transfer. This module measures the
background from each position's *own* frames instead, and takes the illumination
from the difference of the two blanks.

The model behind everything here is

    I(x, y, t) = D + A(t) F(x, y) + S(x, y, t) F(x, y)

    D     camera offset, constant in space and time
    F     flat field: excitation profile x collection efficiency, mean 1.
          A property of the optics, so one per plate, not one per position
    A(t)  brightness of the medium. Depends on the well, the position and the
          time - this is what the blank-well map gets wrong
    S     the fluorophore in the cell, which is what we are after

giving   S = (I - B) / F   with the background   B(x, y, t) = D + A(t) F(x, y).

`B` is measured directly and never decomposed, which is what keeps the camera
offset `D` - the one quantity that cannot be pinned down from this data - off
the critical path entirely.

Usage, the parameters and what to check afterwards are in
SIGNAL_CORRECTION_README.md. Why it is built this way, what was measured, and
which alternatives were tried and rejected are in SIGNAL_CORRECTION_DESIGN.md.

Typical use::

    flat = flatfield_from_blank_pair(tifffile.imread(dmem_blank),
                                     tifffile.imread(fluorobrite_blank))
    with tifffile.TiffFile(position_stack) as tf:
        corr = estimate_position_correction(
            tf.series[0].asarray(out='memmap'),
            labels=tifffile.imread(instance_stack),
            stem='A08_s1', channel='GFP', flatfield=flat)
    corrected = corr.correct_measurements(raw, x, y, frame)

Nothing here imports matplotlib, and nothing writes into the data folder unless
you hand `PositionCorrection.save` a path that points there.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np
import numpy.typing as npt
import scipy.ndimage as ndi

__all__ = [
    # cell-free mask
    "cell_free_mask_from_labels",
    # background, per position and per frame
    "block_reduce_robust",
    "background_surface",
    "background_surfaces",
    "fit_background_to_flatfield",
    # flat field, one per plate
    "flatfield_from_blank_pair",
    # driver and container
    "estimate_position_correction",
    "PositionCorrection",
    # diagnostics
    "upsample",
    "radial_profile",
    "flatness",
    "centre_edge_ratio",
]

# Block size on the full-resolution fluorescence frame. 64 px on a 2048 frame
# gives a 32 x 32 grid: fine enough for a vignette, coarse enough that a block
# still holds several hundred cell-free pixels at 80% confluence.
DEFAULT_BLOCK = 64


# --------------------------------------------------------------------------
# cell-free masks
# --------------------------------------------------------------------------

def cell_free_mask_from_labels(labels: npt.NDArray,
                               shape: tuple[int, int] | None = None,
                               dilation: int = 21) -> npt.NDArray:
    """Cell-free pixels from an instance or semantic segmentation.

    The best detector available in this pipeline, and the cheapest: the
    segmentation is already on disk next to every position.

    Parameters
    ----------
    labels : 2D array
        One frame of the instance or semantic stack. Anything non-zero is a
        cell. cellaap writes these at half the fluorescence resolution.
    shape : (rows, cols), optional
        Shape of the fluorescence frame. `labels` is nearest-neighbour
        upsampled to it when the two differ, so pass the 2048-px frame shape
        for a 1024-px segmentation.
    dilation : int
        Side of the square used to grow the mask, in fluorescence pixels. The
        segmentation follows the cell body; out-of-focus haze reaches past it,
        and 21 px covers that on a 20x field. Set 0 to skip.

    Returns
    -------
    Boolean array the shape of the fluorescence frame, True where no cell is.
    """
    cell = np.asarray(labels) != 0
    if shape is not None and cell.shape != tuple(shape):
        fy = shape[0] / cell.shape[0]
        fx = shape[1] / cell.shape[1]
        if abs(fy - round(fy)) < 1e-9 and abs(fx - round(fx)) < 1e-9 and fy >= 1:
            # integer upsampling: repeat is far cheaper than ndi.zoom
            cell = np.repeat(np.repeat(cell, int(round(fy)), 0), int(round(fx)), 1)
        else:
            cell = ndi.zoom(cell.astype(np.uint8), (fy, fx), order=0).astype(bool)
    if dilation:
        # maximum_filter with a scalar size is the separable form of dilating
        # by a square, and gives a bit-identical result 100x faster - 0.03 s
        # against 4.9 s for a 121 px square on a 2048 frame, which is the
        # difference between this being usable per frame and not.
        cell = ndi.maximum_filter(cell, size=int(dilation), mode="nearest")
    return ~cell
def block_reduce_robust(image: npt.NDArray, block: int, how: str = "mean",
                        mask: npt.NDArray | None = None,
                        min_pixels: int = 64,
                        quantile: float = 0.10,
                        clip_sigma: float = 2.5,
                        clip_iters: int = 3) -> npt.NDArray:
    """Reduce an image to a grid of per-block statistics.

    `how` selects the estimator, and the choice matters more than it looks:

    ``mean``, ``std``, ``median``
        the plain thing, over the pixels `mask` selects.
    ``quantile``
        the `quantile`-th percentile of the block. Rejects cells without any
        mask, but sits low by roughly `z_q` noise sigmas - on the 20250213
        plate the 10th percentile runs 12-25 counts under the true cell-free
        median, which is 10-17% of the background. Use it only where a
        constant offset does not matter.
    ``clipped_mean``
        iterative sigma-clipping, rejecting the high side only. This is the
        one to use without a mask: it converges on the bulk of the medium
        pixels rather than a tail quantile, so it is unbiased where
        ``quantile`` is not, and it tolerates any cell coverage that leaves
        the medium as the majority of the block.

    Blocks with fewer than `min_pixels` usable pixels come back NaN.
    """
    a = np.asarray(image, np.float32)
    h, w = a.shape
    ny, nx = h // block, w // block
    if ny == 0 or nx == 0:
        raise ValueError(f"block {block} is larger than the {h}x{w} image")
    a = a[:ny * block, :nx * block]
    v = a.reshape(ny, block, nx, block).swapaxes(1, 2).reshape(ny, nx, -1)

    if mask is not None:
        m = np.asarray(mask, bool)[:ny * block, :nx * block]
        m = m.reshape(ny, block, nx, block).swapaxes(1, 2).reshape(ny, nx, -1)
        v = np.where(m, v, np.nan)
        enough = m.sum(-1) >= min_pixels
    else:
        enough = np.ones((ny, nx), bool)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN blocks
        if how == "mean":
            out = np.nanmean(v, -1)
        elif how == "median":
            out = np.nanmedian(v, -1)
        elif how == "std":
            out = np.nanstd(v, -1)
        elif how == "quantile":
            out = np.nanquantile(v, quantile, axis=-1)
        elif how == "clipped_mean":
            out = _clipped_mean(v, clip_sigma, clip_iters)
        else:
            raise ValueError(f"unknown reduction {how!r}")

    return np.where(enough, out, np.nan).astype(np.float32)


def _clipped_mean(v: npt.NDArray, sigma: float, iters: int) -> npt.NDArray:
    """Mean of each block after iteratively dropping the bright tail.

    One-sided on purpose: cells only ever add signal, so clipping the low side
    as well would throw away the medium pixels this is trying to measure.
    """
    x = np.where(np.isfinite(v), v, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for _ in range(iters):
            med = np.nanmedian(x, -1, keepdims=True)
            # median absolute deviation scaled to a Gaussian sigma, measured on
            # the low half only so cells cannot inflate the width they are
            # about to be compared against
            lo = np.where(x <= med, x, np.nan)
            mad = np.nanmedian(med - lo, -1, keepdims=True) * 1.4826
            mad = np.where(mad > 0, mad, np.nanstd(x, -1, keepdims=True))
            x = np.where(x <= med + sigma * mad, x, np.nan)
        return np.nanmean(x, -1)


def _fill_and_smooth(grid: npt.NDArray, smooth: float = 1.5) -> npt.NDArray:
    """Fill NaN grid points from their nearest measured neighbour, then smooth.

    Blocks go NaN where a position is too crowded to measure. Nearest-neighbour
    fill keeps the surface defined everywhere without inventing structure, and
    the Gaussian afterwards is what makes the result a smooth field rather than
    a mosaic. Returns a copy; the input is untouched.
    """
    g = np.array(grid, np.float32)
    bad = ~np.isfinite(g)
    if bad.all():
        raise ValueError("no usable blocks: every block was masked out or empty")
    if bad.any():
        idx = ndi.distance_transform_edt(bad, return_distances=False,
                                         return_indices=True)
        g = g[tuple(idx)]
    if smooth:
        g = ndi.gaussian_filter(g, smooth, mode="nearest")
    return g


def background_surface(image: npt.NDArray,
                       cell_free: npt.NDArray | None = None,
                       block: int = DEFAULT_BLOCK,
                       how: str | None = None,
                       smooth: float = 1.5,
                       full_resolution: bool = False) -> npt.NDArray:
    """The background under one frame: `D + A(t) F(x, y)`, cells excluded.

    Parameters
    ----------
    image : 2D array
        One fluorescence frame.
    cell_free : bool array, optional
        From `cell_free_mask_from_labels`. With a mask the estimator defaults
        to a per-block median; without one it defaults to `clipped_mean`,
        which rejects cells on its own.
    block, smooth
        Grid coarseness, and the Gaussian applied to the grid afterwards in
        units of grid points.
    full_resolution : bool
        Return a surface the size of `image` instead of the coarse grid.
        Bilinear from the grid - the surface has no structure finer than
        `block` either way, so keep the grid unless you need to subtract the
        surface from the frame pixel by pixel.
    """
    if how is None:
        how = "median" if cell_free is not None else "clipped_mean"
    grid = _fill_and_smooth(block_reduce_robust(image, block, how, mask=cell_free),
                            smooth)
    if not full_resolution:
        return grid
    return upsample(grid, np.asarray(image).shape)


def fit_background_to_flatfield(grid: npt.NDArray, flatfield: npt.NDArray
                                ) -> tuple[float, float]:
    """Fit `B = offset + level * F` to one measured background grid.

    The constrained alternative to smoothing the grid, and the better estimator
    when the cells are faint. Two free numbers instead of ~1000 grid points, so
    the fit averages over every block rather than tracking each one, which
    matters when the thing being subtracted is ten times the size of the signal
    left behind: a 1% wobble in the surface is a 10% error in the cell.

    It also removes two biases that a smoothed grid carries. Gaussian smoothing
    pulls the corners of the grid towards the brighter interior, and blocks too
    crowded to measure get filled from their neighbours - both land at the
    field edge, where the flat field is furthest from 1 and the error is
    multiplied.

    Requires a flat field from somewhere independent - `flatfield_from_blank_pair`
    is the intended source. Returns `(offset, level)`; `offset` is an estimate
    of the camera offset, and comparing it across frames and positions is a
    good check that the model is holding.
    """
    g = np.asarray(grid, np.float32)
    f = np.asarray(flatfield, np.float32)
    if f.shape != g.shape:
        f = ndi.zoom(f, np.array(g.shape) / np.array(f.shape), order=1,
                     mode="nearest")
    ok = np.isfinite(g)
    if ok.sum() < 10:
        raise ValueError("fewer than 10 usable blocks to fit a background to")
    A = np.stack([np.ones(int(ok.sum()), np.float32), f[ok]], axis=1)
    coef, *_ = np.linalg.lstsq(A, g[ok], rcond=None)
    return float(coef[0]), float(coef[1])


def background_surfaces(frames: Iterable[npt.NDArray],
                        cell_free: Iterable[npt.NDArray] | None = None,
                        block: int = DEFAULT_BLOCK,
                        how: str | None = None,
                        smooth: float = 1.5) -> npt.NDArray:
    """`background_surface` over a sequence of frames -> (n, ny, nx) grids."""
    masks = iter(cell_free) if cell_free is not None else None
    out = []
    for frame in frames:
        m = next(masks) if masks is not None else None
        out.append(background_surface(frame, m, block, how, smooth))
    return np.stack(out)
def flatfield_from_blank_pair(bright: npt.NDArray, dim: npt.NDArray,
                              block: int = DEFAULT_BLOCK,
                              smooth: float = 1.0,
                              min_contrast: float = 20.0) -> npt.NDArray:
    """Flat field from two blank wells of different brightness. **Preferred.**

    Two blanks imaged through the same optics are `D + A_bright F` and
    `D + A_dim F`. Their difference is `(A_bright - A_dim) F` - the camera
    offset cancels exactly, with no need to know it. That makes this the only
    background-derived route to `F` that does not rest on the one quantity
    none of the estimators here can pin down.

    On the 20250213 plate (DMEM as `bright`, FluoroBrite as `dim`) it gives a
    centre-to-edge ratio of 1.220 in GFP and 1.230 in Texas Red. Two
    independent checks say that is right where the single-blank maps are not:

    * the two channels agree to 1%, as an optical property must, while the
      single-blank maps disagree by 9% (1.190 vs 1.090 in GFP) purely because
      the wells held different amounts of medium;
    * feeding it back through `raw = D + (A + S) F` predicts a raw cell-signal
      falloff of 1.132, against 1.121 measured over 430k cells.

    Parameters
    ----------
    bright, dim : (t, y, x) or (y, x) arrays
        The two blank stacks. Order matters only for the sign; a bigger gap
        between them is better, so pass the brightest and the dimmest blanks
        available.
    min_contrast : float
        Refuse if the mean difference is smaller than this many counts - too
        small a gap and the difference is mostly noise.
    """
    gb = _mean_blank_surface(bright, block)
    gd = _mean_blank_surface(dim, block)
    diff = gb - gd
    contrast = float(np.nanmean(diff))
    if abs(contrast) < min_contrast:
        raise ValueError(
            f"the two blanks differ by only {contrast:.1f} counts on average; "
            f"their difference is mostly noise. Use blanks whose media differ "
            f"more - see SIGNAL_CORRECTION_DESIGN.md section 4 for the "
            f"alternatives when no usable blank pair exists.")
    if contrast < 0:                     # caller swapped them; harmless
        diff = -diff
    return _normalise_flatfield(diff, smooth)


def _mean_blank_surface(blank: npt.NDArray, block: int) -> npt.NDArray:
    """Mean block-reduced surface of a cell-free stack."""
    a = np.asarray(blank, np.float32)
    frames = a if a.ndim == 3 else a[None]
    return np.stack([block_reduce_robust(f, block, "clipped_mean")
                     for f in frames]).mean(axis=0)
def _normalise_flatfield(grid: npt.NDArray, smooth: float) -> npt.NDArray:
    """Fill, smooth, normalise to mean 1, and refuse anything unphysical."""
    g = _fill_and_smooth(grid, smooth)
    m = float(np.mean(g))
    if not np.isfinite(m) or m <= 0:
        raise ValueError("flat field has a non-positive mean; the darkfield is "
                         "probably larger than the background it was removed from")
    g = g / m
    if g.min() <= 0:
        raise ValueError(f"flat field reaches {g.min():.3f}; dividing by it would "
                         "flip the sign of the signal. The darkfield is too large.")
    return g.astype(np.float32)
@dataclass
class PositionCorrection:
    """Everything needed to correct one channel of one position.

    The background is stored as one fixed shape with two numbers per frame::

        B(x, y, t) = background_offset[t] + background_level[t] * background_shape(x, y)

    The offset term is what lets the `"flatfield"` background model be stored
    exactly: there `background_shape` *is* the flat field, `background_level`
    is the medium brightness `A(t)`, and `background_offset` is the camera
    offset the fit found - per frame, never assumed. The `"grid"` model leaves
    the offset at zero and puts everything in the shape. Either way nothing
    downstream has to know which was used.

        corrected = (I - B(t)) / flatfield

    Three small arrays replace the two full-size map stacks the pipeline
    carries today: a 32x32 shape, a 32x32 flat field, and two numbers a frame.
    """

    stem: str
    channel: str
    background_shape: npt.NDArray             # (gy, gx), mean 1
    background_level: npt.NDArray             # (n_frames,), counts
    frame_shape: tuple[int, int]
    flatfield: npt.NDArray | None = None      # (gy, gx), mean 1; None = uncorrected
    background_offset: npt.NDArray | None = None  # (n_frames,), counts
    block: int = DEFAULT_BLOCK
    darkfield: float | None = None            # provenance only, if one was used
    background_grids: npt.NDArray | None = None   # measured, (n_sampled, gy, gx)
    sampled_frames: npt.NDArray | None = None
    diagnostics: dict = field(default_factory=dict)

    def __post_init__(self):
        self.background_shape = np.asarray(self.background_shape, np.float32)
        self.background_level = np.asarray(self.background_level, np.float32)
        self.frame_shape = tuple(int(s) for s in self.frame_shape)
        if self.background_offset is None:
            self.background_offset = np.zeros_like(self.background_level)
        else:
            self.background_offset = np.asarray(self.background_offset, np.float32)
        if self.flatfield is not None:
            self.flatfield = np.asarray(self.flatfield, np.float32)

    # -- the flat field is set later, after pooling the plate --------------
    def set_flatfield(self, flatfield: npt.NDArray,
                      darkfield: float | None = None) -> "PositionCorrection":
        """Attach a flat field, resampling it to this position's grid."""
        f = np.asarray(flatfield, np.float32)
        if f.shape != self.background_shape.shape:
            f = ndi.zoom(f, np.array(self.background_shape.shape) / np.array(f.shape),
                         order=1, mode="nearest")
        self.flatfield = (f / f.mean()).astype(np.float32)
        if darkfield is not None:
            self.darkfield = float(darkfield)
        self.diagnostics["flatfield_centre_edge"] = centre_edge_ratio(self.flatfield)
        return self

    def _flat_grid(self) -> npt.NDArray:
        if self.flatfield is None:
            return np.ones_like(self.background_shape)
        return self.flatfield

    # -- applying it ------------------------------------------------------
    def background(self, frame: int, full_resolution: bool = True) -> npt.NDArray:
        """The background under one frame."""
        b = (float(self.background_offset[frame])
             + float(self.background_level[frame]) * self.background_shape)
        return upsample(b, self.frame_shape) if full_resolution else b

    def flat(self, full_resolution: bool = True) -> npt.NDArray:
        """The flat field, optionally at frame resolution."""
        g = self._flat_grid()
        return upsample(g, self.frame_shape) if full_resolution else g

    def apply(self, image: npt.NDArray, frame: int) -> npt.NDArray:
        """`(I - B(t)) / F` for a whole frame."""
        return ((np.asarray(image, np.float32) - self.background(frame))
                / self.flat())

    def _grid_index(self, x, y):
        gy, gx = self.background_shape.shape
        j = np.clip((np.asarray(x, float) / self.frame_shape[1] * gx).astype(int),
                    0, gx - 1)
        i = np.clip((np.asarray(y, float) / self.frame_shape[0] * gy).astype(int),
                    0, gy - 1)
        return i, j

    def subtract_background(self, raw, x, y, frame) -> npt.NDArray:
        """`raw - B(x, y, t)` for per-cell measurements. Pass 1 of the workflow.

        Useful on its own when no flat field is available: the background is
        the larger of the two corrections by an order of magnitude here.
        """
        i, j = self._grid_index(x, y)
        f = np.clip(np.asarray(frame, int), 0, len(self.background_level) - 1)
        bg = (self.background_offset[f]
              + self.background_level[f] * self.background_shape[i, j])
        return np.asarray(raw, float) - bg

    def correct_measurements(self, raw, x, y, frame) -> npt.NDArray:
        """Fully corrected per-cell signal, without touching the images.

        The route into the existing tables: `cellaap_analysis` already stores a
        raw mean per cell per frame and the centroid it came from, so this
        recomputes the corrected signal from `*_analysis.xlsx` alone.

        Sampling the fields at the centroid rather than averaging them over the
        mask - which is what `<ch>_bkg_corr` and `<ch>_int_corr` do - is fine
        here: both are smooth on a scale of `block` pixels, far larger than a
        cell, so the two agree to well under a percent.
        """
        i, j = self._grid_index(x, y)
        return self.subtract_background(raw, x, y, frame) / self._flat_grid()[i, j]

    # -- persistence ------------------------------------------------------
    def save(self, path: str | Path) -> Path:
        """Write to a `.npz`, with a human-readable `.json` beside it."""
        path = Path(path).with_suffix(".npz")
        path.parent.mkdir(parents=True, exist_ok=True)
        meta = {"stem": self.stem, "channel": self.channel,
                "frame_shape": list(self.frame_shape), "block": self.block,
                "darkfield": self.darkfield, "diagnostics": self.diagnostics}
        np.savez_compressed(
            path,
            background_shape=self.background_shape,
            background_level=self.background_level,
            background_offset=self.background_offset,
            flatfield=(self.flatfield if self.flatfield is not None
                       else np.array([], np.float32)),
            background_grids=(self.background_grids
                              if self.background_grids is not None
                              else np.array([], np.float32)),
            sampled_frames=(self.sampled_frames if self.sampled_frames is not None
                            else np.array([], int)),
            meta=json.dumps(meta))
        path.with_suffix(".json").write_text(json.dumps(
            {**meta,
             "background_level_first": float(self.background_level[0]),
             "background_level_last": float(self.background_level[-1]),
             "flatfield_applied": self.flatfield is not None}, indent=2))
        return path

    @classmethod
    def load(cls, path: str | Path) -> "PositionCorrection":
        z = np.load(Path(path).with_suffix(".npz"), allow_pickle=False)
        meta = json.loads(str(z["meta"]))
        flat, grids, frames = z["flatfield"], z["background_grids"], z["sampled_frames"]
        return cls(stem=meta["stem"], channel=meta["channel"],
                   background_shape=z["background_shape"],
                   background_level=z["background_level"],
                   background_offset=z["background_offset"],
                   frame_shape=tuple(meta["frame_shape"]),
                   flatfield=flat if flat.size else None,
                   block=meta["block"], darkfield=meta.get("darkfield"),
                   background_grids=grids if grids.size else None,
                   sampled_frames=frames if frames.size else None,
                   diagnostics=meta.get("diagnostics", {}))


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def estimate_position_correction(stack,
                                 *,
                                 stem: str = "",
                                 channel: str = "",
                                 labels: npt.NDArray | None = None,
                                 flatfield: npt.NDArray | None = None,
                                 darkfield: float | None = None,
                                 background_model: str = "grid",
                                 n_frames: int = 24,
                                 block: int = DEFAULT_BLOCK,
                                 dilation: int = 121,
                                 min_usable_blocks: float = 0.05,
                                 how: str | None = None,
                                 smooth: float = 1.5,
                                 keep_grids: bool = True) -> PositionCorrection:
    """Estimate the background correction for one position and channel.

    Measures a background surface on `n_frames` evenly spaced frames, splits
    those into one shape and a per-frame level, and interpolates the level over
    the whole movie.

    Parameters
    ----------
    stack : (t, y, x) array or anything indexable per frame
        The fluorescence movie. A `tifffile` memory map is fine and preferred -
        only the sampled frames are read.
    labels : (t, ly, lx) array, optional
        The instance or semantic segmentation, for masking cells out. Strongly
        recommended; without it the estimator falls back to `clipped_mean`,
        which is decent but not as good.
    flatfield : 2D array, optional
        The illumination correction, normally from `flatfield_from_blank_pair`
        and shared by the whole plate. It is used twice: to divide out the
        illumination, and - under `background_model="flatfield"` - as the shape
        the background is fitted to. Without one the result corrects background
        only, using the `"grid"` model.
    darkfield : float, optional
        Recorded for provenance if the `flatfield` you pass was derived using
        one. Nothing in this function divides by it.
    background_model : {"grid", "flatfield"}
        How the measured block values become a surface. `"grid"`, the default,
        fills and smooths the measured grid. `"flatfield"` instead fits
        `offset + level * F` per frame - two numbers rather than ~1000 grid
        points, no smoothing bias, and a fitted camera offset for free.

        The constrained fit ought to win when the cells are faint, and on the
        20250213 plate it does not: measured over 44k cells its residual radial
        trend is 25% against the grid's 9% (GFP). The reason is that the
        background is genuinely not proportional to `F` - out-of-focus haze
        follows local cell density, which has its own shape - and the grid can
        follow that where two parameters cannot. Reach for `"flatfield"` when
        so few blocks survive masking that the grid is mostly interpolation, or
        when you want the fitted offset as a diagnostic.
    dilation : int
        How far from the segmentation to stay when measuring background, in
        fluorescence pixels, and **the single most important parameter here**.
        Out-of-focus haze around a cell is not background, and subtracting it
        removes signal: on the 20250213 plate the medium 10 px from a cell
        reads 132 counts, 30 px out 127, against ~13-16 counts of actual cell
        signal. Widening the exclusion from 21 to 121 px raises the recovered
        GFP signal from 12.0 to 15.9 counts and drops the residual radial
        trend from 34% to 9%.

        The cost is usable blocks - 75% at 21 px, 10% at 121 px on this plate -
        so the dilation is backed off automatically, per frame, when too few
        blocks survive. `min_usable_blocks` sets that floor.
    n_frames : int
        Frames to sample. The background drifts smoothly, so 24 over a
        137-frame movie resolves it comfortably.

    Returns
    -------
    PositionCorrection
    """
    if background_model == "flatfield" and flatfield is None:
        raise ValueError('background_model="flatfield" needs a flatfield to fit to')
    if background_model not in ("flatfield", "grid"):
        raise ValueError(f"unknown background_model {background_model!r}")

    n_total = len(stack)
    idx = np.unique(np.linspace(0, n_total - 1, min(n_frames, n_total)).astype(int))
    shape = np.asarray(stack[int(idx[0])]).shape

    reduction = how or ("median" if labels is not None else "clipped_mean")
    measured, unusable, used_dilation = [], [], []
    for k in idx:
        img = np.asarray(stack[int(k)], np.float32)
        if labels is None:
            grid = block_reduce_robust(img, block, reduction)
            measured.append(grid)
            unusable.append(float(np.mean(~np.isfinite(grid))))
            used_dilation.append(0)
            continue
        # Stay as far from cells as this frame can afford. A wide exclusion is
        # what keeps out-of-focus haze out of the background, but a crowded
        # frame cannot spare it, and a background measured from three surviving
        # blocks is worse than one measured closer in.
        lab = labels[int(k)]
        for d in _dilation_ladder(dilation):
            grid = block_reduce_robust(img, block, reduction,
                                       mask=cell_free_mask_from_labels(lab, shape, d))
            usable = float(np.mean(np.isfinite(grid)))
            if usable >= min_usable_blocks or d == 0:
                break
        measured.append(grid)
        unusable.append(1.0 - usable)
        used_dilation.append(d)
    measured = np.stack(measured)

    if background_model == "flatfield":
        # The shape is fixed and known; only two numbers a frame are fitted.
        bshape = np.asarray(flatfield, np.float32)
        if bshape.shape != measured.shape[1:]:
            bshape = ndi.zoom(bshape, np.array(measured.shape[1:]) / np.array(bshape.shape),
                              order=1, mode="nearest")
        bshape = (bshape / bshape.mean()).astype(np.float32)
        fits = [fit_background_to_flatfield(g, bshape) for g in measured]
        offset_at = np.array([f[0] for f in fits])
        level_at = np.array([f[1] for f in fits])
        grids = np.stack([o + l * bshape for o, l in fits])
    else:
        # The shape is the mean measured surface; the level is the
        # least-squares amplitude that best explains each frame with it, which
        # is steadier than a frame mean when blocks had to be filled in.
        grids = np.stack([_fill_and_smooth(g, smooth) for g in measured])
        mean_surface = grids.mean(axis=0)
        bshape = (mean_surface / mean_surface.mean()).astype(np.float32)
        w = float((bshape ** 2).sum())
        level_at = np.array([float((g * bshape).sum() / w) for g in grids])
        offset_at = np.zeros_like(level_at)

    frames = np.arange(n_total)
    level = np.interp(frames, idx, level_at).astype(np.float32)
    offset = np.interp(frames, idx, offset_at).astype(np.float32)

    total = offset_at + level_at * float(bshape.mean())
    corr = PositionCorrection(
        stem=stem, channel=channel,
        background_shape=bshape, background_level=level,
        background_offset=offset,
        frame_shape=tuple(shape), block=block, darkfield=darkfield,
        background_grids=grids if keep_grids else None,
        sampled_frames=idx if keep_grids else None,
        diagnostics={
            "n_frames_sampled": int(len(idx)),
            "background_model": background_model,
            "background_shape_centre_edge": centre_edge_ratio(bshape),
            "background_mean_range": [float(total.min()), float(total.max())],
            "background_drift_percent": float(100 * (total.max() - total.min())
                                              / max(total.mean(), 1e-9)),
            "fitted_offset_range": [float(offset_at.min()), float(offset_at.max())],
            "unusable_block_fraction": float(np.mean(unusable)),
            "masked": labels is not None,
            "dilation_requested": int(dilation),
            "dilation_used": [int(min(used_dilation)), int(max(used_dilation))],
        })
    if flatfield is not None:
        corr.set_flatfield(flatfield, darkfield)
    return corr


# --------------------------------------------------------------------------
# diagnostics
# --------------------------------------------------------------------------

def _dilation_ladder(dilation: int) -> list[int]:
    """Exclusion widths to try, widest first, ending at no mask at all."""
    d, out = int(dilation), []
    while d > 10:
        out.append(d)
        d //= 2
    return out + [0]


def upsample(grid: npt.NDArray, shape: tuple[int, int]) -> npt.NDArray:
    """Bilinear resize of a coarse grid to a full frame."""
    out = ndi.zoom(grid, (shape[0] / grid.shape[0], shape[1] / grid.shape[1]),
                   order=1, mode="nearest")
    if out.shape != tuple(shape):              # zoom can land a pixel off
        out = np.pad(out, [(0, max(0, s - o)) for s, o in zip(shape, out.shape)],
                     mode="edge")[:shape[0], :shape[1]]
    return out.astype(np.float32)


def radial_profile(values: npt.NDArray, x: npt.NDArray, y: npt.NDArray,
                   shape: tuple[int, int], n_bins: int = 8
                   ) -> tuple[npt.NDArray, npt.NDArray]:
    """Median of `values` against normalised distance from the field centre.

    Radius runs 0 at the centre to 1 at the mid-edge, as in the aggregate
    notebook, so the numbers are comparable with the ones there. Returns
    `(bin_centres, medians)`, the medians divided by the innermost bin.
    """
    r = np.hypot(np.asarray(y, float) / (shape[0] - 1) - .5,
                 np.asarray(x, float) / (shape[1] - 1) - .5) / .5
    edges = np.linspace(0, float(r.max()) + 1e-9, n_bins + 1)
    idx = np.digitize(r, edges) - 1
    v = np.asarray(values, float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.array([np.nanmedian(v[idx == k]) if (idx == k).sum() > 5 else np.nan
                        for k in range(n_bins)])
    return (edges[:-1] + edges[1:]) / 2, med / med[np.isfinite(med)][0]


def flatness(profile: npt.NDArray) -> dict:
    """How flat is a radial profile that should be 1 everywhere?

    RMS deviation from 1 across all bins is the verdict, for the reasons the
    aggregate notebook sets out: `max - min` scores a 12% falloff and a 12%
    rise the same, and edge/centre misses a profile that swings in the middle
    and happens to land near 1 at the last bin. The other two are kept as
    description.
    """
    p = np.asarray(profile, float)
    p = p[np.isfinite(p)]
    d = np.diff(p)
    return {"rms_deviation_percent": float(100 * np.sqrt(np.mean((p - 1) ** 2))),
            "spread_percent": float(100 * (p.max() - p.min())),
            "edge_over_centre": float(p[-1] / p[0]),
            "monotonic": bool((d >= -0.03).all() or (d <= 0.03).all())}


def centre_edge_ratio(grid: npt.NDArray) -> float:
    """Ratio of the middle of a grid to its rim - a one-number vignette depth."""
    g = np.asarray(grid, float)
    if g.ndim == 1:                       # already a radial profile
        return float(g[0] / g[-1])
    ny, nx = g.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    r = np.hypot(yy / (ny - 1) - .5, xx / (nx - 1) - .5) / .5
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return float(np.nanmedian(g[r < .25]) / np.nanmedian(g[r > .85]))
