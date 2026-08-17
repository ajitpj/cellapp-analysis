"""Position-specific background and illumination corrections.

The pipeline's existing corrections come from two blank wells: a FluoroBrite
well supplies the background map, a DMEM well the intensity map, and both are
applied to every position on the plate. This module estimates the same two
quantities from each position's *own* frames instead, because the blank-well
maps do not transfer.

The model
---------
Every measurement in a fluorescence frame is

    I(x, y, t) = D + A(t) F(x, y) + S(x, y, t) F(x, y)

    D     camera offset (dark level), constant in space and time
    F     flat field: excitation profile x collection efficiency, normalised
          to mean 1. A property of the optics, not of the well
    A(t)  brightness of the medium (autofluorescence, stray light). Depends on
          the well, the position, and the time
    S     what we actually want: the fluorophore in the cell

so the corrected signal is

    S = (I - B) / F     with the background     B(x, y, t) = D + A(t) F(x, y)

`B` is the "background map" of the existing pipeline and `F` its "intensity
map". On the 20250213 plate the linear part of that model - a blank well is
`D + A(t) F` and nothing else - reproduces the measured frames to 0.1% of their
mean, so the model is not the weak point.

Why the blank-well maps fail
----------------------------
Three separate errors, all of which this module avoids:

1. **The background level is per-position.** Across the 13 sample positions of
   that plate the cell-free GFP level runs 95-105 counts and the FluoroBrite
   blank sits at 90 - below every one of them. Subtracting the blank leaves a
   position-dependent residual.

2. **The background level drifts, and the blank does not drift with it.** Over
   the 137-frame movie the sample background rises ~11% (GFP) and ~14-21%
   (Texas Red) while the FluoroBrite well stays flat. A per-frame map built
   from the wrong well cannot track that.

3. **The intensity map keeps the camera offset in it, and how much that matters
   depends on how bright the blank was.** `map = blank / blank.max()` estimates
   `(D + A F) / (D + A F).max()`, which is `F` flattened towards 1 by the
   fraction `D / (D + A)`. The dim FluoroBrite well gives a 6% centre-to-edge
   falloff, the bright DMEM well 14%, and the cells themselves report 11-12% -
   three different answers to the same question, differing only in how much
   medium was in the well. Dividing by a map built without removing `D` is why
   the corrected field in the aggregate notebook came out *less* flat than the
   raw one.

Why nothing here depends on the camera offset
---------------------------------------------
Point 3 makes `D` look essential, and it is - for anything that turns a single
*background* into a flat field. It is also the one quantity that cannot be
measured from this data: `darkfield_from_ptc` extrapolates the sensor's
noise-versus-signal line below the measured range and lands anywhere between
-35 and +25 counts, and a flat field derived from a background goes as
`1 / (background - D)`, so that scatter is fatal.

The way out is to not need it:

* **The background** is subtracted as the *measured* surface `B`, never
  reconstructed from `D` and `A(t)`. Nothing is decomposed, so nothing needs
  `D`.
* **The flat field** comes from the *difference of two blank wells*
  (`flatfield_from_blank_pair`). Both are `D + A F` with the same `D`, so
  subtracting one from the other cancels the offset exactly and leaves
  `(A_bright - A_dim) F`. Where there is no second blank and the cells are
  bright enough, `flatfield_from_signal` fits `F` to the cells instead, which
  also carries no additive term.

`darkfield_from_flatfield` then reports `D` as a by-product, and it is worth
reading as a check: on the 20250213 plate the DMEM and FluoroBrite wells imply
53 and 61 counts in GFP, 76 and 79 in Texas Red - two very different wells
agreeing, which is what says the model holds. `flatfield_from_blank` and
`flatfield_from_surfaces` remain for the one case where `D` is genuinely
known: a measured dark frame.

Finding cell-free regions
-------------------------
`cell_free_mask_from_phase` implements the minimum-variance idea and its
relatives, but on this data they barely work: at 64-78% cell coverage the
lowest-variance 10% of 64-pixel blocks are only 8-26% strictly cell-free,
a factor 2-4 over picking blocks at random. Run `benchmark_cell_free_masks`
against the segmentation before trusting any of them on new data.

So the estimators here are built not to need cell-free regions:

* `background_surface` takes a robust low-side location estimate inside each
  block. It needs *some* cell-free pixels per block, not a cell-free region,
  and it degrades smoothly as coverage rises.
* `temporal_background` uses a low quantile through time at each pixel. It
  needs every pixel to be cell-free at some point in the movie, which moving
  cells provide, and needs no cell-free region in any single frame.
* `flatfield_from_signal` fits the flat field to the cells, from the assumption
  that cells are on average equally bright everywhere in the field. It needs no
  background estimate and no cell-free pixel at all - but it does need the
  cells to stand out, and on this plate they do not (see below).

When the segmentation is already on disk - which it is, everywhere this module
runs - `cell_free_mask_from_labels` beats every image-based detector for free,
and is the default.

The recommended workflow
------------------------
::

    # the flat field: one for the plate, from the two blank wells
    flat = flatfield_from_blank_pair(tifffile.imread(dmem_well),
                                     tifffile.imread(fluorobrite_well))

    # the background: per position, per frame, from that position's own frames
    for pos in positions:
        stack  = tifffile.TiffFile(pos.image).series[0].asarray(out='memmap')
        labels = tifffile.imread(pos.instance)
        corrections[pos.key] = estimate_position_correction(
            stack, labels=labels, stem=pos.key, channel='GFP', flatfield=flat)

    corrected = corrections[key].correct_measurements(raw, x, y, frame)

The flat field is pooled over the plate on purpose: it is a property of the
optics, so every position measures the same thing and pooling is what makes it
precise. The background is emphatically not pooled - that is the whole point.

With no blank wells at all, drop the first line and fit the flat field to the
cells instead, in a second pass over `subtract_background` output. That route
needs the cells to be well above their local background; where they are not it
refuses rather than returning noise, and correcting the background alone is
still worth doing - it is the larger of the two errors here by an order of
magnitude.

What it buys, measured
----------------------
Checked on the 20250213 plate against the 431k per-cell measurements already in
the `*_analysis.xlsx` files - 13 positions, both channels. The honest metric is
the residual spatial bias in *counts*: how much a cell's reported signal moves
across the field purely because of where it sat. Percentages are not comparable
between stages whose medians differ by tenfold.

===========================  ===========  =============  =========
stage                        median       radial swing   RMS
===========================  ===========  =============  =========
GFP raw                      139.6         15.4 counts    8.2
GFP blank-well maps           24.9          4.2 counts    2.5
GFP this module               15.9          4.0 counts    1.4
Texas Red raw                164.7         20.5 counts   11.7
Texas Red blank-well maps     14.7          5.7 counts    1.6
Texas Red this module          9.9          2.5 counts    1.9
===========================  ===========  =============  =========

The bigger correction is not the spatial one, though - it is the level. The
blank-well map subtracts 116 counts (GFP) from every position, where the
measured cell-free background of the positions themselves runs 118-128 and
climbs 9% through the movie. So the existing pipeline under-subtracts by 2-12
counts in GFP and by up to 21 in Texas Red, *position by position*, against a
cell signal of 10-16 counts. That is why it reports a GFP median of 24.9 where
this module reports 15.9: most of the difference is unremoved background, and
how much of it there is depends on the position.

Cost
----
Everything works on block-reduced grids (2048 -> 32 x 32 by default) and on a
subsample of frames: ~14 s per position per channel on this plate, most of it
reading frames off the share, and a saved correction is 50 kB against a
full-size map stack.

Nothing here imports matplotlib, and nothing here writes into the data folder
unless you call `PositionCorrection.save` with a path that points there.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import numpy.typing as npt
import scipy.ndimage as ndi

__all__ = [
    # masks
    "cell_free_mask_from_labels",
    "cell_free_mask_from_phase",
    "benchmark_cell_free_masks",
    # background
    "block_reduce_robust",
    "background_surface",
    "background_surfaces",
    "fit_background_to_flatfield",
    "temporal_background",
    # flat field
    "flatfield_from_blank_pair",
    "flatfield_from_signal",
    "flatfield_from_surfaces",
    "flatfield_from_blank",
    # camera offset - a diagnostic here, never on the critical path
    "darkfield_from_flatfield",
    "darkfield_from_ptc",
    "darkfield_from_shape_consistency",
    "level_regression",
    # driver and container
    "PositionCorrection",
    "estimate_position_correction",
    # diagnostics
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


def cell_free_mask_from_phase(phase: npt.NDArray,
                              statistic: str = "std",
                              block: int = DEFAULT_BLOCK,
                              keep: float = 0.25) -> npt.NDArray:
    """Cell-free pixels from the texture of the phase image.

    For when no segmentation exists. Measure its performance on your own data
    with `benchmark_cell_free_masks` before relying on it - on the 20250213
    plate, at 64-78% cell coverage, the blocks it keeps are only 8-26%
    strictly cell-free.

    Parameters
    ----------
    phase : 2D array
        One phase frame.
    statistic : {"std", "gradient", "tophat"}
        What "looks like medium" means.

        ``std``      per-block standard deviation - the minimum-variance idea.
                     Cheapest, and as good as the others here.
        ``gradient`` mean Sobel magnitude per block. Slightly more selective
                     for cell edges, slightly more sensitive to the phase halo,
                     which is bright *next to* cells rather than on them.
        ``tophat``   mean |frame - local mean| per block. Responds to a cell
                     body that is uniformly offset from the medium but flat
                     inside, which the two above miss; the most useful of the
                     three on large spread cells.
    block : int
        Block side in pixels.
    keep : float
        Fraction of blocks to call cell-free, lowest statistic first.

    Returns
    -------
    Boolean array the shape of `phase`, True in the kept blocks.
    """
    a = np.asarray(phase, np.float32)
    stat = _phase_texture(a, statistic, block)
    keep_blocks = stat <= np.quantile(stat, keep)
    tiled = np.repeat(np.repeat(keep_blocks, block, 0), block, 1)
    out = np.zeros(a.shape, bool)              # pad back the trimmed edge blocks
    out[:tiled.shape[0], :tiled.shape[1]] = tiled
    return out


def _phase_texture(a: npt.NDArray, statistic: str, block: int) -> npt.NDArray:
    """Per-block texture statistic used by the phase-based mask."""
    if statistic == "std":
        return block_reduce_robust(a, block, "std")
    if statistic == "gradient":
        from skimage.filters import sobel
        return block_reduce_robust(sobel(a), block, "mean")
    if statistic == "tophat":
        # uniform_filter over ~2 blocks is a cheap stand-in for a grey opening
        # with a large structuring element, which costs far more.
        return block_reduce_robust(np.abs(a - ndi.uniform_filter(a, 2 * block + 1)),
                                   block, "mean")
    raise ValueError(f"unknown statistic {statistic!r}; "
                     "expected 'std', 'gradient' or 'tophat'")


def benchmark_cell_free_masks(phase: npt.NDArray, labels: npt.NDArray,
                              block: int = DEFAULT_BLOCK,
                              keeps: Sequence[float] = (0.10, 0.25, 0.40),
                              dilation: int = 21) -> list[dict]:
    """How well do the phase detectors find the cell-free blocks?

    Scores each `cell_free_mask_from_phase` statistic against the segmentation,
    with picking blocks at random as the floor. Run this once per new cell line
    or magnification; if `precision` is not far above the `chance` row, use
    `cell_free_mask_from_labels`, or an estimator that needs no mask at all.

    Returns one row per (statistic, keep) plus a `chance` row, each with the
    fraction of kept blocks that hold no cell (`precision`) and their mean cell
    occupancy (`occupancy`, lower is better).
    """
    free = cell_free_mask_from_labels(labels, phase.shape, dilation)
    occ = block_reduce_robust((~free).astype(np.float32), block, "mean")
    rows = [{"statistic": "chance", "keep": 1.0,
             "precision": float((occ == 0).mean()),
             "occupancy": float(occ.mean())}]
    a = np.asarray(phase, np.float32)
    for statistic in ("std", "gradient", "tophat"):
        stat = _phase_texture(a, statistic, block)
        for k in keeps:
            sel = stat <= np.quantile(stat, k)
            rows.append({"statistic": statistic, "keep": float(k),
                         "precision": float((occ[sel] == 0).mean()),
                         "occupancy": float(occ[sel].mean())})
    return rows


# --------------------------------------------------------------------------
# background
# --------------------------------------------------------------------------

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
        From `cell_free_mask_from_labels` or `cell_free_mask_from_phase`. With
        a mask the estimator defaults to a per-block median; without one it
        defaults to `clipped_mean`, which rejects cells on its own.
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


def temporal_background(stack: npt.NDArray, quantile: float = 0.10,
                        normalise_level: bool = True) -> npt.NDArray:
    """Background from a low quantile through time at each pixel.

    Needs no cell-free region in any single frame - only that each pixel is
    cell-free at *some* point in the movie, which moving cells provide. The
    natural choice when confluence is too high for `background_surface` to find
    usable blocks.

    Two caveats, both real on this data. Cells that never move leave their
    footprint in the result, so check the output for cell-shaped structure.
    And the background level drifts through the movie (11-21% here), so a
    quantile across raw frames mixes times; `normalise_level` divides each
    frame by its own median first and restores the median level at the end,
    which removes that.

    Parameters
    ----------
    stack : (t, y, x) array
        Frames, ideally the whole movie. Subsample it if memory is tight -
        30 frames is plenty.
    quantile : float
        Low quantile through time. Lower rejects cells harder and sits deeper
        in the noise; 0.05-0.20 is the useful range.
    """
    a = np.asarray(stack, np.float32)
    if normalise_level:
        lv = np.median(a, axis=(1, 2), keepdims=True)
        a = a / lv * float(np.median(lv))
    return np.quantile(a, quantile, axis=0).astype(np.float32)


# --------------------------------------------------------------------------
# flat field
# --------------------------------------------------------------------------

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
            f"more, or fall back to flatfield_from_signal.")
    if contrast < 0:                     # caller swapped them; harmless
        diff = -diff
    return _normalise_flatfield(diff, smooth)


def _mean_blank_surface(blank: npt.NDArray, block: int) -> npt.NDArray:
    """Mean block-reduced surface of a cell-free stack."""
    a = np.asarray(blank, np.float32)
    frames = a if a.ndim == 3 else a[None]
    return np.stack([block_reduce_robust(f, block, "clipped_mean")
                     for f in frames]).mean(axis=0)


def flatfield_from_signal(x: npt.NDArray, y: npt.NDArray, signal: npt.NDArray,
                          shape: tuple[int, int],
                          grid: int = 8, smooth: float = 1.0,
                          min_per_bin: int = 50,
                          min_snr: float = 3.0) -> npt.NDArray:
    """Flat field from the cells themselves. For when there is no blank well.

    Bins background-subtracted cell measurements by where in the field they
    were taken and takes the median per bin. The assumption is that cells are
    on average equally bright everywhere in the field, so any remaining trend
    across the field is the illumination.

    It needs no blank well, no cell-free region, and - because cell signal
    above background carries no additive term - no camera offset. Pool the
    measurements over every position on the plate: the flat field is a property
    of the optics, so all of them measure the same thing.

    **It needs contrast, and that is what usually rules it out.** On the
    20250213 plate the cells sit ~10% above their own local background (mean
    140 counts inside the mask against 127 in the gaps), which leaves ~12
    counts of signal per cell; the 8x8 bin medians then scatter over
    0.30-1.75 and the "flat field" is noise. `min_snr` catches that and
    raises. Prefer `flatfield_from_blank_pair` whenever two blanks exist.

    Two further things to watch, when contrast is not the problem:

    * It is circular if you then judge it by how flat the corrected signal is.
      Fit it on some positions and measure the flatness on the others.
    * It assumes the cell population is spatially homogeneous. Anything that
      makes cells at the edge of a field genuinely different - drift out of
      focus, a density gradient from seeding - lands in the flat field.

    Parameters
    ----------
    x, y : arrays
        Cell centroids in **fluorescence-frame** pixels. cellaap tables store
        them at segmentation scale, so double them for a 2x-upsampled frame.
    signal : array
        Per-cell signal with the background already subtracted - e.g. from
        `PositionCorrection.subtract_background`. Not divided by anything.
    shape : (rows, cols)
        Fluorescence frame shape, so `x`/`y` can be binned.
    grid : int
        Bins per axis. Coarse on purpose: 8 x 8 over a plate's worth of
        measurements is already thousands of cells per bin, and the
        illumination has no fine structure to resolve.
    min_per_bin : int
        Bins with fewer measurements are filled from their neighbours.
    min_snr : float
        Refuse if the fitted field is not smooth enough to be illumination:
        the spread between neighbouring bins has to be this many times smaller
        than the spread across the whole field. Set to 0 to fit anyway.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    s = np.asarray(signal, float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(s)
    x, y, s = x[ok], y[ok], s[ok]
    if x.size == 0:
        raise ValueError("no finite measurements to fit a flat field to")

    ix = np.clip((x / shape[1] * grid).astype(int), 0, grid - 1)
    iy = np.clip((y / shape[0] * grid).astype(int), 0, grid - 1)
    flat_index = iy * grid + ix
    out = np.full(grid * grid, np.nan, np.float32)
    order = np.argsort(flat_index, kind="stable")
    fi, sv = flat_index[order], s[order]
    edges = np.searchsorted(fi, np.arange(grid * grid + 1))
    for k in range(grid * grid):
        lo, hi = edges[k], edges[k + 1]
        if hi - lo >= min_per_bin:
            out[k] = np.median(sv[lo:hi])
    field = out.reshape(grid, grid)

    if min_snr > 0:
        # Illumination is smooth, so neighbouring bins should differ far less
        # than the field spans. When they do not, the bins are measuring noise.
        g = _fill_and_smooth(field, 0.0)
        rough = np.mean([np.std(np.diff(g, axis=0)), np.std(np.diff(g, axis=1))])
        span = float(np.ptp(g))
        snr = span / rough if rough > 0 else np.inf
        if snr < min_snr:
            raise ValueError(
                f"the cell signal is too faint to fit a flat field to: bin-to-bin "
                f"scatter is {rough:.3g} against a {span:.3g} span across the field "
                f"(SNR {snr:.1f} < {min_snr}). Use flatfield_from_blank_pair, or "
                f"correct the background only.")
    return _normalise_flatfield(field, smooth)


def flatfield_from_surfaces(surfaces: npt.NDArray, darkfield: float,
                            smooth: float = 1.0) -> npt.NDArray:
    """Flat field from background surfaces, once the camera offset is known.

    `mean(surfaces) - darkfield` is `A_bar F`; normalising to mean 1 gives `F`.
    This is exactly the step the blank-well intensity map skips: with
    `darkfield=0` you get back the too-flat map the existing pipeline divides
    by.

    **Only use this with a `darkfield` from a measured dark frame.** The
    background here is ~130 counts, so `D = 20` amplifies the measured falloff
    by 130/110 = 1.2x while `D = 100` amplifies it by 130/30 = 4.3x - and the
    estimators in this module bracket `D` no better than +-30 counts on this
    plate. Prefer `flatfield_from_signal`, which does not have this problem.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean_surface = np.nanmean(np.asarray(surfaces, np.float32), axis=0)
    return _normalise_flatfield(mean_surface - float(darkfield), smooth)


def flatfield_from_blank(blank: npt.NDArray, darkfield: float,
                         block: int = DEFAULT_BLOCK,
                         smooth: float = 1.0) -> npt.NDArray:
    """Flat field from a blank-well stack, with the camera offset removed.

    The blank-well route done properly. The shape is an optical property, so a
    blank well is a legitimate source for it - unlike the *level*, which is not
    transferable between wells. Use the brightest blank available: the higher
    `A / D`, the less the offset subtraction has to do and the less its error
    matters. Carries the same `darkfield` warning as `flatfield_from_surfaces`.
    """
    a = np.asarray(blank, np.float32)
    frames = a if a.ndim == 3 else a[None]
    grids = np.stack([block_reduce_robust(f, block, "clipped_mean") for f in frames])
    return flatfield_from_surfaces(grids, darkfield, smooth)


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


# --------------------------------------------------------------------------
# camera offset
# --------------------------------------------------------------------------

def darkfield_from_flatfield(surface: npt.NDArray,
                             flatfield: npt.NDArray) -> dict:
    """Camera offset from a measured surface and a known flat field.

    A cell-free surface is `D + A F`, so with `F` already in hand from
    `flatfield_from_blank_pair` the offset is just the intercept of a straight
    line through `surface` against `F`. Unlike `level_regression`, which
    extrapolates to zero medium brightness from a 20% lever arm, this fits
    inside the data: `F` spans ~40% across the field, so the intercept is
    pinned by real measurements.

    Nothing in this module divides by the result - the correction subtracts the
    measured surface and divides by `F`, and neither step needs `D`. It is
    reported because a sane, stable value is good evidence the model holds, and
    because a `D` that comes out negative or above the dimmest blank means the
    flat field is wrong.

    On the 20250213 plate: 53 counts (GFP, from the DMEM blank) and 61 (from
    FluoroBrite); 76 and 79 in Texas Red. Agreement between two very different
    wells is the check worth reading.

    Returns
    -------
    dict with `darkfield`, `medium_amplitude` (the `A` of that surface), and
    the fit's `r_squared`.
    """
    s = np.asarray(surface, np.float32)
    f = np.asarray(flatfield, np.float32)
    if f.shape != s.shape:
        f = ndi.zoom(f, np.array(s.shape) / np.array(f.shape), order=1,
                     mode="nearest")
    ok = np.isfinite(s) & np.isfinite(f)
    a, b = np.polyfit(f[ok].ravel(), s[ok].ravel(), 1)
    resid = s[ok] - (a * f[ok] + b)
    var = float(np.var(s[ok]))
    return {"darkfield": float(b), "medium_amplitude": float(a),
            "r_squared": float(1 - np.var(resid) / var) if var > 0 else float("nan")}


def darkfield_from_ptc(stack: npt.NDArray, window: int = 5,
                       n_bins: int = 40, max_drift_percent: float = 5.0) -> dict:
    """Camera offset from the photon-transfer relation of a blank stack.

    For a shot-noise-limited sensor the temporal variance of a pixel is
    `var = g (mu - D) + read^2`, so variance against mean is a line whose
    x-intercept is the offset. Temporal variance over a short window is used
    rather than spatial variance, because spatial variance would pick up the
    flat field and the sensor's fixed pattern instead of the shot noise.

    Windows are fitted separately and combined by median, with any window whose
    fitted gain is non-positive or far from the median gain thrown out: on the
    GFP DMEM blank, which bleaches ~70% over the movie, individual windows
    otherwise return gains of -0.4 and offsets in the thousands.

    Even fitted carefully this is an extrapolation below the measured range.
    On the 20250213 plate it lands between -35 and +25 counts across wells and
    channels - enough to say the offset is small, not enough to divide by.
    Prefer a real dark frame (shutter closed, same exposure and gain), and
    prefer `flatfield_from_signal`, which needs no offset at all.

    Parameters
    ----------
    stack : (t, y, x) array
        A **cell-free** stack: a blank well, or a movie of an empty field.
        Cells break the relation - they change a pixel's mean and its variance
        together, off the line.
    window : int
        Consecutive frames per variance estimate. Short enough that the level
        does not drift much inside it; residual drift is divided out, and
        windows drifting more than `max_drift_percent` are skipped.

    Returns
    -------
    dict with `darkfield`, `gain` (counts per electron), `read_variance`, the
    per-window spread, and how many windows survived.
    """
    a = np.asarray(stack, np.float32)
    if a.shape[0] < window + 1:
        raise ValueError(f"need at least {window + 1} frames for a PTC fit")

    fits = []
    for s0 in range(0, a.shape[0] - window + 1, max(1, window * 3)):
        w = a[s0:s0 + window]
        lv = w.mean(axis=(1, 2))
        if 100 * (lv.max() - lv.min()) / lv.mean() > max_drift_percent:
            continue
        w = w * (lv[0] / lv)[:, None, None]        # take the slow drift out
        mu = w.mean(0).ravel()
        var = w.var(0, ddof=1).ravel()

        lo, hi = np.percentile(mu, [0.5, 99.5])
        edges = np.linspace(lo, hi, n_bins + 1)
        idx = np.digitize(mu, edges) - 1
        bm, bv = [], []
        for k in range(n_bins):
            sel = idx == k
            if sel.sum() > 500:
                bm.append(float(mu[sel].mean()))
                # the mean of a variance estimate is unbiased but a hot pixel
                # or a cosmic ray moves it, so trim the top 2%
                v = np.sort(var[sel])
                bv.append(float(v[: max(1, int(0.98 * v.size))].mean()))
        if len(bm) >= 5:
            gain, c = np.polyfit(bm, bv, 1)
            fits.append((float(gain), float(c), min(bm), max(bm)))

    if not fits:
        raise ValueError("no usable windows: every one drifted more than "
                         f"{max_drift_percent}% or was too sparse")

    gains = np.array([f[0] for f in fits])
    med_gain = float(np.median(gains))
    good = [f for f in fits
            if f[0] > 0 and abs(f[0] - med_gain) <= 0.5 * abs(med_gain)]
    if not good:
        raise ValueError("no window gave a physical gain; this stack is not "
                         "shot-noise limited (bleaching? cells? averaging?)")

    darks = np.array([-c / g for g, c, _, _ in good])
    return {"darkfield": float(np.median(darks)),
            "darkfield_spread": float(np.percentile(darks, 84)
                                      - np.percentile(darks, 16)),
            "gain": float(np.median([g for g, _, _, _ in good])),
            "read_variance": float(np.median([c for _, c, _, _ in good])
                                   + np.median([g for g, _, _, _ in good])
                                   * float(np.median(darks))),
            "mean_range": (min(f[2] for f in good), max(f[3] for f in good)),
            "n_windows": len(good), "n_windows_tried": len(fits)}


def darkfield_from_shape_consistency(background_shape: npt.NDArray,
                                     signal_shape: npt.NDArray,
                                     background_mean: float) -> float:
    """Camera offset from the gap between two profiles that should match.

    The cells report the flat field directly (`S F`), while the cell-free
    background reports `D + A F` - the same shape, flattened towards uniform by
    the offset. The amount of flattening fixes `D`::

        (B_centre - D) / (B_edge - D) = F_centre / F_edge

    Both inputs are radial profiles or grids normalised to mean 1;
    `background_mean` restores the background to counts. Uses only data already
    in hand - no blank well, no dark frame.

    On the 20250213 plate this gives ~34 counts (GFP) and ~52 (Texas Red),
    the same order as `darkfield_from_ptc`. Treat it as a cross-check, not a
    calibration: it inherits every assumption `flatfield_from_signal` makes
    about cells being equally bright everywhere.
    """
    b = centre_edge_ratio(np.asarray(background_shape, float))
    s = centre_edge_ratio(np.asarray(signal_shape, float))
    if not np.isfinite(b) or not np.isfinite(s) or abs(s - b) < 1e-9:
        return float("nan")
    # B = D + A F, normalised to mean 1 means A = background_mean - D, and
    # solving (D + a f_c)/(D + a f_e) = s for D with b as the measured ratio
    f_c, f_e = s, 1.0
    frac = (b - 1) / (f_c - f_e) if abs(f_c - f_e) > 1e-9 else np.nan
    return float(background_mean * (1 - frac))


def level_regression(surfaces: npt.NDArray) -> tuple[npt.NDArray, npt.NDArray, float]:
    """Split a set of background surfaces into a shape and an offset.

    Each surface is `D + A_k F`, so regressing every grid point on the surface
    means `M_k` gives

        slope     = F / mean(F)      the flat field, mean-normalised
        intercept = D (1 - slope)    if D is spatially uniform

    and the slope of intercept against `1 - slope` is `D`. Differences in level
    between surfaces are what drives it, so it needs a real lever arm: a blank
    well whose brightness drifts through the movie is ideal.

    Returns `(flatfield, intercept, darkfield)`.

    Kept for the blank-well case, where it works: on the DMEM blank, whose
    level moves 70%, it recovers the flat field to a model residual of 0.1% of
    the mean. **It does not work on sample positions**, where the lever arm is
    only 8-21% and rising background comes partly from out-of-focus cell haze,
    which is flatter than `F`: that biases the slope towards 1 and inflates the
    extrapolated darkfield, giving 106 counts (GFP) and 143 (Texas Red) on this
    plate - the latter above the FluoroBrite blank's own level, so impossible.
    """
    S = np.asarray(surfaces, np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        M = np.nanmean(S, axis=(1, 2))
    X = M - M.mean()
    if (X ** 2).sum() <= 0:
        raise ValueError("all surfaces are at the same level; no lever arm to "
                         "regress against")

    W = np.isfinite(S)
    Sf = np.where(W, S, 0.0)
    n = W.sum(0)
    sx = (W * X[:, None, None]).sum(0)
    sy = Sf.sum(0)
    sxx = (W * (X ** 2)[:, None, None]).sum(0)
    sxy = (Sf * X[:, None, None]).sum(0)
    den = n * sxx - sx ** 2
    ok = den > 0
    slope = np.where(ok, (n * sxy - sx * sy) / np.where(ok, den, 1.0), np.nan)
    mean_s = np.where(n > 0, sy / np.maximum(n, 1), np.nan)
    intercept = mean_s - slope * M.mean()

    good = np.isfinite(slope) & np.isfinite(intercept) & (np.abs(1 - slope) > 0.02)
    darkfield = (float(np.polyfit((1 - slope)[good], intercept[good], 1)[0])
                 if good.sum() > 100 else float("nan"))
    return slope, intercept, darkfield


# --------------------------------------------------------------------------
# the per-position result
# --------------------------------------------------------------------------

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

        This is what `flatfield_from_signal` wants as its `signal`.
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
