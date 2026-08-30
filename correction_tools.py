"""Reading, inspecting and re-applying corrections after the pipeline has run.

`signal_correction` estimates a correction while a position is being analyzed.
This module is the other half: everything you might want to do with one
*afterwards* - look at the background that was subtracted, apply a correction to
numbers already sitting in a workbook, or repair a position whose own
background could not be measured.

It is separate so that neither concern clutters the other. Nothing in
`signal_correction` imports this; this imports `signal_correction`.

What the pipeline leaves on disk
--------------------------------
::

    <root>/pipeline/state/flatfield/<channel>.npz          one per plate
    <root>/pipeline/state/surfaces/<stem>_<channel>_bkg.tif  one per position

Those two files are together enough to reconstruct the whole correction, which
is what `Surface.from_files` does - no images are re-read and nothing is
re-estimated.

Usage and worked examples are in CORRECTION_TOOLS_README.md.
"""

from __future__ import annotations

import json
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import numpy.typing as npt

from signal_correction import PositionCorrection, upsample

__all__ = [
    # the saved background surface
    "save_background_stack",
    "read_background_stack",
    "upsample_stack",
    "surface_diagnostics",
    # applying a correction to numbers already measured
    "Surface",
    "correct_cell_table",
    "track_means",
    # repairing a position that could not measure its own background
    "average_shape",
    "borrow_shape",
]

# cellaap_analysis._load_maps walks the entire root - the pipeline directory
# included - and treats any file whose name matches this as a plate-wide
# correction map, imreading whatever it finds. A per-position surface caught
# that way would be applied to every position on the plate, silently
# reinstating the blank-well bug signal_correction exists to remove.
_RESERVED_IN_FILENAMES = re.compile(r"background|intensity")


# --------------------------------------------------------------------------
# the saved background surface
# --------------------------------------------------------------------------

def save_background_stack(correction: PositionCorrection, path: str | Path,
                          compression: str | None = "zlib") -> Path:
    """Write the per-frame background map as a TIFF stack, for inspection.

    `(n_frames, gy, gx)` float32 **in counts**, at the correction's own grid
    resolution - 32 x 32 for the default 64-pixel block on a 2048 frame. It is
    not upsampled on the way out because the surface genuinely has no structure
    finer than one block; `read_background_stack(..., shape=...)` expands it
    when something needs frame-sized pixels.

    About 490 kB per position per channel compressed (137 frames), so ~12 MB
    for a 13-position two-channel plate - 0.03% of that plate's raw stacks.

    **The file name must not contain "background" or "intensity"**, for the
    reason given at the top of this module. This refuses such a name rather
    than trusting the caller to remember.
    """
    path = Path(path).with_suffix(".tif")
    if _RESERVED_IN_FILENAMES.search(path.name):
        raise ValueError(
            f"{path.name!r} contains 'background' or 'intensity'; "
            f"cellaap_analysis._load_maps would pick it up as a plate-wide "
            f"correction map and apply this one position's surface to every "
            f"position. Name it '<stem>_<channel>_bkg.tif' or similar.")
    path.parent.mkdir(parents=True, exist_ok=True)

    stack = correction.background_stack()
    meta = {"stem": correction.stem, "channel": correction.channel,
            "units": "counts", "block": correction.block,
            "frame_shape": list(correction.frame_shape),
            "grid_shape": list(correction.background_shape.shape),
            "n_frames": int(stack.shape[0]),
            "what": "per-frame background B(x,y,t) that was subtracted; "
                    "upsample to frame_shape to overlay on the raw stack",
            "diagnostics": correction.diagnostics}
    import tifffile
    tifffile.imwrite(path, stack, compression=compression,
                     description=json.dumps(meta))
    return path


def read_background_stack(path: str | Path,
                          shape: tuple[int, int] | None = None
                          ) -> tuple[npt.NDArray, dict]:
    """Read a stack written by `save_background_stack`.

    Returns `(stack, metadata)`. Pass `shape` - or `metadata["frame_shape"]` -
    to get it back at frame resolution, ready to subtract from or overlay on
    the raw images.
    """
    import tifffile
    with tifffile.TiffFile(path) as fh:
        stack = fh.series[0].asarray().astype(np.float32)
        try:
            meta = json.loads(fh.pages[0].description)
        except Exception:
            meta = {}
    if shape is not None:
        stack = upsample_stack(stack, shape)
    return stack, meta


def upsample_stack(stack: npt.NDArray, shape: tuple[int, int]) -> npt.NDArray:
    """`upsample` every frame of an `(n, gy, gx)` stack to `shape`.

    Bilinear, frame by frame, so a 137 x 32 x 32 stack expanded to 2048 x 2048
    becomes 2.3 GB - materialise a slice rather than the whole movie unless you
    mean it.
    """
    a = np.asarray(stack, np.float32)
    if a.ndim == 2:
        return upsample(a, shape)
    return np.stack([upsample(f, shape) for f in a])


def surface_diagnostics(root: str | Path, channel: str | None = None,
                        subdir: str = "pipeline/state/surfaces"):
    """What the background estimator had to work with, position by position.

    Reads only the TIFF headers of the saved surfaces, so a whole plate costs
    a fraction of a second and no image data is touched.

    This is the companion to `cellaap_aggregate.baseline_offsets`: that measures
    that a position's zero is off, this says why. The columns to read first are

    ``usable_blocks``   fraction of the grid that had enough cell-free pixels
                        to measure. It falls as the field fills up.
    ``dilation_used``   how far from the cells the estimator managed to stay,
                        widest and narrowest over the sampled frames. It backs
                        off the requested value when a frame is too crowded,
                        and every step down lets more out-of-focus halo into
                        the background. That is the mechanism: the background
                        comes out too high, so the corrected signal - and the
                        position's zero with it - comes out too LOW.
    ``drift_percent``   how much the mean background moved over the movie. The
                        medium does drift, but tens of percent on a plate whose
                        uncrowded positions drift 5% is cells accumulating, not
                        medium.

    On the 20260826 CycB plate, over the ten positions of the two wells with no
    GFP induced, the three sort together exactly as that story predicts. The
    D01 positions hold 121 px and drift 5-6%; the A01 positions fall to 60 or
    30 px and drift 23-34%, and their floors are the low ones. Against
    `cellaap_aggregate.baseline_offsets`' floor: r = -0.83 for the peak
    background, -0.80 for the drift, +0.79 for the narrowest dilation held.

    Read it only over positions expected to hold the same fluorophore. A well
    that is genuinely brighter has a genuinely higher floor, and this table
    cannot tell you which you are looking at - it tells you whether the
    estimator was in trouble.

    Parameters
    ----------
    root : path
        The plate's root folder - the one holding `pipeline/`.
    channel : str, optional
        Only this channel's surfaces, e.g. `"GFP"`.
    subdir : str
        Where the pipeline put them, if it is not the default.

    Returns
    -------
    DataFrame, one row per position and channel. Empty if the folder holds no
    surfaces - a plate analyzed before they were written, or one whose
    `pipeline/` directory was not copied along with the results.
    """
    import pandas as pd
    import tifffile

    folder = Path(root) / subdir
    pattern = f"*_{channel}_bkg.tif" if channel else "*_bkg.tif"
    rows = []
    for path in sorted(folder.glob(pattern)):
        with tifffile.TiffFile(path) as fh:
            try:
                meta = json.loads(fh.pages[0].description)
            except Exception:
                warnings.warn(f"{path.name} carries no readable metadata; skipped")
                continue
        diagnostics = meta.get("diagnostics", {})
        stub = re.search(r"([A-H]\d{2})_s(\d+)", meta.get("stem", path.name))
        dilation = diagnostics.get("dilation_used", [None, None])
        rows.append({
            "stem": meta.get("stem", ""),
            "well": stub.group(1) if stub else "",
            "position": f"s{int(stub.group(2))}" if stub else "",
            "channel": meta.get("channel", ""),
            "model": diagnostics.get("background_model", ""),
            "usable_blocks": diagnostics.get("usable_block_fraction", float("nan")),
            "dilation_min": min(dilation) if dilation else None,
            "dilation_max": max(dilation) if dilation else None,
            "drift_percent": diagnostics.get("background_drift_percent", float("nan")),
            "background_min": (diagnostics.get("background_mean_range") or [None, None])[0],
            "background_max": (diagnostics.get("background_mean_range") or [None, None])[-1],
            "shape_centre_edge": diagnostics.get("background_shape_centre_edge",
                                                 float("nan")),
            "borrowed_from": ", ".join(diagnostics.get("shape_borrowed_from", [])),
            "reason": diagnostics.get("background_model_reason", ""),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# applying a correction to numbers already measured
# --------------------------------------------------------------------------

@dataclass
class Surface:
    """A correction rebuilt from what the pipeline left on disk.

    Deliberately thinner than `PositionCorrection`: it holds the background
    exactly as it was applied - one grid per frame, no shape/level split - plus
    the flat field, which is all that re-deriving a corrected number needs.
    """

    background: npt.NDArray          # (n_frames, gy, gx), counts
    flatfield: npt.NDArray | None    # (gy, gx), mean 1
    frame_shape: tuple[int, int]
    meta: dict

    @classmethod
    def from_files(cls, surface_tif: str | Path,
                   flatfield_npz: str | Path | None = None) -> "Surface":
        """Rebuild from `<stem>_<channel>_bkg.tif` and the plate's flat field.

        Pass `flatfield_npz=None` to correct background only, which is what the
        position got if the plate declared no blank pair.
        """
        stack, meta = read_background_stack(surface_tif)
        flat = None
        if flatfield_npz is not None:
            with np.load(Path(flatfield_npz), allow_pickle=False) as z:
                flat = z["flatfield"].astype(np.float32)
        return cls(background=stack, flatfield=flat,
                   frame_shape=tuple(meta.get("frame_shape", (0, 0))), meta=meta)

    @classmethod
    def from_correction(cls, correction: PositionCorrection) -> "Surface":
        """Same thing straight from a live `PositionCorrection`, no file needed."""
        return cls(background=correction.background_stack(),
                   flatfield=correction.flatfield,
                   frame_shape=correction.frame_shape,
                   meta={"stem": correction.stem, "channel": correction.channel,
                         "diagnostics": dict(correction.diagnostics)})

    def _index(self, x, y, frame):
        gy, gx = self.background.shape[1:]
        j = np.clip((np.asarray(x, float) / self.frame_shape[1] * gx).astype(int), 0, gx - 1)
        i = np.clip((np.asarray(y, float) / self.frame_shape[0] * gy).astype(int), 0, gy - 1)
        f = np.clip(np.asarray(frame, int), 0, self.background.shape[0] - 1)
        return f, i, j

    def background_at(self, x, y, frame) -> npt.NDArray:
        """The background that was subtracted under these cells, in counts."""
        f, i, j = self._index(x, y, frame)
        return self.background[f, i, j]

    def correct(self, raw, x, y, frame) -> npt.NDArray:
        """`(raw - B) / F` for per-cell measurements.

        `x` and `y` must be in **fluorescence-frame** pixels; cellaap writes
        centroids at segmentation scale, so double them for a 2x-upsampled
        frame. `correct_cell_table` does that for you.
        """
        f, i, j = self._index(x, y, frame)
        out = np.asarray(raw, float) - self.background[f, i, j]
        if self.flatfield is not None:
            fl = self.flatfield
            if fl.shape != self.background.shape[1:]:
                fl = upsample(fl, self.background.shape[1:])
                fl = fl / fl.mean()
            out = out / fl[i, j]
        return out


def correct_cell_table(cells, surfaces: Mapping[str, Surface],
                       xy_scale: float | None = None,
                       x_col: str = "x", y_col: str = "y",
                       frame_col: str = "frame", suffix: str = "_corrected"):
    """Add `<channel><suffix>` columns to a `cell_data` table.

    Parameters
    ----------
    cells : DataFrame
        The `cell_data` sheet of a `*_analysis.xlsx`, or the cached parquet.
    surfaces : {channel: Surface}
        One per channel to correct. Channels absent from `cells` are skipped.
    xy_scale : float, optional
        Segmentation-to-fluorescence pixel scale. Inferred from the frame shape
        and the centroid range when left out - cellaap writes centroids at half
        the fluorescence resolution, so this is 2 on a 2048 frame. It is
        checked either way, because getting it wrong samples the wrong part of
        the field and says nothing.
    """
    out = cells.copy()
    for channel, surface in surfaces.items():
        if channel not in out.columns:
            warnings.warn(f"{channel!r} is not a column of this table; skipped")
            continue
        scale = xy_scale if xy_scale is not None else _infer_xy_scale(out, surface,
                                                                     x_col, y_col)
        x = out[x_col].to_numpy(float) * scale
        y = out[y_col].to_numpy(float) * scale
        _check_in_frame(x, y, surface)
        out[f"{channel}{suffix}"] = surface.correct(
            out[channel].to_numpy(float), x, y, out[frame_col].to_numpy(int))
    return out


def _infer_xy_scale(cells, surface: Surface, x_col: str, y_col: str) -> float:
    """Nearest power-of-two scale that puts the centroids inside the frame."""
    span = max(float(cells[x_col].max()), float(cells[y_col].max()))
    if span <= 0:
        return 1.0
    scale = surface.frame_shape[1] / (2 ** np.ceil(np.log2(span + 1)))
    return float(max(scale, 1.0))


def _check_in_frame(x, y, surface: Surface) -> None:
    h, w = surface.frame_shape
    if x.max() > w or y.max() > h:
        raise ValueError(
            f"centroids reach ({x.max():.0f}, {y.max():.0f}) but the frame is "
            f"{w} x {h}; the xy_scale is wrong")
    if x.max() < 0.4 * w and y.max() < 0.4 * h:
        warnings.warn(
            f"centroids only reach ({x.max():.0f}, {y.max():.0f}) in a "
            f"{w} x {h} frame - is xy_scale too small? The correction would be "
            f"sampled from one corner of the field.")


def track_means(cells, summary, columns: Sequence[str],
                delta_t: float = 10.0, key_cols: Sequence[str] = ("particle",),
                start_col: str = "mitotic_start_frame",
                duration_col: str = "corrected_frames_in_mitosis",
                frame_col: str = "frame", collision_suffix: str = "_recomputed"):
    """Average per-frame columns over each track's mitotic window.

    Reproduces what `summarize_data` does, so a column added to `cell_data`
    after the fact can be carried to the summary without re-running the
    analysis. Validated against the summary's own raw column on the 20250213
    plate: r = 0.995, median difference -0.5%.

    `duration_col` is in **minutes** once `cellaap_aggregate.load_experiment`
    has been through the summary; `delta_t` converts it back to frames. Pass
    `delta_t=1` if your summary still holds frames.

    A requested column that the summary already carries - asking for `GFP` when
    the summary has its own `GFP` - comes back as `GFP<collision_suffix>` so the
    original is never silently replaced. Recomputing the raw column that way is
    the check worth running before trusting a corrected one: on the 20250213
    plate it reproduces the summary to r = 0.999 within -0.6%.
    """
    import pandas as pd

    # A raw summary from summarize_data holds frames; load_experiment converts
    # it to minutes. Getting this backwards shortens every window by delta_t
    # and still returns plausible numbers - r falls to 0.970 with a -4% bias
    # where the right value gives 0.9996 and +0.00% - so it is worth a look.
    longest = float(np.nanmax(summary[duration_col])) if len(summary) else 0.0
    if delta_t > 1 and longest and longest < 60:
        warnings.warn(
            f"{duration_col} tops out at {longest:g}, which looks like frames "
            f"rather than minutes, but delta_t={delta_t:g} assumes minutes. "
            f"Pass delta_t=1 for a summary straight from summarize_data.")

    win = summary.set_index(list(key_cols))[[start_col, duration_col]]
    joined = cells.join(win, on=list(key_cols), how="inner")
    inside = ((joined[frame_col] >= joined[start_col])
              & (joined[frame_col] < joined[start_col] + joined[duration_col] / delta_t))
    joined = joined[inside]
    agg = joined.groupby(list(key_cols)).agg(
        n_frames=(frame_col, "size"),
        **{c: (c, "mean") for c in columns if c in joined.columns}).reset_index()
    return summary.merge(agg, on=list(key_cols), how="left",
                         suffixes=("", collision_suffix))


# --------------------------------------------------------------------------
# repairing a position that could not measure its own background
# --------------------------------------------------------------------------

def average_shape(donors: Iterable[PositionCorrection]) -> npt.NDArray:
    """Mean background shape of several positions, normalised to mean 1."""
    shapes = [np.asarray(d.background_shape, np.float32) for d in donors]
    if not shapes:
        raise ValueError("no donors given")
    if len({s.shape for s in shapes}) != 1:
        raise ValueError(f"donors disagree on grid size: {[s.shape for s in shapes]}")
    m = np.mean(shapes, axis=0)
    return (m / m.mean()).astype(np.float32)


def borrow_shape(target: PositionCorrection,
                 donors: Iterable[PositionCorrection],
                 note: str = "") -> PositionCorrection:
    """Give a position a donor background **shape**, keeping its **own level**.

    For a position too crowded to measure a usable background of its own. Only
    the shape is borrowed, and that asymmetry is the whole point:

    * the *shape* transfers. Across the 20250213 plate each position's
      background shape differs from the plate mean by 1.6% (GFP) / 2.1% (Texas
      Red), so borrowing it costs ~2-3 counts.
    * the *level* does not. It varies 4.1 (GFP) / 7.3 (Texas Red) counts
      between positions, and - the part that matters - it **correlates with
      confluence** (r = +0.51, ~26 counts per unit confluence). A position that
      cannot measure its own background is by definition the most crowded one,
      so its true level is the highest on the plate while every donor is
      sparser and therefore darker. Borrowing the level would under-subtract by
      ~5 counts against a ~13-count signal, systematically, and in exactly the
      wells where cells grew densest.

    The level is therefore re-derived locally: the position's own measured
    block grids are projected onto the borrowed shape, which is the same
    least-squares step `estimate_position_correction` uses. Where those grids
    were not kept, the existing level is rescaled so the mean background is
    unchanged - weaker, and flagged in the diagnostics.

    Returns a new `PositionCorrection`; `target` is not modified. The result
    records the borrow under `diagnostics["shape_borrowed_from"]`, so a summary
    made from it says so.
    """
    import copy

    donors = list(donors)
    shape = average_shape(donors)
    if shape.shape != target.background_shape.shape:
        raise ValueError(
            f"donor grid {shape.shape} does not match this position's "
            f"{target.background_shape.shape}; they must share `block`")

    out = copy.deepcopy(target)
    out.background_shape = shape

    grids = target.background_grids
    if grids is not None and len(grids):
        w = float((shape ** 2).sum())
        level_at = np.array([float(((g - o) * shape).sum() / w)
                             for g, o in zip(grids, np.zeros(len(grids)))])
        idx = (target.sampled_frames if target.sampled_frames is not None
               else np.linspace(0, len(target.background_level) - 1, len(grids)).astype(int))
        out.background_level = np.interp(
            np.arange(len(target.background_level)), idx, level_at).astype(np.float32)
        how = "re-projected this position's own measured grids onto the borrowed shape"
    else:
        # No grids kept: preserve the mean background this position measured.
        before = (target.background_offset
                  + target.background_level * float(target.background_shape.mean()))
        out.background_level = (before / float(shape.mean())).astype(np.float32)
        out.background_offset = np.zeros_like(out.background_level)
        how = ("rescaled the existing level to preserve the mean background "
               "(background_grids were not kept, so it could not be re-projected)")

    out.diagnostics = dict(target.diagnostics)
    out.diagnostics.update({
        "shape_borrowed_from": [d.stem for d in donors],
        "shape_borrow_n_donors": len(donors),
        "level_kept_local": True,
        "level_method": how,
        "shape_borrow_note": note,
    })
    return out
