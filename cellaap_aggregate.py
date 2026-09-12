"""Compiling a plate's per-position summaries into per-condition tables.

`cellaap_analysis` turns one position into a `*_summary.xlsx`. This module is
the step after: it collects many of those into one DataFrame keyed by what was
in the well, and fits and plots the result.

The plate layout is not restated here. It is read from the same `platemap.csv`
that `pipeline.py` ran the plate from, through `pipeline`'s own parsing, so a
compiled group cannot disagree with what was actually segmented and analyzed.
That is the duplicated-truth problem the platemap exists to prevent, and it
applies to compilation exactly as it applies to the two pipeline stages.

Three things follow from reading the platemap rather than a hand-written well
list, and each fixes a way a hand-written list quietly got the wrong answer:

* **Well ranges and normalization.** `B01-B06` expands and `g3` becomes `G03`,
  through `pipeline.expand_well_token`. Splitting a `well_ids` cell on commas
  and whitespace alone leaves a range as the literal string `B01-B06`, which
  matches no folder at all.
* **`skip` and blank-media rows drop out.** A `fluorobrite` well has no
  celltype, transfection or drug, so it forms a junk `('', '', '')` group that
  then finds no summaries, because blank wells are never segmented.
* **Position overrides stop double-counting.** With a `G03` row and a `G03_s9`
  row under a different drug, matching the substring `_G03_` against folder
  names puts site 9 in *both* groups. Precedence is resolved per position
  here, as the pipeline resolves it, so each position lands in exactly one
  group.

Entry points, in decreasing order of how much you have to say:

    import cellaap_aggregate as agg

    # the platemap does everything
    df = agg.load_experiment(root, expt_length=150, delta_t=10)

    # or the pieces, to filter or regroup in between
    positions = agg.platemap_positions(root)
    groups    = agg.group_positions(positions)
    raw       = agg.compile_positions(groups[("HeLa", "pEN2", "DMSO")])

A compiled table still carries one defect: where a field is crowded,
`signal_correction` subtracts too much background, and cells with no
fluorophore read below zero. See "Putting each well's zero back where it
belongs" below - `well_offsets` measures one offset per well from its most
negative cells, `plot_well_offsets` shows it, and `apply_well_offsets`
subtracts it. There is deliberately no per-position version.

    # or the whole thing, which also writes `<signal>_zeroed` back into each
    # *_summary.xlsx
    df, offsets = agg.correct_wells(root)

The well-list API that predated `platemap_positions` - `create_wellmap_dict`,
`wellmap_from_platemap`, `compile_summaries`, `import_filter_data_for_wells`
and `import_whole_expt_data` - has been removed. `load_experiment` replaces the
whole chain and resolves per-position precedence, which that path could not.

This module imports matplotlib and seaborn at module scope. `cellaap_utils`
deliberately does not: it is imported by `cellaap_analysis` and therefore by
every analysis array task on the cluster, and none of them plot.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit

# pipeline.py imports only the standard library at module scope - it has to be
# runnable under both conda environments - so this is a cheap import, and it is
# the whole point: the platemap is parsed by the code that wrote it.
import pipeline

__all__ = [
    "PositionRef",
    "platemap_positions",
    "group_positions",
    "compile_positions",
    "filter_summary",
    "load_experiment",
    "read_platemap",
    "negative_tail_offset",
    "channel_columns",
    "well_offsets",
    "apply_well_offsets",
    "plot_well_offsets",
    "correct_wells",
    "export_to_excel_by_col",
    "fit_model",
    "sigmoid_4par",
]

# Columns `summarize_data` writes that the mitosis filter reads.
START_COL = "mitotic_start_frame"
DURATION_COL = "corrected_frames_in_mitosis"

# The data sheet in a *_summary.xlsx. It is written first, so an unqualified
# read_excel happens to land on it, but naming it means a later sheet reorder
# does not silently compile the parameters sheet instead.
SUMMARY_SHEET = "Summary"

# `<stem>_<model>_<confluency>_<threshold>_inference`, read from the right so
# that model names containing an underscore (HeLa_focal) do not eat the stem.
INFERENCE_DIR_RE = re.compile(
    r"_(?P<model>[^_]+(?:_focal)?)_(?P<confluency>\d+)_(?P<threshold>[0-9.]+)_inference$")


# ---------------------------------------------------------------------------
# The platemap, resolved down to positions
# ---------------------------------------------------------------------------

@dataclass
class PositionRef:
    """One sample position, and the condition the platemap puts it under.

    `inference_dir` is what the platemap *asks* for - the folder name encodes
    the model, confluency and threshold, so it is the folder the pipeline would
    have written. `summary_path()` is what is actually on disk.
    """
    stem: str               # 20251009_HeLa-siRNA_G03_s8_phs
    well: str               # G03
    site: int               # 8
    celltype: str
    transfection: str
    drug: str
    inference_dir: Path
    root: Path
    mapped: bool            # False when the well fell through to the default

    @property
    def stub(self) -> str:
        return f"{self.well}_s{self.site}"

    @property
    def position(self) -> str:
        """The site as the old `position` column spelled it: s8."""
        return f"s{self.site}"

    @property
    def key(self) -> tuple[str, str, str]:
        return (self.celltype, self.transfection, self.drug)

    @property
    def code(self) -> str:
        return f"{self.celltype}_{self.transfection}_{self.drug}"

    def summary_path(self, suffix: str = "") -> Path | None:
        """The summary spreadsheet for this position, or None if there is none.

        The platemap's folder is tried first. If it is not there - a plate
        re-compiled after the platemap's model or threshold was edited, say -
        any `<stem>_*_inference` folder is accepted instead, since for
        compilation the segmentation parameters do not change what a summary
        means. Two candidate folders is ambiguous and returns None rather than
        picking one.
        """
        candidates = [self.inference_dir] if self.inference_dir.is_dir() else [
            d for d in sorted(self.root.glob(f"{self.stem}_*_inference")) if d.is_dir()]
        if len(candidates) != 1:
            return None
        hits = sorted(candidates[0].glob(f"*_summary{suffix}.xlsx"))
        return hits[0] if len(hits) == 1 else (hits[0] if hits else None)


def _resolve_root(source) -> Path:
    """Accept a path, a string, or a `cellaap_analysis.analysis` session.

    The old API took the analysis object and reached into `.root_folder`;
    the platemap API only ever needs the folder, and most callers have one.
    """
    if hasattr(source, "root_folder"):
        return Path(source.root_folder)
    return Path(source)


def _positions_from_inference_dirs(root: Path) -> list[pipeline.Position]:
    """Fall back to the `*_inference` folder names when the phase stacks are gone.

    Compiling a plate whose raw stacks have been archived is a normal thing to
    want, and `pipeline.discover_positions` globs `*phs.tif`. The folder names
    carry the same stub, so the positions can be recovered from them.
    """
    seen: dict[str, pipeline.Position] = {}
    for directory in sorted(root.glob("*_inference")):
        if not directory.is_dir():
            continue
        trim = INFERENCE_DIR_RE.search(directory.name)
        stem = directory.name[:trim.start()] if trim else directory.name
        m = pipeline.STUB_RE.search(stem)
        if not m or stem in seen:
            continue
        seen[stem] = pipeline.Position(
            stem=stem, phs=root / f"{stem}.tif",
            well=pipeline.normalize_well(m.group("well")),
            site=int(m.group("site")))
    return list(seen.values())


def platemap_positions(source, pattern: str = "*phs.tif",
                       platemap: Path | str | None = None,
                       verbose: bool = True) -> list[PositionRef]:
    """Every sample position on the plate, tagged with its condition.

    Joins the positions in the folder to the platemap rows through
    `pipeline.build_tasks`, which is the same call the pipeline's own `check`,
    `submit` and `analyze` make. Precedence (exact position beats well beats
    the HeLa default), `skip` rows and the blank-media roles are therefore
    handled identically to the run itself.

    Parameters
    ----------
    source : Path | str | cellaap_analysis.analysis
        The root folder holding the stacks and inference folders, or an
        analysis session to take `.root_folder` from.
    pattern : str
        Glob for the phase stacks. If it matches nothing, the positions are
        recovered from the `*_inference` folder names instead.
    platemap : Path | str, optional
        A platemap elsewhere. Defaults to `<root>/platemap.csv`.
    verbose : bool
        Print the warnings `build_tasks` produces (unmapped wells, duplicate
        rows, platemap wells with no data).

    Returns
    -------
    list[PositionRef]
        One entry per sample position, in folder order. Blank-media wells and
        `skip` rows are not included.
    """
    root = _resolve_root(source)
    if not root.is_dir():
        raise FileNotFoundError(f"{root} is not a folder")

    path = Path(platemap) if platemap else pipeline.platemap_path(root)
    if not path.exists():
        raise FileNotFoundError(
            f"No platemap at {path}. Run `pipeline.py init --root {root}` first, "
            f"or pass platemap=... if it lives elsewhere.")
    rows = pipeline.read_platemap(path)

    positions = pipeline.discover_positions(root, pattern)
    if not positions:
        positions = _positions_from_inference_dirs(root)
        if positions and verbose:
            print(f"no {pattern} in {root}; took {len(positions)} position(s) "
                  f"from the inference folder names")
    if not positions:
        raise FileNotFoundError(
            f"no positions in {root}: neither {pattern} nor *_inference folders")

    # Channels are an inference/analysis concern; autodetecting them here would
    # be a filesystem glob per position for something no caller reads.
    tasks, _map_wells, warnings = pipeline.build_tasks(
        root, positions, rows, autodetect_channels=False)

    if verbose:
        for warning in warnings:
            print(f"  ! {warning}")

    return [PositionRef(stem=t.pos.stem, well=t.pos.well, site=t.pos.site,
                        celltype=t.celltype, transfection=t.transfection,
                        drug=t.drug, inference_dir=t.inference_dir,
                        root=root, mapped=t.mapped)
            for t in tasks]


def group_positions(positions: list[PositionRef],
                    by: tuple[str, ...] = ("celltype", "transfection", "drug"),
                    ) -> dict[tuple, list[PositionRef]]:
    """Bucket positions by condition, preserving the order they were found in.

    `by` names attributes of `PositionRef`, so grouping by a subset works too:
    `by=("drug",)` pools every cell type under each drug.
    """
    groups: dict[tuple, list[PositionRef]] = {}
    for pos in positions:
        key = tuple(getattr(pos, field) for field in by)
        groups.setdefault(key, []).append(pos)
    return groups


# ---------------------------------------------------------------------------
# Reading the summaries
# ---------------------------------------------------------------------------

def compile_positions(positions: list[PositionRef], suffix: str = "",
                      verbose: bool = True) -> pd.DataFrame:
    """Read and concatenate the summary spreadsheets for these positions.

    Each position's summary is read from its own inference folder rather than
    found by matching a well id against every folder name, so a site that the
    platemap moved to a different condition cannot also be counted under its
    well's condition.

    Parameters
    ----------
    positions : list[PositionRef]
        From `platemap_positions`, usually one group of `group_positions`.
    suffix : str
        Which summary variant to read: `""` for `*_summary.xlsx`, or e.g.
        `"_dead"` for the copies `augment_dead_label.py --suffix _dead` writes.
    verbose : bool
        Report each position loaded, and each one with no summary.

    Returns
    -------
    pandas.DataFrame
        The concatenated summaries with `well`, `position`, `stem`, `celltype`,
        `transfection`, `drug`, `code` and `storage_location` added. Empty (and
        with no columns) if none of the positions has a summary.
    """
    frames = []
    for pos in positions:
        path = pos.summary_path(suffix)
        if path is None:
            if verbose:
                print(f"  ! {pos.stub}: no *_summary{suffix}.xlsx "
                      f"(expected in {pos.inference_dir.name})")
            continue
        df = pd.read_excel(path, sheet_name=SUMMARY_SHEET)
        df["well"] = pos.well
        df["position"] = pos.position
        df["stem"] = pos.stem
        df["celltype"] = pos.celltype
        df["transfection"] = pos.transfection
        df["drug"] = pos.drug
        df["code"] = pos.code
        # As the well-list API defined it before this one - the folder the root
        # sits in, which identifies the experiment when several are pooled.
        df["storage_location"] = str(pos.root.parent)
        frames.append(df)
        if verbose:
            print(f"{pos.stub} loaded ({len(df)} tracks)")

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def filter_summary(df: pd.DataFrame, expt_length: int, delta_t: int | float,
                   ) -> pd.DataFrame:
    """Drop unobserved and unfinished mitoses, and put the duration in real time.

    Two filters, both reading `mitotic_start_frame` as an absolute movie frame,
    which is what `summarize_data` writes:

    * mitotic entry was actually observed (`mitotic_start_frame > 0`);
    * the whole episode finished inside the movie
      (`mitotic_start_frame + corrected_frames_in_mitosis < expt_length`).

    `corrected_frames_in_mitosis` is then multiplied by `delta_t`.

    Summary files written before `summarize_data` changed stored a within-track
    row index in `mitotic_start_frame`; against those the second filter is too
    permissive for tracks that begin late in the movie.
    """
    if df.empty:
        return df
    missing = {START_COL, DURATION_COL} - set(df.columns)
    if missing:
        raise KeyError(f"summary is missing column(s): {', '.join(sorted(missing))}")

    out = df[df[START_COL] > 0]
    out = out[out[START_COL] + out[DURATION_COL] < expt_length].copy()
    out[DURATION_COL] = out[DURATION_COL] * delta_t
    return out


def load_experiment(source, expt_length: int, delta_t: int | float,
                    pattern: str = "*phs.tif",
                    platemap: Path | str | None = None,
                    suffix: str = "", by: tuple[str, ...] = ("celltype", "transfection", "drug"),
                    filter_rows: bool = True, verbose: bool = True) -> pd.DataFrame:
    """A whole plate, compiled and filtered, from its platemap alone.

    The one call that compiles a plate:

        df = load_experiment(root, expt_length=150, delta_t=10)

    Its output still has each position on its own zero; see "Aligning the zero
    across positions" for putting them on a common one.

    Parameters
    ----------
    source : Path | str | cellaap_analysis.analysis
        Root folder, or an analysis session to take `.root_folder` from.
    expt_length : int
        Frames in the movie; mitoses that would end past it are dropped.
    delta_t : int | float
        Time per frame, multiplied into `corrected_frames_in_mitosis`.
    pattern, platemap, suffix, by, verbose
        As for `platemap_positions`, `compile_positions` and `group_positions`.
    filter_rows : bool
        Apply `filter_summary`. Set False for the unfiltered table.

    Returns
    -------
    pandas.DataFrame
        Every position's summary, concatenated, with the identity columns
        `compile_positions` adds. `code` is the group label, built from `by`.
    """
    positions = platemap_positions(source, pattern=pattern, platemap=platemap,
                                   verbose=verbose)
    groups = group_positions(positions, by=by)

    frames = []
    for key, members in groups.items():
        label = "_".join(str(k) for k in key)
        if verbose:
            print(f"\n{label}: {len(members)} position(s)")
        raw = compile_positions(members, suffix=suffix, verbose=verbose)
        if raw.empty:
            if verbose:
                print(f"  ! {label}: no summaries found, skipped")
            continue
        # `by` may be a subset of the three grouping columns, in which case the
        # per-position `code` is finer than the group; make them agree.
        raw["code"] = label
        frames.append(filter_summary(raw, expt_length, delta_t) if filter_rows else raw)

    if not frames:
        if verbose:
            print("\nNo summaries found anywhere on this plate.")
        return pd.DataFrame()

    whole = pd.concat(frames, ignore_index=True)
    if verbose:
        print(f"\n{len(whole)} rows from {whole['stem'].nunique()} position(s) "
              f"in {whole['code'].nunique()} group(s)")
    return whole


def read_platemap(source, platemap: Path | str | None = None) -> pd.DataFrame:
    """The platemap as a DataFrame.

    Thin wrapper over `pipeline.read_platemap_df` that also accepts a root
    folder rather than the file itself. The comment header needs
    `comment='#'`, which is what this exists to remember.
    """
    if platemap is not None:
        path = Path(platemap)
    else:
        path = Path(source)
        if path.is_dir():
            path = pipeline.platemap_path(path)
    return pipeline.read_platemap_df(path)


# ---------------------------------------------------------------------------
# Putting each well's zero back where it belongs
#
# `signal_correction` measures each position's background from that position's
# own cell-free pixels, and a crowded field does not have many. The estimator
# backs its exclusion ring off the cells to keep enough blocks measurable,
# counts more of their out-of-focus halo as medium, and subtracts too much. A
# cell with no fluorophore then reads BELOW zero - and those cells, the
# non-expressing and barely expressing ones, are the ones that define the
# bottom of a dose-response curve.
#
# The correction reads the over-subtraction off the negative values
# themselves. A well with more than `min_negative` cells below zero is shifted
# right by the median of its most negative `fraction` of those cells:
#
#     negatives = the well's values < 0, sorted
#     tail      = the lowest ceil(fraction * len(negatives)) of them
#     offset    = median(tail)              (a negative number)
#     zeroed    = value - offset
#
# A well with `min_negative` negative cells or fewer is left where it is: that
# is scatter around a zero that is already about right, and a median of a
# handful of cells is not a number to move a well by.
#
# The tail, and not all the negatives, is a choice and not an estimate. The
# negatives are the left half of the non-expressing peak, so shifting by their
# median leaves much of that peak below zero, and a Hill curve is only defined
# for x >= 0. A fit whose non-expressing cells sit below zero takes its base
# from where the rise starts rather than from the basal duration. On the
# pPS18 20250402 plate (Texas Red_corrected, these defaults) HeLa's fitted
# base went from 82 min to 59 and RPE1's from 65 to 43 - both onto the basal
# duration of their dimmest cells - with under 5% of cells left negative.
# What that costs is that the zero is a convention: the shift includes about
# one width of the non-expressing peak, and EC50 moves with it (HeLa
# 12.5 -> 25.2, RPE1 9.3 -> 12.4; U2OS, which already had a long flat foot,
# 40.2 -> 37.6 with its top poorly determined either way). A Hill fit
# with a free x-offset could not choose the zero - its fit is flat anywhere
# below the non-expressing peak - so EC50s compare only between data zeroed
# the same way, with the same `fraction`.
#
# One offset per well, never per position. On the HeLa well, per-position
# offsets were tested against the dose-response itself: one shared curve
# shape, a free shift per position, and the question whether the correction
# brings those shifts together. Offsets from each position's non-expressing
# peak widened their scatter (SD 2.8 -> 4.3 a.u.); per-position negative-tail
# offsets left it where it was (2.7). Neither tracked the shift a position's
# curve needed. A position holds a few hundred cells and its most negative
# quarter a few dozen, so a per-position offset adds noise to every cell and
# removes nothing. The sites of a well are pooled and moved together.
#
# Two things it cannot do:
#
#   * repair a position over-subtracted by different amounts across its own
#     field. A constant moves a distribution; a tail from one corner of the
#     field stays a tail.
#   * make two acquisitions comparable. A different exposure rescales the
#     signal as well as shifting it.
# ---------------------------------------------------------------------------

# The share of a well's negative cells the offset is read from - the most
# negative quarter, by default.
NEGATIVE_TAIL_FRACTION = 0.25

# A well is shifted only when MORE than this many of its cells are negative.
MIN_NEGATIVE_CELLS = 10

# Appended to a signal's name for the zeroed column.
ZEROED_SUFFIX = "_zeroed"

# The sheet each summary gets recording what was subtracted from it.
OFFSET_SHEET = "offset_correction"

# What the per-position alignment this replaced wrote into the summaries. A
# re-run removes them, so a file never carries two differently-zeroed columns.
LEGACY_SUFFIX = "_well_aligned"
LEGACY_SHEET = "well_alignment"

# Columns that name a position rather than a well. They are refused as keys:
# a per-position offset is the thing this section exists not to apply.
_POSITION_COLUMNS = ("position", "stem", "site")

# The channels `cellaap_analysis.summarize_data` knows how to measure.
CHANNELS = ("GFP", "Texas Red", "Cy5")

# Per-channel signal columns, best first. `<ch>_corrected` is
# signal_correction's number; `<ch>` is raw. `_std` columns describe scatter
# within a track and are not shifted by an offset, so they never appear here.
SIGNAL_SUFFIXES = ("_corrected", "_bkg_corr", "")


def channel_columns(df: pd.DataFrame, prefer: tuple[str, ...] = SIGNAL_SUFFIXES,
                    ) -> list[str]:
    """The best available signal column for each channel this table carries.

    One column per channel, not all of them: zeroing `GFP` and `GFP_corrected`
    separately would put two differently-zeroed numbers in one file under names
    that look like variants of each other. A column that is constant - all
    zeros for `_bkg_corr` when no legacy correction map was found - carries no
    signal and is skipped.
    """
    picked = []
    for channel in CHANNELS:
        for suffix in prefer:
            column = f"{channel}{suffix}"
            if column in df.columns and df[column].nunique(dropna=True) > 1:
                picked.append(column)
                break
    return picked


def _check_fraction(fraction: float) -> float:
    fraction = float(fraction)
    if not 0 < fraction <= 1:
        raise ValueError(f"fraction must be in (0, 1], got {fraction}")
    return fraction


def negative_tail_offset(values, fraction: float = NEGATIVE_TAIL_FRACTION,
                         min_negative: int = MIN_NEGATIVE_CELLS) -> float:
    """The offset one well's cells are shifted by. See the section notes above.

    Parameters
    ----------
    values : array-like
        One well's per-cell signal. NaNs are dropped.
    fraction : float
        Share of the negative values, taken from the most negative end, whose
        median is the offset. 1.0 is the median of every negative value.
    min_negative : int
        The offset is 0 unless more than this many values are negative.

    Returns
    -------
    float
        The median of the lowest `ceil(fraction * n_negative)` negative values
        - a negative number, to be subtracted - or 0.0 when the well has
        `min_negative` negative values or fewer.
    """
    fraction = _check_fraction(fraction)
    v = np.asarray(values, dtype=float)
    negatives = np.sort(v[np.isfinite(v) & (v < 0)])
    if negatives.size <= min_negative:
        return 0.0
    tail = negatives[:max(1, int(np.ceil(fraction * negatives.size)))]
    return float(np.median(tail))


def _well_keys(df: pd.DataFrame, well_keys) -> list[str]:
    keys = [well_keys] if isinstance(well_keys, str) else list(well_keys)
    if "well" not in keys:
        raise ValueError(f"well_keys must include 'well', got {keys}")
    refused = [k for k in keys if k in _POSITION_COLUMNS]
    if refused:
        raise ValueError(f"{', '.join(refused)} would give one offset per "
                         f"position; offsets are measured per well only")
    missing = [k for k in keys if k not in df.columns]
    if missing:
        raise KeyError(f"this table has no {', '.join(missing)} column; "
                       f"compile it with `compile_positions` or "
                       f"`load_experiment`, which add `well`")
    return keys


def well_offsets(df: pd.DataFrame, columns=None,
                 fraction: float = NEGATIVE_TAIL_FRACTION,
                 min_negative: int = MIN_NEGATIVE_CELLS,
                 well_keys=("well",), verbose: bool = True) -> pd.DataFrame:
    """Measure every channel's offset, one per well. **Step 1.**

    Nothing is changed and nothing is written.

    Parameters
    ----------
    df : DataFrame
        A compiled plate, from `compile_positions` or `load_experiment`.
    columns : sequence of str, optional
        Signal columns to measure. Defaults to `channel_columns(df)` - the best
        column of each channel present.
    fraction, min_negative
        Passed to `negative_tail_offset`.
    well_keys : sequence of str
        What identifies a well. `("well",)` for one plate; add the plate's own
        column, e.g. `("experiment", "well")`, when several plates are pooled
        in one table so that their A01s are not measured together. `position`,
        `stem` and `site` are refused.
    verbose : bool
        Report each column's offset range and the wells left unshifted.

    Returns
    -------
    DataFrame, one row per well per signal, with

    ============== ==========================================================
    `well_keys`    the well's identity
    signal         the column the row's offset belongs to
    positions      imaging positions pooled into the well, where known
    n_cells        cells with a finite value
    n_negative     of those, how many are below zero
    n_tail         how many of the most negative values the offset is the
                   median of; 0 when the well was not shifted
    offset         what `apply_well_offsets` subtracts (<= 0)
    applied        True when the well had more than `min_negative` negatives
    negative_after cells still below zero once the offset is subtracted
    ============== ==========================================================

    `.attrs` carries the settings, so `apply_well_offsets` and
    `plot_well_offsets` can be called with the table alone.
    """
    fraction = _check_fraction(fraction)
    keys = _well_keys(df, well_keys)
    columns = list(columns) if columns is not None else channel_columns(df)
    if not columns:
        raise ValueError(
            f"no usable signal column in this table; looked for "
            f"{', '.join(f'<ch>{s}' for s in SIGNAL_SUFFIXES)} for "
            f"{', '.join(CHANNELS)}")
    absent = [c for c in columns if c not in df.columns]
    if absent:
        raise KeyError(f"no column(s) {', '.join(absent)} in this table")

    rows = []
    for column in columns:
        for key, block in df.groupby(keys, sort=True, observed=True):
            key = key if isinstance(key, tuple) else (key,)
            v = block[column].to_numpy(dtype=float)
            v = v[np.isfinite(v)]
            n_negative = int((v < 0).sum())
            applied = n_negative > min_negative
            offset = negative_tail_offset(v, fraction, min_negative)
            rows.append(dict(zip(keys, key)) | {
                "signal": column,
                "positions": (int(block["position"].nunique())
                              if "position" in block.columns else np.nan),
                "n_cells": int(v.size),
                "n_negative": n_negative,
                "n_tail": (int(max(1, np.ceil(fraction * n_negative)))
                           if applied else 0),
                "offset": offset,
                "applied": applied,
                "negative_after": int((v - offset < 0).sum()),
            })
    offsets = pd.DataFrame(rows)

    if verbose:
        for column in columns:
            block = offsets[offsets["signal"] == column]
            label = block[keys].astype(str).agg("_".join, axis=1)
            idle = sorted(label[~block["applied"]])
            note = (f"; not shifted (<= {min_negative} negative cells): "
                    f"{', '.join(idle)}" if idle else "")
            print(f"{column}: {len(block)} well(s), offsets "
                  f"{block['offset'].min():+.2f} to "
                  f"{block['offset'].max():+.2f}{note}")

    offsets.attrs.update({"columns": columns, "well_keys": keys,
                          "fraction": fraction, "min_negative": min_negative})
    return offsets


def apply_well_offsets(df: pd.DataFrame, offsets: pd.DataFrame, columns=None,
                       suffix: str = ZEROED_SUFFIX) -> pd.DataFrame:
    """Subtract `well_offsets`' table, one new column per signal. **Step 2.**

    Adds `<column><suffix>` and leaves every original column alone, which is
    what makes a re-run safe: the offsets are always measured from the
    original, never from an already-zeroed one.

    A row whose well has no offset gets NaN in the new column and is reported
    - leaving it at its unshifted value would put two different zeros in one
    column.
    """
    keys = list(offsets.attrs.get("well_keys") or ["well"])
    columns = list(columns) if columns is not None else list(
        offsets.attrs.get("columns") or offsets["signal"].unique())

    out = df.copy()
    for column in columns:
        table = offsets.loc[offsets["signal"] == column, keys + ["offset"]]
        if table.empty or column not in out.columns:
            warnings.warn(f"no offsets for {column!r}, or no such column; "
                          f"not zeroed")
            continue
        if table.duplicated(keys).any():
            raise ValueError(f"the offsets table has more than one {column!r} "
                             f"row per {keys}")
        # A left merge keeps the rows of `out` in order, so the offsets line up
        # by position and not by index - a compiled table's index is not
        # guaranteed unique.
        offset = out[keys].merge(table, on=keys, how="left",
                                 validate="many_to_one")["offset"].to_numpy()
        missing = np.isnan(offset)
        if missing.any():
            labels = sorted(out.loc[missing, keys].astype(str)
                            .agg("_".join, axis=1).unique())
            warnings.warn(
                f"{int(missing.sum())} row(s) in {len(labels)} well(s) have no "
                f"offset and are NaN in {column}{suffix}: "
                f"{', '.join(labels[:5])}{' ...' if len(labels) > 5 else ''}")
        out[f"{column}{suffix}"] = out[column].to_numpy(dtype=float) - offset
    return out


def plot_well_offsets(df: pd.DataFrame, offsets: pd.DataFrame,
                      column: str | None = None, suffix: str = ZEROED_SUFFIX,
                      show_points: bool = True,
                      ylim_quantiles: tuple[float, float] = (0.005, 0.95),
                      figsize=None):
    """One signal, well by well, before and after the offset.

    Two panels down one well axis. The top one is the signal as compiled, with
    each well's offset marked; the bottom one is the zeroed signal. Each well's
    x label carries how many of its cells were negative before and after.

    The y-window is cut at `ylim_quantiles` of the pooled signal and shared by
    both panels, so the change between them is a shift and not a rescale.
    """
    column = column or (offsets.attrs.get("columns") or [None])[0]
    if column is None:
        raise ValueError("pass `column`")
    keys = list(offsets.attrs.get("well_keys") or ["well"])
    zeroed = f"{column}{suffix}"
    data = df if zeroed in df.columns else apply_well_offsets(
        df, offsets, columns=[column], suffix=suffix)
    data = data.copy()
    data["_well"] = data[keys].astype(str).agg("_".join, axis=1)

    table = offsets[offsets["signal"] == column].copy()
    table["_well"] = table[keys].astype(str).agg("_".join, axis=1)
    table = table.set_index("_well")
    labels = sorted(table.index)
    ticks = [f"{w}\n{table.loc[w, 'n_negative']}->{table.loc[w, 'negative_after']}"
             for w in labels]

    pooled = pd.concat([data[column], data[zeroed]]).to_numpy(dtype=float)
    pooled = pooled[np.isfinite(pooled)]
    low, high = np.quantile(pooled, ylim_quantiles)
    low = min(low, float(table["offset"].min()))
    pad = 0.06 * max(high - low, 1.0)

    fig, axes = plt.subplots(2, 1, sharex=True, sharey=True,
                             figsize=figsize or (0.7 * len(labels) + 5, 8))
    for ax, value, title in ((axes[0], column, "as compiled"),
                             (axes[1], zeroed, "after the well offset")):
        sns.boxplot(data=data, x="_well", y=value, order=labels, ax=ax,
                    color="0.85", showfliers=False, width=0.7, linewidth=1)
        if show_points:
            sns.stripplot(data=data, x="_well", y=value, order=labels, ax=ax,
                          color="0.25", size=2, alpha=0.35, jitter=0.28)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(f"{value} - {title}", fontsize=11)
        ax.set_ylabel("signal (a.u.)")
        ax.set_xlabel("")
    axes[0].plot(range(len(labels)), [table.loc[w, "offset"] for w in labels],
                 "o", color="crimson", mfc="white", mew=1.6, ms=7, zorder=5,
                 label=f"offset: median of the most negative "
                       f"{offsets.attrs.get('fraction', NEGATIVE_TAIL_FRACTION):.0%} "
                       f"of negative cells")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    axes[0].set_ylim(low - pad, high + pad)
    axes[1].set_xticks(range(len(labels)))
    axes[1].set_xticklabels(ticks, fontsize=8)
    axes[1].set_xlabel(f"well ({', '.join(keys)}); negative cells before->after")
    fig.tight_layout()
    return fig


def _write_summary_columns(path: Path, frame: pd.DataFrame, columns: list[str],
                           record: pd.DataFrame) -> None:
    """Add `columns` to a summary workbook's Summary sheet, in place.

    Every other sheet - `file_data`, `parameters`, `quality`, and whatever a
    later tool added - is left untouched, which is why this edits the workbook
    through openpyxl instead of reading it into pandas and writing it back out:
    a round trip through `read_excel`/`ExcelWriter` re-types every cell of
    every sheet and turns each one's index into an `Unnamed: 0` column.

    The two things the per-position alignment wrote, the `*_well_aligned`
    columns and the `well_alignment` sheet, are removed, so a file never
    carries a column zeroed the old way beside one zeroed this way.

    `frame` is this position's rows in file order, carrying `particle` and the
    new columns. The particle column is checked against the sheet rather than
    trusted, so a summary rewritten between compiling and writing fails loudly
    instead of having another position's numbers pasted into it.
    """
    import openpyxl

    book = openpyxl.load_workbook(path)
    if SUMMARY_SHEET not in book.sheetnames:
        raise KeyError(f"{path.name} has no {SUMMARY_SHEET!r} sheet")
    sheet = book[SUMMARY_SHEET]

    if sheet.max_row - 1 != len(frame):
        raise ValueError(f"{path.name}: {sheet.max_row - 1} rows in "
                         f"{SUMMARY_SHEET}, {len(frame)} compiled - the file "
                         f"changed since it was read")
    header = [cell.value for cell in sheet[1]]
    if "particle" in header:
        on_disk = [sheet.cell(row=r, column=header.index("particle") + 1).value
                   for r in range(2, sheet.max_row + 1)]
        if list(frame["particle"]) != list(on_disk):
            raise ValueError(f"{path.name}: the particle column does not match "
                             f"what was compiled; not written")

    # Right to left, so deleting one does not move the next one's index.
    for index in sorted((i + 1 for i, name in enumerate(header)
                         if isinstance(name, str)
                         and name.endswith(LEGACY_SUFFIX)), reverse=True):
        sheet.delete_cols(index)
    if LEGACY_SHEET in book.sheetnames:
        del book[LEGACY_SHEET]
    header = [cell.value for cell in sheet[1]]

    for column in columns:
        index = (header.index(column) + 1 if column in header
                 else sheet.max_column + 1)
        if column not in header:
            header.append(column)
        sheet.cell(row=1, column=index, value=column)
        for offset, value in enumerate(frame[column], start=2):
            # `.value =`, not `cell(..., value=...)`: openpyxl reads a None
            # there as "no value supplied" and leaves the cell alone, which on
            # a re-run would keep a stale number under a name that now
            # promises a NaN.
            sheet.cell(row=offset, column=index).value = (
                None if pd.isna(value) else float(value))

    # The provenance sheet is rebuilt rather than appended to: it describes the
    # columns now in the file, and two runs' worth of it would not say which.
    if OFFSET_SHEET in book.sheetnames:
        del book[OFFSET_SHEET]
    note = book.create_sheet(OFFSET_SHEET)
    note.append(list(record.columns))
    for row in record.itertuples(index=False):
        note.append([None if pd.isna(v) else
                     (v if isinstance(v, (int, float, str)) else str(v))
                     for v in row])
    book.save(path)


def correct_wells(source, columns=None,
                  fraction: float = NEGATIVE_TAIL_FRACTION,
                  min_negative: int = MIN_NEGATIVE_CELLS,
                  suffix: str = ZEROED_SUFFIX, write: bool = True,
                  file_suffix: str = "", pattern: str = "*phs.tif",
                  platemap: Path | str | None = None, verbose: bool = True,
                  ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Shift every well of a plate onto its zero, and write it into the summaries.

    The whole correction in one call:

        df, offsets = agg.correct_wells(root)

    which compiles the plate, measures one negative-tail offset per well per
    channel, and adds the shifted signal to every `*_summary.xlsx` as
    `<signal>_zeroed` beside the signal it came from.

    The rows are the summaries as written - `filter_summary` is not applied -
    so every particle in every file gets a value, and the offset is measured
    from every particle the well has.

    Parameters
    ----------
    source : Path | str | cellaap_analysis.analysis
        The plate's root folder, as for `load_experiment`.
    columns : sequence of str, optional
        Signal columns to correct. Defaults to the best column of each channel
        present - `<ch>_corrected` where signal_correction ran.
    fraction, min_negative
        Passed to `well_offsets`; see `negative_tail_offset`.
    suffix : str
        Appended to each column's name for the zeroed one.
    write : bool
        Add the columns to the summary workbooks. False computes and returns
        everything without touching disk, which is how to look at the offsets
        before believing them.
    file_suffix : str
        Which summary variant to read and write: `""` for `*_summary.xlsx`, or
        e.g. `"_dead"` for the copies `augment_dead_label.py --suffix` writes.
    pattern, platemap, verbose
        As for `platemap_positions`.

    Returns
    -------
    (df, offsets) : (DataFrame, DataFrame)
        `df` is the compiled plate with the zeroed columns added - the same
        numbers that went into the files. `offsets` is `well_offsets`' table.

    Notes
    -----
    Written in place. Each summary gets an `offset_correction` sheet naming
    what was subtracted from its well and under what settings, so a file says
    on its own what its `_zeroed` columns mean. Re-running overwrites both; the
    original signal columns are never touched, so a second run measures the
    same offsets as the first. `*_well_aligned` columns and the
    `well_alignment` sheet left by the per-position alignment are removed.
    """
    positions = platemap_positions(source, pattern=pattern, platemap=platemap,
                                   verbose=verbose)
    df = compile_positions(positions, suffix=file_suffix, verbose=verbose)
    if df.empty:
        raise ValueError(f"no *_summary{file_suffix}.xlsx found under "
                         f"{_resolve_root(source)}")

    columns = list(columns) if columns is not None else channel_columns(df)
    if verbose:
        print(f"\nzeroing {', '.join(columns)} in {df['well'].nunique()} "
              f"well(s), one offset per well")
    offsets = well_offsets(df, columns=columns, fraction=fraction,
                           min_negative=min_negative, verbose=verbose)
    zeroed = apply_well_offsets(df, offsets, columns=columns, suffix=suffix)
    new_columns = [f"{c}{suffix}" for c in columns
                   if f"{c}{suffix}" in zeroed.columns]

    if not write:
        if verbose:
            print(f"\nnothing written (write=False); "
                  f"{', '.join(new_columns)} are on the returned table only")
        return zeroed, offsets

    settings = {"method": "median of the most negative fraction of negative "
                          "cells, per well",
                "fraction": fraction, "min_negative": min_negative,
                "zeroed_suffix": suffix}
    written = 0
    for pos in positions:
        rows = zeroed[zeroed["stem"] == pos.stem]
        if rows.empty:
            continue
        path = pos.summary_path(file_suffix)
        if path is None:
            if verbose:
                print(f"  ! {pos.stub}: summary vanished, not written")
            continue
        record = offsets[offsets["well"] == pos.well].copy()
        for name, value in settings.items():
            record[name] = value
        _write_summary_columns(path, rows[["particle"] + new_columns],
                               new_columns, record)
        written += 1
        if verbose:
            print(f"{pos.stub}: {', '.join(new_columns)} written "
                  f"({len(rows)} particles)")

    if verbose:
        print(f"\n{written} summary file(s) updated in place; each carries an "
              f"{OFFSET_SHEET!r} sheet saying what was subtracted")
    return zeroed, offsets


# ---------------------------------------------------------------------------
# Fitting and plotting
# ---------------------------------------------------------------------------

def fit_model(xy_data: pd.DataFrame, plot: bool = True, quant_fraction=None,
              bin_size=None):
    '''
    Function to fit the dose-response data with a 4-parameter sigmoid.
    Bin range is determined by quantiles. Default is 0.025 and 0.85. The data
    typically contain outliers on the high side, but not the low side. Hence the
    default values are aysmmetric. For the model to be applicable, the fluorescence
    signal must be background subtracted. Subtracting the smallest value in the
    column is the crude version and it takes its zero from one cell, which on a
    plate with a crowded position is the most over-corrected cell there is.
    `well_offsets` then `apply_well_offsets` (or `correct_wells`) does the job
    per well, from the tail of its negative cells rather than from one point -
    see "Putting each well's zero back where it belongs".

    Inputs:
    xy_data        - dataframe w/ dose as the first column and response as
                     the second column
    plot           - Boolean to enable plotting
    quant_fraction - quantiles to determine bin range;
    bin_size       - size of each bin, default is 2.5 (empirical)
    Outputs:
    xy_data        - the input dataframe with bin labels added as a new column
    bin_means      - per-bin means
    bin_stderrs    - per-bin standard errors
    bin_sizes      - number of points per bin
    fit_pars       - dictionary containing fit parameters
    '''

    xy_data.dropna(inplace=True)
    if quant_fraction is None:
        quant_fraction = [0.025, 0.85]
    quants = np.round(xy_data.iloc[:,0].quantile(quant_fraction)).tolist()

    #
    if bin_size is None:
        bin_size = 2.5
    bins   = np.arange(0.5*quants[0], 1.5*quants[-1], bin_size).tolist()

    labels, _ = pd.cut(xy_data.iloc[:, 0], bins, retbins=True)
    xy_data["bins"] = labels

    bin_means = xy_data.groupby("bins").mean()
    bin_sizes = xy_data.groupby("bins").size()
    bin_stderrs = xy_data.groupby("bins").std()
    bin_stderrs.iloc[:,0] /= bin_sizes**0.5
    bin_stderrs.iloc[:,1] /= bin_sizes**0.5
    bin_means.dropna(inplace=True) # Some of the bins may not have any data
    bin_stderrs.dropna(inplace=True)


    fits, _ = curve_fit(sigmoid_4par, bin_means.iloc[:,0], bin_means.iloc[:,1],
                        p0 = [bin_means.iloc[:,1].min(), bin_means.iloc[:,1].max(),
                              5, (quants[0] + quants[-1])/ 4
                             ],
                        # sigma = bin_stderrs.iloc[:,1].to_numpy(),
                        maxfev = 10000
                       )

    fit_values = { 'min_duration' : fits[0],
                   'max_duration' : fits[1],
                   'Hill_exponent': fits[2],
                   'EC50'         : fits[3]
                 }
    if plot:
        fig, ax = plt.subplots(1,1, figsize=(8,6))
        sns.scatterplot(x=xy_data.iloc[:,0], y=xy_data.iloc[:,1],
                        ax=ax, alpha=0.1,
                        color="gray", edgecolor="None", size=1,
                        )

        sns.scatterplot(x = bin_means.iloc[:,0], y = bin_means.iloc[:,1],
                        color='w', edgecolor="blue", marker='s', linewidth=1,
                        label = "binned mean values")

        x_range = np.arange(0,1.5*quants[-1])
        sns.lineplot(x=x_range, y=sigmoid_4par(x_range,
                                               fit_values['min_duration'],
                                               fit_values['max_duration'],
                                               fit_values["Hill_exponent"],
                                               fit_values["EC50"]),
                                               ax = ax,
                                               markers='',
                                               color='b',
                                               label="Hill sigmoid fit")
        ax.set_xlabel("eSAC dosage (a.u.)")
        ax.set_ylabel("Time in mitosis (x 10 min)")
        ax.set_xlim(xmax=x_range[-1], xmin=x_range[0])
        y_quant = float(np.round(xy_data.iloc[:,1].quantile(0.99)))
        ax.set_ylim(0.0, y_quant)

    return xy_data, bin_means, bin_stderrs, bin_sizes, fit_values


def sigmoid_4par(x, base, top, exponent, ec50):

    """
    4-parameter Hill (sigmoid) function.

    This implements a common dose-response parameterization where
    `base` and `top` are the lower and upper asymptotes, `exponent`
    is the Hill coefficient, and `ec50` is the x value at half-max.

    Parameters
    ----------
    x : array-like or float
        Independent variable(s).
    base : float
        Minimum (baseline) value of the function.
    top : float
        Maximum (top) value of the function.
    exponent : float
        Hill coefficient (controls slope/steepness).
    ec50 : float
        Half-maximal effective concentration (EC50).

    Returns
    -------
    array-like or float
        The evaluated sigmoid at `x`.
    """

    return base + (top - base)*(x**exponent)/(x**exponent+ec50**exponent)


def export_to_excel_by_col(df, output_path, by_column="code", root_folder=None,
                           wellmap_path=None):
    """
    Export a dataframe to Excel with separate sheets for each code value and a summary sheet.

    Parameters:
    -----------
    df : pd.DataFrame
        The dataframe to export
    output_path : str or Path
        Path where the Excel file will be saved
    by_column : str, default="code"
        Column name to group by for sheet names
    root_folder : str or Path, optional
        Path to root folder to include in summary
    wellmap_path : str or Path, optional
        Path to wellmap file to include in summary
    """
    # Create summary data
    summary_data = []
    summary_data.append({"Field": "Root Folder", "Value": str(root_folder) if root_folder else "Not provided"})
    summary_data.append({"Field": "Wellmap Path", "Value": str(wellmap_path) if wellmap_path else "Not provided"})

    summary_df = pd.DataFrame(summary_data)

    # Build code summary table
    code_summary = []
    for code_value in sorted(df[by_column].unique()):
        code_df = df[df[by_column] == code_value]
        unique_wells = ", ".join(sorted(code_df["well"].unique().astype(str)))
        unique_positions = ", ".join(sorted(code_df["position"].unique().astype(str)))
        code_summary.append({
            "Code": str(code_value),
            "Unique Wells": unique_wells,
            "Unique Positions": unique_positions
        })

    code_summary_df = pd.DataFrame(code_summary)

    # Write to Excel
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="Summary", index=False)
        code_summary_df.to_excel(writer, sheet_name="Code Summary", index=False)

        # Create sheets for each code
        for code_value in sorted(df[by_column].unique()):
            sheet_df = df[df[by_column] == code_value]
            code_value = str(code_value).replace("/", " ")  # Replace slashes to avoid Excel sheet name issues
            if len(code_value) > 31:
                print(f"Warning: Code value '{code_value}' exceeds Excel sheet name limit. Truncating to 31 characters.")
                code_value = code_value[:31]
            sheet_df.to_excel(writer, sheet_name=code_value, index=False)

    print(f"Excel file saved to {output_path}")
    print(f"Created {len(df[by_column].unique()) + 2} sheets (1 metadata + 1 code summary + {len(df[by_column].unique())} data sheets)")

    return
