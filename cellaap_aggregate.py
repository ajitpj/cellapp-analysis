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

A compiled table still carries one plate-level defect: every position had its
background estimated on its own, so the positions do not share a zero. See
"Aligning the zero across positions" below - `baseline_offsets` measures the
disagreement, `plot_baseline_offsets` shows it, and `apply_baseline_offsets`
removes it. Three steps rather than one, so the correction can be looked at
before it is believed.

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
    "baseline_floor",
    "baseline_offsets",
    "apply_baseline_offsets",
    "plot_baseline_offsets",
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
# Aligning the zero across positions
#
# `signal_correction` measures each position's background from that position's
# own cell-free pixels, which is the right thing to do and is why
# `<ch>_corrected` beats the blank-well maps. It has one failure mode, and it
# is systematic:
# there have to BE cell-free pixels. As a field fills up, the estimator backs
# its exclusion ring off the cells (`_dilation_ladder`) to keep enough blocks
# measurable, and the closer it measures to a cell the more of that cell's
# out-of-focus halo it counts as medium. The background comes out too high, the
# corrected signal too low, and the error grows with confluence - so it differs
# between positions in one well, between wells, and between days.
#
# The size of it, on the 20260826 CycB plate. Across the four A01 positions -
# one well, one treatment, no GFP induced, so their dim cells are the same
# cells - the RAW floor spans 0.9 counts and the per-position-corrected floor
# spans 5.1, against a median signal of 23. The correction added
# between-position spread rather than removing it. Across the ten positions of
# the two uninduced wells the floor tracks the estimator's difficulty
# precisely: r = -0.83 against the position's peak background, -0.80 against
# its apparent background drift, and +0.79 against how wide an exclusion ring
# the estimator could hold (`correction_tools.surface_diagnostics` reports all
# three). One position, C01_s6 - the densest on the plate, ring backed off from
# 121 px to 30 - over-subtracts hard enough to put its dimmest cells at -300.
#
# What is done about it here is deliberately not a better background estimate -
# that would have to happen back at the images. It is the observation that a
# cell with no fluorophore reads the same number everywhere, because that
# number is a property of the microscope and not of the well. So the bottom of
# each position's distribution is a landmark that SHOULD line up, and how far
# it fails to line up is the residual background error, directly measured.
# Subtracting that per-position constant is the correction.
#
# The estimate and the subtraction are separate calls on purpose. The offsets
# come back as an ordinary DataFrame that can be read, plotted, edited or
# thrown away before anything touches the data.
#
# It removes an OFFSET and nothing else. Three things it therefore cannot do:
#
#   * fix a within-position error. C01_s6 above is over-subtracted by a
#     different amount in different parts of the field and at different times;
#     one constant cannot straighten that, and `plot_baseline_offsets` will
#     still show it with a long low tail after alignment. Such a position is
#     to be dropped, not aligned.
#   * survive a unit with no dim cells in it. The floor is only a landmark
#     where some of the cells really are at zero. In a well where every cell
#     expresses, the bottom of the distribution is a biological number, and
#     aligning on it flattens a real difference into nothing. This is the one
#     way to do actual damage with these functions, which is why nothing is
#     applied until you have looked.
#   * make two experiments comparable on its own. A different exposure or
#     laser power rescales the signal as well as shifting it; aligning the zero
#     is necessary for pooling repeats but is not by itself sufficient.
# ---------------------------------------------------------------------------

# Between the 2nd and 15th percentile: high enough to clear the over-corrected
# tail a crowded position leaves at the very bottom, low enough to stay inside
# the non-expressing population. Averaging that band rather than reading one
# order statistic is what makes it steady at the ~100 tracks a position has -
# bootstrap SE 1.2 counts against 1.5 for a bare 5th percentile on the
# 20260826 plate, and a smaller within-well spread on every well of it.
BASELINE_TRIM = (0.02, 0.15)

# Reasons a unit's offset should be looked at before it is used.
FLAG_FEW_CELLS = "few cells"
FLAG_OUTLIER = "outlier floor"
FLAG_LOW_TAIL = "low tail"


def _half_sample_mode(values: np.ndarray, min_n: int = 8) -> float:
    """Densest point of a distribution, by recursive half-sample shrinking.

    Bickel & Fruehwirth's estimator: repeatedly keep the half of the sorted
    sample that spans the smallest range. Unlike a histogram mode it needs no
    bin width, and unlike a quantile it ignores the shape of the tails
    entirely - which is what makes it the estimator to use when a unit holds a
    large non-expressing population and a long bright one.
    """
    x = np.sort(values)
    while len(x) > min_n:
        n = len(x)
        half = n // 2
        i = int(np.argmin(x[half:] - x[:n - half]))
        x = x[i:i + half + 1]
    return float(np.median(x))


def baseline_floor(values, estimator: str = "trimmed", q: float = 0.05,
                   trim: tuple[float, float] = BASELINE_TRIM) -> float:
    """Where the dim cells of one unit sit: the landmark that should line up.

    Parameters
    ----------
    values : array-like
        One unit's per-cell signal. NaNs are dropped.
    estimator : {"trimmed", "quantile", "mode"}
        ``trimmed``
            mean of the values between the two `trim` quantiles. The default,
            for the reasons on `BASELINE_TRIM`.
        ``quantile``
            the `q`-th percentile. One order statistic, so noisier, but it is
            the number people already read off a boxplot and it makes the
            offsets easy to check by eye.
        ``mode``
            `_half_sample_mode`. Use it when most cells in a unit are
            non-expressing, so the peak of the distribution IS the zero;
            it ignores both tails, including an over-corrected one.

    Returns
    -------
    float, or NaN for an empty unit.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan")
    if estimator == "quantile":
        return float(np.quantile(v, q))
    if estimator == "trimmed":
        lo, hi = np.quantile(v, trim)
        band = v[(v >= lo) & (v <= hi)]
        return float(band.mean()) if band.size else float(hi)
    if estimator == "mode":
        return _half_sample_mode(v)
    raise ValueError(f"unknown estimator {estimator!r}; expected 'trimmed', "
                     f"'quantile' or 'mode'")


def baseline_offsets(df: pd.DataFrame, column: str, by="stem",
                     estimator: str = "trimmed", q: float = 0.05,
                     trim: tuple[float, float] = BASELINE_TRIM,
                     reference="median", within=None,
                     min_cells: int = 20, n_boot: int = 200,
                     seed: int = 0, flag_z: float = 3.0,
                     max_tail_drop: float = 1.0) -> pd.DataFrame:
    """Measure how far each unit's zero sits from the common one. **Step 1.**

    Nothing is changed here. The result is a table you look at - and, when a
    unit turns out to have no dim cells to measure, edit - before
    `apply_baseline_offsets` uses it.

    Parameters
    ----------
    df : DataFrame
        A compiled summary, from `load_experiment` or `compile_positions`.
    column : str
        The signal to align, e.g. `"GFP_corrected"`. One channel per call; the
        tables from two calls concatenate, since each carries a `signal` column.
    by : str or sequence of str
        What a unit is. `"stem"` (the default) is one imaging position, which
        is the level the background was estimated at and therefore the level
        the error lives at. `"well"` pools the sites of a well; `["experiment",
        "code"]` pools a whole condition of a whole plate. Coarser units give a
        steadier floor and leave more residual error behind.
    estimator, q, trim
        Passed to `baseline_floor`.
    reference : {"median", "min", "zero"}, a unit label, or a number
        The zero everything is moved onto.

        ``median``  the median floor over units that clear `min_cells`. The
                    default: it corrects the disagreement between units without
                    claiming to know the absolute zero, so the plate's overall
                    level - and any comparison against an earlier analysis of
                    it - is left where it was.
        ``min``     the lowest floor. Assumes the least-corrected unit is the
                    most trustworthy, which is true when the errors are all
                    over-subtraction of background.
        ``zero``    put every floor at 0. Only when a genuinely non-expressing
                    population is present in every unit; then the corrected
                    number really is fluorophore, and dose-response fits that
                    need a positive zero-dose value (see `fit_model`) can use
                    it directly.
        a label     align on one named unit, e.g. an untreated control well.
                    Matched against the unit label - the `by` columns joined
                    with `_` when there is more than one.
        a number    that value, whatever the data say.
    within : str or sequence of str, optional
        Compute a separate reference inside each of these groups instead of one
        for the whole table. This is how you choose what the alignment is
        allowed to touch, and it is the argument to think about:

        ``None``        one zero for everything. The strongest correction, and
                        the right one when every unit really should read the
                        same at zero - repeats of a plate, or wells that differ
                        only in a drug that does not touch the reporter.
        ``"code"``      align positions within each condition and leave the
                        conditions where they are. The conservative choice, and
                        the one to reach for when a treatment induces the
                        reporter, because it cannot flatten the induction. On
                        the 20260826 plate it leaves the three condition
                        medians within 0.2 counts and pulls the
                        position-to-position spread of the uninduced wells from
                        7.0 and 5.0 counts down to 4.3 and 4.2.
        ``"experiment"`` align each plate to its own median and leave the
                        plates' levels alone - for plates that are NOT expected
                        to share a zero, e.g. a changed exposure.
    min_cells : int
        Below this a unit is flagged `"few cells"`, and it is left out of a
        `"median"` or `"min"` reference. Its offset is still computed.
    n_boot : int
        Bootstrap resamples behind `floor_se`. 0 skips it and returns NaN.
    seed : int
        Bootstrap seed, so the same table gives the same standard errors.
    flag_z : float
        Flag a unit `"outlier floor"` when its floor is this many robust SDs
        (MAD-scaled, over the units in its `within` group) from the reference.
    max_tail_drop : float
        Flag a unit `"low tail"` when `tail_drop` exceeds this **and**
        `tail_z` puts it more than `flag_z` robust SDs above the other units.
        Both are needed: the ratio alone reads high on every unit of a channel
        whose bulk is tight, however healthy they are. With only two or three
        units there is nothing to be an outlier against, so nothing is flagged.

    Returns
    -------
    DataFrame, one row per unit, with

    ============== ==========================================================
    `by` columns   the unit's identity, plus any `within` columns
    signal         `column`, so tables for several channels concatenate
    n_cells        rows behind the floor
    floor          `baseline_floor` of this unit
    floor_se       bootstrap standard error of `floor`
    reference_floor the zero this unit is being moved onto
    offset         `floor - reference_floor`: the excess baseline this unit
                   carries. `apply_baseline_offsets` subtracts it.
    z              `(floor - reference_floor)` in robust SDs of the floors
    tail_drop      how far the unit's 1st percentile falls below its own
                   floor, in interquartile ranges. This is the column that
                   catches a position no constant can fix: on the 20260826
                   plate every healthy GFP position sits under 0.6 and the one
                   broken one at 3.4
    tail_z         `tail_drop` in robust SDs of the other units' `tail_drop`.
                   What makes the flag work across channels - see
                   `max_tail_drop`
    flag           `""`, or the reasons this row deserves a look, comma-joined
    ============== ==========================================================

    The call's settings are on `.attrs`, which is what lets
    `apply_baseline_offsets` and `plot_baseline_offsets` be called with just
    the table.
    """
    by = [by] if isinstance(by, str) else list(by)
    within = ([] if within is None else
              [within] if isinstance(within, str) else list(within))
    missing = [c for c in by + within + [column] if c not in df.columns]
    if missing:
        raise KeyError(f"{df.__class__.__name__} has no column(s) "
                       f"{', '.join(missing)}")

    rng = np.random.default_rng(seed)
    keys = by + [c for c in within if c not in by]

    rows = []
    for key, group in df.groupby(keys, sort=False, observed=True):
        key = key if isinstance(key, tuple) else (key,)
        v = group[column].to_numpy(dtype=float)
        v = v[np.isfinite(v)]
        floor = baseline_floor(v, estimator, q, trim)
        se = float("nan")
        if n_boot and v.size >= 5:
            draws = [baseline_floor(rng.choice(v, v.size, replace=True),
                                    estimator, q, trim)
                     for _ in range(n_boot)]
            se = float(np.std(draws))
        rows.append(dict(zip(keys, key)) |
                    {"signal": column, "n_cells": int(v.size),
                     "floor": floor, "floor_se": se, "_values": v})

    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError(f"no rows to measure a floor from in {column!r}")
    out["label"] = out[by].astype(str).agg("_".join, axis=1)

    # The reference, once per `within` group.
    group_keys = [c for c in within if c in out.columns]
    grouped = ([("", out)] if not group_keys
               else list(out.groupby(group_keys, sort=False, observed=True)))
    pieces = []
    for _, block in grouped:
        block = block.copy()
        block["reference_floor"] = _baseline_reference(block, reference,
                                                       min_cells)
        block["offset"] = block["floor"] - block["reference_floor"]
        # MAD over the floors of this block, as the yardstick for "far".
        spread = float(np.median(np.abs(block["floor"]
                                        - np.median(block["floor"])))) * 1.4826
        block["z"] = (block["offset"] / spread if spread > 0
                      else np.where(block["offset"] == 0, 0.0, np.inf))
        pieces.append(block)
    out = pd.concat(pieces, ignore_index=True)

    # How far the very bottom of a unit falls below its own floor, in units of
    # that unit's interquartile range. The one thing a constant offset cannot
    # repair is a unit over-subtracted by different amounts in different parts
    # of the field, and that shows up here and nowhere else: the floor moves a
    # little, the tail underneath it collapses. On the 20260826 plate the
    # healthy positions sit at 0.25-0.60 and the one broken position at 3.4.
    tail_drop = []
    for floor, v in zip(out["floor"], out["_values"]):
        if v.size < 5:
            tail_drop.append(float("nan"))
            continue
        bottom, q1, q3 = np.quantile(v, [0.01, 0.25, 0.75])
        iqr = float(q3 - q1)
        tail_drop.append((floor - bottom) / iqr if iqr > 0 else float("nan"))
    out["tail_drop"] = tail_drop

    # `tail_drop` is a ratio to the unit's own interquartile range, so a
    # channel whose bulk is tight reads high on every unit with nothing wrong:
    # a near-saturated stain has a narrow IQR and a few dim cells under it, and
    # on the 20260826 plate that put 10 of 15 Cy5 positions over an absolute
    # threshold that caught exactly one GFP position. What marks a position no
    # offset can repair is a tail unlike the OTHER units of the same channel,
    # so the flag needs both: over `max_tail_drop`, and an outlier among its
    # peers. The broken GFP position sits 17 robust SDs out; the worst Cy5 one
    # sits at 1.7 and is left alone. `tail_z` is that second number, reported
    # so a large `tail_drop` with no flag explains itself.
    finite = np.asarray([t for t in tail_drop if np.isfinite(t)], dtype=float)
    tail_mid = float(np.median(finite)) if finite.size else np.nan
    tail_mad = (float(np.median(np.abs(finite - tail_mid))) * 1.4826
                if finite.size else 0.0)
    out["tail_z"] = [((t - tail_mid) / tail_mad if tail_mad > 0 else 0.0)
                     if np.isfinite(t) else np.nan for t in tail_drop]

    out["flag"] = [
        ", ".join(f for f in (
            FLAG_FEW_CELLS if n < min_cells else "",
            FLAG_OUTLIER if np.isfinite(z) and abs(z) > flag_z else "",
            FLAG_LOW_TAIL if (np.isfinite(t) and t > max_tail_drop
                              and np.isfinite(tz) and tz > flag_z) else "",
        ) if f)
        for n, z, t, tz in zip(out["n_cells"], out["z"], out["tail_drop"],
                               out["tail_z"])]

    out = out.drop(columns="_values")
    order = (by + [c for c in within if c not in by] +
             ["label", "signal", "n_cells", "floor", "floor_se",
              "reference_floor", "offset", "z", "tail_drop", "tail_z",
              "flag"])
    out = out[[c for c in order if c in out.columns]]
    out.attrs.update({"column": column, "by": by, "within": within,
                      "estimator": estimator, "q": q, "trim": tuple(trim),
                      "reference": reference, "min_cells": min_cells,
                      "flag_z": flag_z, "max_tail_drop": max_tail_drop})
    return out


def _baseline_reference(block: pd.DataFrame, reference, min_cells: int) -> float:
    """The zero for one `within` group. See `baseline_offsets`'s `reference`."""
    if isinstance(reference, (int, float)) and not isinstance(reference, bool):
        return float(reference)
    usable = block[block["n_cells"] >= min_cells]
    if usable.empty:
        usable = block
    if reference == "median":
        return float(np.median(usable["floor"]))
    if reference == "min":
        return float(np.min(usable["floor"]))
    if reference == "zero":
        return 0.0
    hit = block[block["label"] == str(reference)]
    if hit.empty:
        raise ValueError(
            f"reference {reference!r} is neither 'median', 'min', 'zero', a "
            f"number, nor one of the units {sorted(block['label'])}")
    return float(hit["floor"].iloc[0])


def apply_baseline_offsets(df: pd.DataFrame, offsets: pd.DataFrame,
                           column: str | None = None, by=None,
                           suffix: str = "_aligned",
                           columns=None,
                           skip_flagged: bool = False) -> pd.DataFrame:
    """Subtract the measured offsets. **Step 3** (step 2 is looking at them).

    Adds `<column><suffix>`; the original column is never touched, so the two
    can be plotted against each other and the alignment undone by dropping a
    column.

    Parameters
    ----------
    df, offsets
        The compiled summary, and `baseline_offsets`' table for it. `column`
        and `by` default to what that call used.
    columns : sequence of str, optional
        Extra columns to shift by the same offsets. The obvious use is a
        per-cell standard deviation or a second summary of the same channel;
        a column of a DIFFERENT channel has its own baseline error and needs
        its own `baseline_offsets` call.
    skip_flagged : bool
        Treat a flagged unit's offset as 0 rather than applying it. The unit
        stays in the table, uncorrected. Off by default: a flag is an
        instruction to look, and if the look says the offset is wrong the row
        should be edited or dropped rather than silently neutralized.

    Returns
    -------
    A copy of `df`. A row whose unit has no offset gets NaN in the new column
    and is reported - leaving it at its unaligned value would put two different
    zeros in one column, which is the thing this exists to prevent.
    """
    column = column or offsets.attrs.get("column")
    by = by or offsets.attrs.get("by")
    if column is None or by is None:
        raise ValueError("pass `column` and `by`; this offsets table carries "
                         "no .attrs (a round trip through a spreadsheet drops "
                         "them)")
    by = [by] if isinstance(by, str) else list(by)
    targets = [column] + [c for c in (columns or []) if c != column]

    table = offsets[by + ["offset"]].copy()
    if skip_flagged and "flag" in offsets.columns:
        table.loc[offsets["flag"].astype(bool).to_numpy(), "offset"] = 0.0
    if table.duplicated(by).any():
        raise ValueError(f"the offsets table has more than one row per "
                         f"{by}; is it two channels concatenated? Filter it "
                         f"to one `signal` first")

    out = df.merge(table, on=by, how="left", validate="many_to_one")
    unmatched = out["offset"].isna()
    if unmatched.any():
        labels = sorted(out.loc[unmatched, by].astype(str)
                        .agg("_".join, axis=1).unique())
        warnings.warn(
            f"{int(unmatched.sum())} row(s) in {len(labels)} unit(s) have no "
            f"offset and are NaN in {column}{suffix}: "
            f"{', '.join(labels[:5])}{' ...' if len(labels) > 5 else ''}")
    for target in targets:
        if target not in out.columns:
            warnings.warn(f"{target!r} is not a column of this table; skipped")
            continue
        out[f"{target}{suffix}"] = out[target] - out["offset"]
    return out.drop(columns="offset")


def _strip_common_prefix(labels: list[str]) -> dict[str, str]:
    """Shorten `20260826_Hela CycB Oe BubR1 kd_A01_s2` to `A01_s2` for an axis.

    Position stems on one plate share everything but the last two fields, and
    the shared part is what makes the tick labels unreadable. Only whole
    underscore-separated fields are dropped, and only when every label keeps
    at least one; otherwise the labels come back untouched.
    """
    parts = [label.split("_") for label in labels]
    if len(labels) < 2:
        return {label: label for label in labels}
    n = 0
    while n < min(len(p) for p in parts) - 1 and len({p[n] for p in parts}) == 1:
        n += 1
    return {label: "_".join(p[n:]) for label, p in zip(labels, parts)}


def plot_baseline_offsets(df: pd.DataFrame, offsets: pd.DataFrame,
                          column: str | None = None, by=None,
                          suffix: str = "_aligned", hue: str | None = None,
                          order: str = "floor", palette: str = "colorblind",
                          show_points: bool = True, ylim_quantile: float = 0.95,
                          figsize=None):
    """Before, after, and the floors themselves. **Step 2.**

    Three panels down one shared unit axis:

    1. the signal as it stands, one box per unit, with each unit's measured
       floor as a marker and the reference as a dashed line. Units whose boxes
       are at different heights may just be different; units whose FLOORS are
       at different heights are misaligned, and this panel separates the two.
    2. the same after `apply_baseline_offsets`. The floors should now sit on
       the line. A unit that still spills below it is one no constant fixes.
    3. floor +/- bootstrap SE against the reference, flagged units in red.
       The size of the correction, against the noise in measuring it.

    `order="floor"` sorts units by floor, which puts the misaligned ones at
    the ends; `order="label"` keeps them in name order for reading off a plate.

    The top two panels are cut off at `ylim_quantile` of the pooled signal,
    always keeping every floor in view. Left at full range a handful of bright
    cells set the scale and the few counts this figure is about are invisible;
    the bright cells are not what is being judged here.
    """
    column = column or offsets.attrs.get("column")
    by = by or offsets.attrs.get("by")
    by = [by] if isinstance(by, str) else list(by)
    aligned = f"{column}{suffix}"

    data = df if aligned in df.columns else apply_baseline_offsets(
        df, offsets, column=column, by=by, suffix=suffix)
    data = data.copy()
    data["label"] = data[by].astype(str).agg("_".join, axis=1)
    table = offsets.copy()
    if "label" not in table.columns:
        table["label"] = table[by].astype(str).agg("_".join, axis=1)

    labels = (table.sort_values("floor")["label"].tolist() if order == "floor"
              else sorted(table["label"]))
    floor = table.set_index("label")["floor"]
    ref = table.set_index("label")["reference_floor"]
    err = table.set_index("label")["floor_se"]
    flagged = table.set_index("label")["flag"].astype(bool)

    # One y-window for both signal panels, so the shift between them is a
    # shift and not a rescale.
    pooled = pd.concat([data[column], data[aligned]]).to_numpy(dtype=float)
    pooled = pooled[np.isfinite(pooled)]
    top = float(np.quantile(pooled, ylim_quantile))
    bottom = float(min(floor.min(), ref.min()))
    pad = 0.08 * max(top - bottom, 1.0)
    window = (bottom - 2 * pad, top + pad)

    short = _strip_common_prefix(labels)
    fig, axes = plt.subplots(3, 1, sharex=True,
                             figsize=figsize or (0.55 * len(labels) + 5, 11),
                             gridspec_kw={"height_ratios": [3, 3, 2]})
    for ax, value, title, marker in (
            (axes[0], column, "as corrected per position", "measured floor"),
            (axes[1], aligned, "after aligning the zero", "aligned floor")):
        sns.boxplot(data=data, x="label", y=value, order=labels, ax=ax,
                    hue=hue, palette=palette if hue else None,
                    color=None if hue else "0.85", dodge=False,
                    showfliers=False, width=0.7, linewidth=1)
        if show_points:
            sns.stripplot(data=data, x="label", y=value, order=labels, ax=ax,
                          color="0.25", size=2, alpha=0.35, jitter=0.28)
        # The floors, drawn where they are measured: on the raw column for the
        # top panel, on the reference for the bottom one - after alignment
        # every unit's floor IS the reference, which is the claim being made.
        marks = ([floor[k] for k in labels] if value == column
                 else [ref[k] for k in labels])
        ax.plot(range(len(labels)), marks, "o", color="crimson", ms=6,
                mfc="white", mew=1.6, label=marker, zorder=5)
        for k, label in enumerate(labels):
            ax.hlines(ref[label], k - 0.45, k + 0.45, color="crimson", lw=1.2,
                      ls="--", zorder=4,
                      label="reference zero" if k == 0 else None)
        ax.set_title(f"{value} - {title}", fontsize=11)
        ax.set_ylabel("signal (a.u.)")
        ax.set_xlabel("")
        ax.set_ylim(*window)
        ax.legend(frameon=False, fontsize=8, loc="upper left", ncol=3)

    ax = axes[2]
    colors = ["crimson" if flagged[k] else "tab:blue" for k in labels]
    ax.errorbar(range(len(labels)), [floor[k] for k in labels],
                yerr=[err[k] if np.isfinite(err[k]) else 0 for k in labels],
                fmt="none", ecolor="0.4", capsize=3, zorder=1)
    ax.scatter(range(len(labels)), [floor[k] for k in labels], c=colors,
               s=45, zorder=2)
    ax.plot(range(len(labels)), [ref[k] for k in labels], color="crimson",
            lw=1.2, ls="--", label="reference zero")
    for k, label in enumerate(labels):
        if flagged[label]:
            ax.annotate(table.set_index("label")["flag"][label],
                        (k, floor[label]), textcoords="offset points",
                        xytext=(0, 9), ha="center", fontsize=7,
                        color="crimson")
    ax.set_ylabel("floor (a.u.)")
    ax.set_xlabel("unit (" + ", ".join(by) + ")")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([short[k] for k in labels], rotation=90)
    fig.tight_layout()
    return fig


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
    `baseline_offsets(df, column, reference="zero")` then
    `apply_baseline_offsets` does the same job per position, from the dim
    population rather than from one point - see "Aligning the zero across
    positions".

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
