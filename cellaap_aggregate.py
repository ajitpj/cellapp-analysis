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
list, and each fixes a way the old well-list path quietly got the wrong answer:

* **Well ranges and normalization.** `B01-B06` expands, `g3` becomes `G03`.
  The old `create_wellmap_dict` split only on commas and whitespace, so a range
  survived as the literal string `B01-B06` and matched no folder at all.
* **`skip` and blank-media rows drop out.** A `fluorobrite` well has no
  celltype, transfection or drug, so it used to form a junk `('', '', '')`
  group that then found no summaries, because blank wells are never segmented.
* **Position overrides stop double-counting.** With a `G03` row and a `G03_s9`
  row under a different drug, matching the substring `_G03_` puts site 9 in
  *both* groups. Precedence is resolved per position here, as the pipeline
  resolves it, so each position lands in exactly one group.

Entry points, in decreasing order of how much you have to say:

    import cellaap_aggregate as agg

    # the platemap does everything
    df = agg.load_experiment(root, expt_length=150, delta_t=10)

    # or the pieces, to filter or regroup in between
    positions = agg.platemap_positions(root)
    groups    = agg.group_positions(positions)
    raw       = agg.compile_positions(groups[("HeLa", "pEN2", "DMSO")])

The older well-list API (`create_wellmap_dict` -> `import_whole_expt_data`) is
still here and still takes the same arguments, so existing notebooks need only
their import line changed. It now understands the platemap's `well_ids` syntax,
which it previously did not.

This module imports matplotlib and seaborn at module scope. `cellaap_utils`
deliberately does not: it is imported by `cellaap_analysis` and therefore by
every analysis array task on the cluster, and none of them plot.
"""

from __future__ import annotations

import re
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
    "create_wellmap_dict",
    "wellmap_from_platemap",
    "compile_summaries",
    "import_filter_data_for_wells",
    "import_whole_expt_data",
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
        # Kept as the old compile_summaries defined it - the folder the root
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

    The one call that replaces `create_wellmap_dict` -> `import_whole_expt_data`:

        df = load_experiment(root, expt_length=150, delta_t=10)

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


# ---------------------------------------------------------------------------
# The older well-list API
#
# Kept working because notebooks use it. The behaviour differences are all
# fixes: well ranges expand, blank-media and skipped rows are dropped, and an
# empty result is an empty DataFrame rather than an exception.
# ---------------------------------------------------------------------------

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


def create_wellmap_dict(imported_wellmap: pd.DataFrame,
                        wellid_col_name: str = "well_ids") -> dict:
    """Map (celltype, transfection, drug) -> list of well ids.

    Takes the platemap as a DataFrame - see `read_platemap` - and groups its
    rows. Wells are expanded with `pipeline.expand_well_token`, so the full
    `well_ids` syntax works: single wells, comma/space lists, ranges within one
    plate row (`B01-B06`), and exact positions (`G03_s9`). Casing and zero
    padding are normalized, so `g3` and `G3` both become `G03`.

    Rows are dropped when `skip` is set or `role` names a blank-media well,
    if those columns are present: neither was ever analyzed, so including them
    only produces groups with no data behind them.

    Returns
    -------
    dict
        Keyed by (celltype, transfection, drug), values are ordered unique
        well ids.

    Notes
    -----
    A well-level entry matches every site in that well, so if the platemap also
    gives one site of that well its own row, the site appears in both groups.
    `platemap_positions` resolves that precedence properly; prefer it, or
    `wellmap_from_platemap`, which returns positions rather than wells.
    """
    for column in ("celltype", "transfection", "drug"):
        if column not in imported_wellmap.columns:
            raise KeyError(f"Required column not found in imported_wellmap: {column}")
    if wellid_col_name not in imported_wellmap.columns:
        raise KeyError(f"{wellid_col_name} column not found in imported_wellmap")

    df = imported_wellmap
    if "skip" in df.columns:
        df = df[~df["skip"].map(lambda v: pipeline.as_bool("" if pd.isna(v) else v))]
    if "role" in df.columns:
        blank = df["role"].map(
            lambda v: pipeline.ROLES.get(str("" if pd.isna(v) else v).strip().lower())
            in pipeline.MAP_ROLES)
        df = df[~blank]

    def parse_wells(value) -> list[str]:
        if isinstance(value, (list, tuple)):
            tokens = [str(v) for v in value]
        elif pd.isna(value):
            return []
        else:
            tokens = re.split(r"[,;\s]+", str(value))
        wells: list[str] = []
        for token in tokens:
            wells.extend(pipeline.expand_well_token(token))
        return wells

    wellmap_dict: dict = {}
    for key, group in df.groupby(["celltype", "transfection", "drug"]):
        wells: list[str] = []
        for value in group[wellid_col_name]:
            wells.extend(parse_wells(value))
        seen = set()
        wellmap_dict[key] = [w for w in wells if not (w in seen or seen.add(w))]

    return wellmap_dict


def wellmap_from_platemap(source, pattern: str = "*phs.tif",
                          platemap: Path | str | None = None,
                          verbose: bool = True) -> dict:
    """`create_wellmap_dict`'s result, but resolved to positions on disk.

    Same keys, but each value is a list of position stubs (`G03_s8`) that
    actually exist, rather than well ids. That makes the mapping unambiguous:
    a position appears under exactly one condition even when the platemap
    overrides a single site, and wells with no data do not appear at all.

    Feed it to `import_whole_expt_data` in place of `create_wellmap_dict`.
    """
    groups = group_positions(platemap_positions(
        source, pattern=pattern, platemap=platemap, verbose=verbose))
    return {key: [p.stub for p in members] for key, members in groups.items()}


def compile_summaries(cellapp_expt, wells: list, suffix: str = "") -> pd.DataFrame:
    """Concatenate the summary spreadsheets for a list of wells or positions.

    Parameters
    ----------
    cellapp_expt : cellaap_analysis.analysis | Path | str
        The analysis session whose `root_folder` holds the inference folders,
        or the folder itself.
    wells : list
        Well ids (`G03`) or position stubs (`G03_s8`). A well id matches every
        site in that well.
    suffix : str
        Summary variant, as for `compile_positions`.

    Returns
    -------
    pandas.DataFrame
        With `well`, `position` and `storage_location` added, as before. Empty
        if `wells` is empty or nothing matched - the previous version raised
        `UnboundLocalError` and `ValueError` respectively in those two cases.
    """
    root = _resolve_root(cellapp_expt)
    if not wells:
        return pd.DataFrame()

    # Rebuilt per call rather than cached on the session: folders appear as a
    # plate finishes, and a cached list from the first call goes stale.
    folders = [f for f in sorted(root.glob("*_inference")) if f.is_dir()]
    storage_location = str(root.parent)

    frames = []
    for well in wells:
        needle = f"_{well}_"
        matched = [f for f in folders if needle in f.name]
        if not matched:
            print(f"No inference folder found for {well}")
            continue
        for folder in matched:
            hits = sorted(folder.glob(f"*_summary{suffix}.xlsx"))
            if not hits:
                print(f"No summary file found for {folder.name}")
                continue
            df = pd.read_excel(hits[0], sheet_name=SUMMARY_SHEET)
            # Read the well and site back out of the folder name rather than
            # from `well`, which may be either a well id or a position stub.
            # Otherwise a stub lands in the `well` column and the site is
            # repeated: G03_s9 / s9.
            stub = pipeline.STUB_RE.search(folder.name)
            df["well"] = pipeline.normalize_well(stub.group("well")) if stub else well
            df["position"] = f"s{int(stub.group('site'))}" if stub else ""
            df["storage_location"] = storage_location
            frames.append(df)
            print(f"{stub.group(0) if stub else folder.name} loaded")

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def import_filter_data_for_wells(analysis_object, expt_label: str, expt_length: int,
                                 delta_t: int, well_list: list) -> pd.DataFrame:
    """Import the summaries for a list of wells, filter them, and tag them.

    `compile_summaries` followed by `filter_summary`, with a `code` column set
    to `expt_label`. See `filter_summary` for what the filters do.
    """
    well_data = compile_summaries(analysis_object, well_list)
    if well_data.empty:
        return well_data
    well_data = filter_summary(well_data, expt_length, delta_t)
    well_data["code"] = expt_label
    return well_data


def import_whole_expt_data(wellmap_dict: dict, analysis_object, expt_length: int,
                           delta_t: int) -> pd.DataFrame:
    """Import and concatenate every well group in a wellmap.

    Iterates the mapping from `create_wellmap_dict` (or, better,
    `wellmap_from_platemap`), calling `import_filter_data_for_wells` per group
    and tagging each with a `code` built from the key.

    Groups with no wells are skipped and a group that raises is reported and
    stepped over, so one bad condition does not lose the rest of the plate.
    """
    frames = []
    for key, wells in (wellmap_dict or {}).items():
        try:
            if not wells:
                print(f"Skipping {key} - no wells")
                continue
            expt_name = "_".join(str(k) for k in key)
            temp_df = import_filter_data_for_wells(
                analysis_object, expt_name, expt_length, delta_t, wells)
            if temp_df.empty:
                print(f"{key} -> no data")
                continue
            frames.append(temp_df)
            print(f"Loaded {key} -> {len(temp_df)} rows")
        except Exception as exc:
            print(f"Failed to load {key}: {type(exc).__name__} {exc}")

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


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
    signal must be background subtracted. A simple method is to subtract the smallest
    signal value from all values.

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
