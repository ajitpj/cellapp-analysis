"""napari browser for the particles listed in a *_summary.xlsx workbook.

    conda run -n img-env python particle_browser.py [top_level_folder]

The top-level folder holds the image stacks (`*_<well>_s<n>_<channel>.tif`) and
one `*_inference` folder per position; each inference folder holds the tracking
table (`*_analysis*.xlsx`) and the per-particle summary (`*_summary*.xlsx`).
Stacks and inference folders are paired on the well_site key (`A12_s2`) rather
than the file stem, because the two can carry different dates.

Selecting a particle pulls a 100x100 ROI that follows its tracked centroid
through every channel, and plots `semantic`, a fluorescence column and
`dead_proba` on a shared 0-1 axis. Particles can be excluded and the surviving
ones written back out as a new summary/analysis pair.

Coordinates: the analysis table's `x`/`y` are centroid row/column on the
half-resolution segmentation grid, so they are doubled for the raw stacks; the
`frame` column indexes the stacks directly (no offset).
"""

from __future__ import annotations

import json
import re
import sys
from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

ROI = 100
COORD_SCALE = 2  # analysis x/y are half-resolution; raw stacks are full
PHASE_OFFSET = 0  # `frame` indexes the raw stacks directly
ROI_CACHE_SIZE = 8

DEFAULT_ROOT = Path(
    "/Volumes/SharedHITSX/cdb-Joglekar-Lab-GL/Anish_Virdi/CycB_dynamics/"
    "ixn/20260624/CycB oe/2026-06-24/20576"
)

# Machine-local and fully regenerable, so it lives outside the repo.
_STATE = Path.home() / ".cache" / "particle_browser"
EXCLUSIONS_FILE = _STATE / "exclusions.json"
TABLE_CACHE = _STATE / "tables"

# `_A12_s2_` in a stack name, an inference folder name or a workbook name.
WELL_SITE = re.compile(r"_([A-H]\d{1,2}_s\d{1,2})(?=[_.])")
# Trailing channel token of a stack name. Underscores are excluded so that
# derived maps (`..._GFP_background_map.tif`) are not mistaken for channels.
CHANNEL = re.compile(r"_[A-H]\d{1,2}_s\d{1,2}_([A-Za-z0-9 ]+)\.tif$")

PLOT_COLORS = {"semantic": "tab:orange", "fluorescence": "tab:green",
               "dead_proba": "tab:red"}
# Preference order for the fluorescence trace; the first one present wins.
FLUOR_PREF = ("GFP", "GFP_bkg_corr", "GFP_int_corr")
# Summary columns worth showing next to the selected particle, when present.
INFO_COLS = ("track_length", "mito_start", "mitosis", "dead_cell_score",
             "fate_label", "death_frame", "n_true_mitotic")


# ---------------------------------------------------------------- data access


@dataclass
class Position:
    """One imaged position: its stacks plus its inference tables."""

    well_site: str
    inference: Path
    summary_path: Path
    analysis_path: Path
    stacks: dict = field(default_factory=dict)

    @property
    def well(self) -> str:
        return self.well_site.split("_s")[0]

    @property
    def site(self) -> str:
        return "s" + self.well_site.split("_s")[1]


def _workbooks(inference: Path):
    """(summary, analysis) paths, preferring the pooled pair when it exists."""
    def usable(pattern):
        # `~$...` are Excel lock files, not data.
        return sorted(p for p in inference.glob(pattern)
                      if not p.name.startswith("~$"))

    summaries = usable("*summary*.xlsx")
    analyses = usable("*analysis*.xlsx")
    if not summaries or not analyses:
        return None, None
    pooled_sum = [p for p in summaries if "pooled" in p.name.lower()]
    summary = (pooled_sum or summaries)[0]
    want_pooled = "pooled" in summary.name.lower()
    matched = [p for p in analyses
               if ("pooled" in p.name.lower()) == want_pooled]
    return summary, (matched or analyses)[0]


def scan(inference_root: Path, image_root: Path) -> dict:
    """well_site -> Position for every inference folder under the roots.

    Both roots are searched one level deep as a fallback, so pointing at a
    parent of the dataset folder still works. Deeper recursion is deliberately
    avoided: these live on an SMB share where rglob is expensive.
    """
    stacks = defaultdict(dict)
    tifs = sorted(image_root.glob("*.tif")) or sorted(image_root.glob("*/*.tif"))
    for p in tifs:
        ws, ch = WELL_SITE.search(p.name), CHANNEL.search(p.name)
        if ws and ch:
            stacks[ws.group(1)][ch.group(1)] = p

    dirs = (sorted(inference_root.glob("*_inference"))
            or sorted(inference_root.glob("*/*_inference")))
    positions = {}
    for inf in dirs:
        m = WELL_SITE.search(inf.name)
        if not m:
            continue
        summary, analysis = _workbooks(inf)
        if summary is None:
            continue
        ws = m.group(1)
        positions[ws] = Position(ws, inf, summary, analysis, dict(stacks.get(ws, {})))
    return positions


def _cached_table(path: Path, **read_kw) -> pd.DataFrame:
    """Read a workbook sheet, memoising it as parquet.

    The analysis tables are ~30 MB of XML; re-parsing one every time a user
    revisits a position costs tens of seconds, while the parquet round-trip is
    well under a second. The stat signature invalidates the cache when the
    workbook is regenerated.
    """
    st = path.stat()
    key = f"{path.stem}_{int(st.st_mtime)}_{st.st_size}"
    cache = TABLE_CACHE / f"{re.sub(r'[^A-Za-z0-9_.-]', '_', key)}.parquet"
    if cache.exists():
        try:
            return pd.read_parquet(cache)
        except Exception:
            pass
    df = pd.read_excel(path, **read_kw)
    try:
        TABLE_CACHE.mkdir(parents=True, exist_ok=True)
        df.to_parquet(cache, index=False)
    except Exception:
        pass  # a cache miss is never fatal
    return df


class Store:
    """Positions plus lazily-loaded tables and the exclusion set."""

    def __init__(self):
        self.positions: dict = {}
        self.inference_root: Path | None = None
        self.image_root: Path | None = None
        self._summaries: dict = {}
        self._analyses: dict = {}
        self.excluded: dict = defaultdict(set)

    # -- discovery -----------------------------------------------------
    def rescan(self, inference_root: Path, image_root: Path) -> None:
        self.positions = scan(Path(inference_root), Path(image_root))
        self.inference_root, self.image_root = Path(inference_root), Path(image_root)
        self._summaries.clear()
        self._analyses.clear()
        self.load_exclusions()

    def wells(self) -> list:
        return sorted({p.well for p in self.positions.values()})

    def sites(self, well: str) -> list:
        return sorted((p.well_site for p in self.positions.values()
                       if p.well == well),
                      key=lambda ws: int(ws.split("_s")[1]))

    # -- tables --------------------------------------------------------
    def summary(self, well_site: str) -> pd.DataFrame:
        if well_site not in self._summaries:
            pos = self.positions[well_site]
            self._summaries[well_site] = _cached_table(
                pos.summary_path, sheet_name="Summary")
        return self._summaries[well_site]

    def analysis(self, well_site: str) -> pd.DataFrame:
        if well_site not in self._analyses:
            pos = self.positions[well_site]
            self._analyses[well_site] = _cached_table(
                pos.analysis_path, sheet_name=0)
        return self._analyses[well_site]

    def particles(self, well_site: str) -> list:
        return [int(p) for p in self.summary(well_site)["particle"].tolist()]

    def summary_row(self, well_site: str, particle: int):
        s = self.summary(well_site)
        hit = s[s["particle"] == particle]
        return None if hit.empty else hit.iloc[0]

    def track(self, well_site: str, particle: int) -> pd.DataFrame:
        a = self.analysis(well_site)
        return a[a["particle"] == particle].sort_values("frame")

    # -- exclusions ----------------------------------------------------
    def _exclusion_key(self) -> str:
        return str(self.inference_root)

    def load_exclusions(self) -> None:
        self.excluded = defaultdict(set)
        if not EXCLUSIONS_FILE.exists():
            return
        try:
            blob = json.loads(EXCLUSIONS_FILE.read_text())
        except (json.JSONDecodeError, OSError):
            return
        for ws, particles in blob.get(self._exclusion_key(), {}).items():
            self.excluded[ws] = {int(p) for p in particles}

    def save_exclusions(self) -> None:
        blob = {}
        if EXCLUSIONS_FILE.exists():
            try:
                blob = json.loads(EXCLUSIONS_FILE.read_text())
            except (json.JSONDecodeError, OSError):
                blob = {}
        blob[self._exclusion_key()] = {
            ws: sorted(int(p) for p in ps)
            for ws, ps in self.excluded.items() if ps
        }
        EXCLUSIONS_FILE.parent.mkdir(parents=True, exist_ok=True)
        EXCLUSIONS_FILE.write_text(json.dumps(blob, indent=1))

    def is_excluded(self, well_site: str, particle: int) -> bool:
        return int(particle) in self.excluded[well_site]

    def toggle(self, well_site: str, particle: int) -> bool:
        s = self.excluded[well_site]
        particle = int(particle)
        if particle in s:
            s.discard(particle)
        else:
            s.add(particle)
        self.save_exclusions()
        return particle in s

    def kept(self, well_site: str) -> list:
        return [p for p in self.particles(well_site)
                if not self.is_excluded(well_site, p)]


# ------------------------------------------------------------------- ROI reads


def roi_movie(path: Path, track: pd.DataFrame, size: int = ROI):
    """(T, size, size) crop centred on the tracked centroid, plus its frames.

    Zero-pads rather than clamping at the image border so the cell stays
    centred. The stacks are uncompressed and contiguous, so we memory-map and
    slice out only the crop: over SMB that fetches ~20 kB per frame instead of
    the whole 8 MB plane. Falls back to per-page reads if mapping is refused.
    """
    half = size // 2
    try:
        stack = tifffile.memmap(path, mode="r")
    except (ValueError, MemoryError, OSError):
        stack = None

    with tifffile.TiffFile(path) as tf:
        series = tf.series[0]
        n_planes, height, width = series.shape[-3:]
        cols = track[["frame", "x", "y"]]
        idx = cols["frame"] + PHASE_OFFSET
        rows = cols[(idx >= 0) & (idx < n_planes)]
        out = np.zeros((len(rows), size, size), dtype=series.dtype)
        for i, rec in enumerate(rows.itertuples(index=False)):
            cr = int(round(rec.x * COORD_SCALE))
            cc = int(round(rec.y * COORD_SCALE))
            r0, c0 = cr - half, cc - half
            rs, re_ = max(r0, 0), min(r0 + size, height)
            cs, ce = max(c0, 0), min(c0 + size, width)
            if rs >= re_ or cs >= ce:
                continue
            f = int(rec.frame) + PHASE_OFFSET
            patch = (stack[f, rs:re_, cs:ce] if stack is not None
                     else tf.pages[f].asarray()[rs:re_, cs:ce])
            out[i, rs - r0:re_ - r0, cs - c0:ce - c0] = patch

    if stack is not None:
        del stack
    return out, rows["frame"].to_numpy()


def channel_order(stacks) -> list:
    """Phase first so it sits at the bottom of the napari layer stack.

    Phase is drawn opaque; if it were added last it would hide the additive
    fluorescence layers entirely.
    """
    return sorted(stacks, key=lambda ch: (ch.lower() not in ("phs", "phase"), ch))


class RoiCache:
    """LRU of per-particle ROI stacks, keyed by (well_site, particle).

    Each entry is a few MB (2 channels x 100x100 x <=450 frames, uint16) but
    costs seconds to re-read over SMB, and users revisit particles constantly.
    """

    def __init__(self, maxsize: int = ROI_CACHE_SIZE):
        self.maxsize = maxsize
        self._d: OrderedDict = OrderedDict()

    def load(self, pos: Position, particle: int, track: pd.DataFrame):
        key = (pos.well_site, int(particle))
        if key in self._d:
            self._d.move_to_end(key)
            return self._d[key]
        movies, frames = {}, np.array([], dtype=int)
        for ch in channel_order(pos.stacks):
            movies[ch], frames = roi_movie(pos.stacks[ch], track)
        self._d[key] = (movies, frames)
        self._d.move_to_end(key)
        while len(self._d) > self.maxsize:
            self._d.popitem(last=False)
        return movies, frames

    def has(self, pos: Position, particle: int) -> bool:
        return (pos.well_site, int(particle)) in self._d

    def clear(self) -> None:
        self._d.clear()


# ----------------------------------------------------------------- trace prep


def rescale(values: np.ndarray, robust: bool = False):
    """Map a trace onto 0-1 for overlaying, plus its original (lo, hi).

    Fluorescence is scaled on 1st-99th percentiles and clipped so that a single
    bright frame does not flatten the rest of the track; the categorical and
    probability traces use their true range.
    """
    v = np.asarray(values, dtype=float)
    finite = v[np.isfinite(v)]
    if finite.size == 0:
        return v, (np.nan, np.nan)
    if robust and finite.size > 4:
        lo, hi = np.percentile(finite, [1, 99])
    else:
        lo, hi = float(finite.min()), float(finite.max())
    if not np.isfinite(hi - lo) or hi - lo <= 0:
        return np.where(np.isfinite(v), 0.5, np.nan), (lo, hi)
    return np.clip((v - lo) / (hi - lo), 0, 1), (lo, hi)


def fluor_columns(track: pd.DataFrame) -> list:
    """Numeric columns that plausibly hold a fluorescence measurement."""
    skip = {"frame", "x", "y", "area", "eccentricity", "label", "particle",
            "semantic", "semantic_smoothed", "dead_flag", "mitotic", "index",
            "offset", "mitotic_proba", "dead_proba", "Unnamed: 0"}
    cols = [c for c in track.columns
            if c not in skip and not str(c).startswith("bbox-")
            and pd.api.types.is_numeric_dtype(track[c])
            and track[c].notna().any()]
    cols.sort(key=lambda c: (c not in FLUOR_PREF,
                             FLUOR_PREF.index(c) if c in FLUOR_PREF else 0))
    return cols


# ------------------------------------------------------------------- exporting


def _sheet(writer, name: str, df: pd.DataFrame) -> None:
    """Write one sheet, restoring an unnamed first column as a written index.

    The pipeline's workbooks store the row index as a headerless first column;
    round-tripping it as a column would leave a stray `Unnamed: 0` header.
    """
    indexed = len(df.columns) and str(df.columns[0]).startswith("Unnamed")
    if indexed:
        df = df.set_index(df.columns[0])
        df.index.name = None
    df.to_excel(writer, sheet_name=name, index=indexed)


def _write_workbook(src: Path, out: Path, replace: dict, extra: dict) -> None:
    """Copy a workbook, swapping in replacement sheets and appending extras."""
    sheets = pd.read_excel(src, sheet_name=None)
    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        for name, df in sheets.items():
            _sheet(writer, name, replace.get(name, df))
        for name, df in extra.items():
            if name not in sheets:
                df.to_excel(writer, sheet_name=name, index=False)


def export_position(store: Store, well_site: str, summary_out: Path,
                    analysis_out: Path) -> tuple:
    """Write the kept particles of one position to a new workbook pair."""
    pos = store.positions[well_site]
    keep = set(store.kept(well_site))
    dropped = sorted(store.excluded[well_site])

    summary = store.summary(well_site)
    kept_summary = summary[summary["particle"].isin(keep)]
    analysis = store.analysis(well_site)
    kept_analysis = analysis[analysis["particle"].isin(keep)]

    note = pd.DataFrame({"excluded_particle": dropped})
    _write_workbook(pos.summary_path, summary_out,
                    {"Summary": kept_summary}, {"excluded": note})

    # The analysis workbook is a single large sheet, so it is rewritten from
    # the cached table rather than re-parsed (~30 MB of XML) just to copy it.
    try:
        with pd.ExcelFile(pos.analysis_path) as xl:
            sheet_name = xl.sheet_names[0]
    except Exception:
        sheet_name = "Sheet1"
    with pd.ExcelWriter(analysis_out, engine="openpyxl") as writer:
        _sheet(writer, sheet_name, kept_analysis)
    return len(kept_summary), len(dropped), len(kept_analysis)


# ------------------------------------------------------------------------ GUI


def build(root: Path):
    import napari
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
    from matplotlib.figure import Figure
    from qtpy.QtWidgets import (
        QCheckBox,
        QComboBox,
        QFileDialog,
        QFrame,
        QHBoxLayout,
        QLabel,
        QLineEdit,
        QMessageBox,
        QPushButton,
        QVBoxLayout,
        QWidget,
    )

    store = Store()
    roi_cache = RoiCache()
    state = {"loading": False, "frames": None, "particle": None}

    viewer = napari.Viewer(title="Particle browser")
    panel = QWidget()
    layout = QVBoxLayout(panel)

    # --- folders --------------------------------------------------------
    inf_box, img_box = QLineEdit(str(root)), QLineEdit(str(root))
    inf_btn = QPushButton("Inference folder...")
    img_btn = QPushButton("Image folder...")
    same_box = QCheckBox("image folder = inference folder")
    same_box.setChecked(True)
    img_box.setEnabled(False)
    img_btn.setEnabled(False)
    for btn, box in ((inf_btn, inf_box), (img_btn, img_box)):
        row = QHBoxLayout()
        row.addWidget(btn)
        row.addWidget(box, 1)
        layout.addLayout(row)
    layout.addWidget(same_box)

    rule0 = QFrame()
    rule0.setFrameShape(QFrame.HLine)
    layout.addWidget(rule0)

    # --- selection ------------------------------------------------------
    well_box, site_box, part_box = QComboBox(), QComboBox(), QComboBox()
    for label, box in (("Well", well_box), ("Site", site_box),
                       ("Particle", part_box)):
        row = QHBoxLayout()
        row.addWidget(QLabel(label))
        row.addWidget(box, 1)
        layout.addLayout(row)

    nav = QHBoxLayout()
    prev_btn, next_btn = QPushButton("< prev"), QPushButton("next >")
    nav.addWidget(prev_btn)
    nav.addWidget(next_btn)
    layout.addLayout(nav)

    fluor_box = QComboBox()
    frow = QHBoxLayout()
    frow.addWidget(QLabel("Fluorescence"))
    frow.addWidget(fluor_box, 1)
    layout.addLayout(frow)

    info = QLabel("")
    info.setWordWrap(True)
    layout.addWidget(info)

    fig = Figure(figsize=(4.2, 2.8), tight_layout=True)
    ax = fig.add_subplot(111)
    canvas = FigureCanvasQTAgg(fig)
    layout.addWidget(canvas)

    exclude_btn = QPushButton("Exclude this particle")
    export_btn = QPushButton("Export selected particles...")
    layout.addWidget(exclude_btn)
    layout.addWidget(export_btn)

    status = QLabel("")
    status.setWordWrap(True)
    layout.addWidget(status)
    layout.addStretch(1)

    # --- helpers --------------------------------------------------------
    def current():
        return site_box.currentData(), part_box.currentData()

    def part_label(well_site, particle):
        row = store.summary_row(well_site, particle)
        bits = [str(particle)]
        if row is not None and "fate_label" in row.index and pd.notna(row["fate_label"]):
            bits.append(str(row["fate_label"]))
        elif row is not None and "track_length" in row.index:
            bits.append(f"len {int(row['track_length'])}")
        text = "  -  ".join(bits)
        return f"[excluded]  {text}" if store.is_excluded(well_site, particle) else text

    def refresh_part_label():
        ws, particle = current()
        if particle is None:
            return
        part_box.setItemText(part_box.currentIndex(), part_label(ws, particle))

    def show_traces(well_site, particle, track, fluor_col):
        """Overlay the rescaled traces; each keeps its native range in the key."""
        ax.clear()
        frames = track["frame"].to_numpy()
        drawn = []

        specs = [("semantic", "semantic", False),
                 (fluor_col, "fluorescence", True),
                 ("dead_proba", "dead_proba", False)]
        for col, kind, robust in specs:
            if col is None or col not in track.columns:
                continue
            values = track[col].to_numpy(dtype=float)
            if not np.isfinite(values).any():
                continue
            scaled, (lo, hi) = rescale(values, robust=robust)
            label = f"{col} [{lo:.3g}, {hi:.3g}]"
            ax.plot(frames, scaled, lw=1.2, color=PLOT_COLORS[kind], label=label)
            drawn.append(col)

        row = store.summary_row(well_site, particle)
        if row is not None:
            for col, color in (("mito_start", "tab:blue"),
                               ("death_frame", "black")):
                if col in row.index and pd.notna(row[col]):
                    ax.axvline(float(row[col]), color=color, ls="--", lw=1.0)
                    ax.annotate(col, (float(row[col]), 1.03),
                                fontsize=6.5, color=color, ha="center")

        cursor = ax.axvline(frames[0] if len(frames) else 0,
                            color="0.4", lw=0.8, alpha=0.8)
        state["cursor"] = cursor
        # Headroom above 1.0 keeps the legend and the event labels clear of the
        # traces, which use the full 0-1 range.
        ax.set_ylim(-0.05, 1.45)
        ax.set_xlabel("frame")
        ax.set_ylabel("rescaled 0-1")
        ax.set_title(f"{well_site}  particle {particle}", fontsize=9)
        if drawn:
            ax.legend(fontsize=5.5, loc="upper center", ncol=len(drawn),
                      framealpha=0.7, borderpad=0.3, columnspacing=0.9,
                      handlelength=1.2)
        canvas.draw_idle()

    def move_cursor(event=None):
        """Keep the trace cursor on the frame napari is showing."""
        frames = state.get("frames")
        cursor = state.get("cursor")
        if frames is None or cursor is None or not len(frames):
            return
        i = int(viewer.dims.current_step[0]) if viewer.dims.ndim else 0
        i = max(0, min(i, len(frames) - 1))
        cursor.set_xdata([frames[i], frames[i]])
        canvas.draw_idle()

    def on_plot_click(event):
        """Clicking the trace jumps the viewer to that frame."""
        frames = state.get("frames")
        if event.inaxes is not ax or frames is None or not len(frames):
            return
        i = int(np.argmin(np.abs(frames - event.xdata)))
        step = list(viewer.dims.current_step)
        step[0] = i
        viewer.dims.current_step = tuple(step)

    def show_info(well_site, particle, track, frames):
        row = store.summary_row(well_site, particle)
        bits = []
        if row is not None:
            for col in INFO_COLS:
                if col in row.index and pd.notna(row[col]):
                    v = row[col]
                    bits.append(f"{col}: {v if isinstance(v, str) else f'{v:g}'}")
        span = f"{frames.min()}-{frames.max()}" if len(frames) else "no frames"
        excl = " <b>[EXCLUDED]</b>" if store.is_excluded(well_site, particle) else ""
        info.setText(f"<b>{well_site} / particle {particle}</b>{excl}<br>"
                     f"{len(track)} rows, frames {span}<br>"
                     + " &nbsp; ".join(bits))

    def load_particle():
        well_site, particle = current()
        if state["loading"] or well_site is None or particle is None:
            return
        pos = store.positions[well_site]
        track = store.track(well_site, particle)
        if track.empty:
            status.setText(f"<b>particle {particle} has no rows in "
                           f"{pos.analysis_path.name}</b>")
            return

        # Keep the fluorescence choice stable across particles when possible.
        cols = fluor_columns(track)
        wanted = fluor_box.currentText()
        fluor_box.blockSignals(True)
        if [fluor_box.itemText(i) for i in range(fluor_box.count())] != cols:
            fluor_box.clear()
            fluor_box.addItems(cols)
        if wanted in cols:
            fluor_box.setCurrentText(wanted)
        fluor_box.blockSignals(False)
        fluor_col = fluor_box.currentText() or (cols[0] if cols else None)

        if not pos.stacks:
            status.setText(f"<b>no image stacks found for {well_site}</b> "
                           f"under {store.image_root}")
            movies, frames = {}, track["frame"].to_numpy()
        else:
            status.setText(
                f"{well_site} / particle {particle} - "
                + ("from cache..." if roi_cache.has(pos, particle)
                   else f"reading ROI from {len(pos.stacks)} stack(s)...")
            )
            panel.repaint()
            movies, frames = roi_cache.load(pos, particle, track)

        state["frames"], state["particle"] = frames, particle
        for ch, movie in movies.items():
            colormap = {"GFP": "green", "Texas Red": "magenta"}.get(ch, "gray")
            lo, hi = np.percentile(movie, [1, 99.5]) if movie.size else (0, 1)
            if ch in viewer.layers:
                layer = viewer.layers[ch]
                layer.data = movie
            else:
                layer = viewer.add_image(
                    movie, name=ch, colormap=colormap,
                    blending="translucent" if ch == "phs" else "additive")
            layer.contrast_limits = (float(lo), float(max(hi, lo + 1)))
        for layer in list(viewer.layers):
            if layer.name not in movies:
                viewer.layers.remove(layer)

        show_traces(well_site, particle, track, fluor_col)
        show_info(well_site, particle, track, frames)
        move_cursor()
        exclude_btn.setText("Include this particle"
                            if store.is_excluded(well_site, particle)
                            else "Exclude this particle")
        n_excl = len(store.excluded[well_site])
        status.setText(f"{well_site}: {len(store.particles(well_site))} particles, "
                       f"{n_excl} excluded. Channels: "
                       f"{', '.join(sorted(pos.stacks)) or 'none'}")

    def refill_particles():
        well_site = site_box.currentData()
        part_box.blockSignals(True)
        part_box.clear()
        if well_site:
            try:
                particles = store.particles(well_site)
            except Exception as exc:  # unreadable workbook
                particles = []
                status.setText(f"<b>could not read summary for {well_site}: {exc}</b>")
            for p in particles:
                part_box.addItem(part_label(well_site, p), p)
        part_box.blockSignals(False)
        load_particle()

    def refill_sites():
        well = well_box.currentData()
        site_box.blockSignals(True)
        site_box.clear()
        for ws in store.sites(well) if well else []:
            pos = store.positions[ws]
            text = pos.site if pos.stacks else f"{pos.site}  (no images)"
            site_box.addItem(text, ws)
        site_box.blockSignals(False)
        refill_particles()

    def rescan():
        inference_root = Path(inf_box.text().strip())
        image_root = (inference_root if same_box.isChecked()
                      else Path(img_box.text().strip()))
        if not inference_root.is_dir():
            status.setText(f"<b>not a folder: {inference_root}</b>")
            return
        status.setText(f"scanning {inference_root}...")
        panel.repaint()
        roi_cache.clear()
        store.rescan(inference_root, image_root)
        state["loading"] = True
        well_box.clear()
        for w in store.wells():
            well_box.addItem(w, w)
        state["loading"] = False
        if not store.positions:
            status.setText(f"<b>no *_inference folders under {inference_root}</b>")
            for box in (site_box, part_box):
                box.clear()
            return
        n_img = sum(1 for p in store.positions.values() if p.stacks)
        status.setText(f"{len(store.positions)} positions "
                       f"({n_img} with image stacks) in {inference_root.name}")
        refill_sites()

    def pick_folder(box, title):
        d = QFileDialog.getExistingDirectory(panel, title, box.text().strip())
        if not d:
            return
        box.setText(d)
        if box is inf_box and same_box.isChecked():
            img_box.setText(d)
        rescan()

    def on_same_toggled(checked):
        img_box.setEnabled(not checked)
        img_btn.setEnabled(not checked)
        if checked:
            img_box.setText(inf_box.text())
        rescan()

    def on_exclude():
        well_site, particle = current()
        if particle is None:
            return
        now_excluded = store.toggle(well_site, particle)
        exclude_btn.setText("Include this particle" if now_excluded
                            else "Exclude this particle")
        refresh_part_label()
        track = store.track(well_site, particle)
        show_info(well_site, particle, track, state["frames"]
                  if state["frames"] is not None else track["frame"].to_numpy())
        status.setText(f"particle {particle} "
                       f"{'excluded from' if now_excluded else 'restored to'} "
                       f"{well_site} ({len(store.kept(well_site))} kept)")

    def on_export():
        well_site, _ = current()
        if well_site is None:
            return
        touched = [ws for ws, ps in store.excluded.items() if ps]
        targets = [well_site]
        others = [ws for ws in touched if ws != well_site]
        if others:
            answer = QMessageBox.question(
                panel, "Export scope",
                f"{len(others)} other position(s) also have exclusions "
                f"({', '.join(sorted(others))}).\n\n"
                "Yes: export every position with exclusions.\n"
                f"No: export only {well_site}.",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if answer == QMessageBox.Yes:
                targets = sorted(set(touched) | {well_site})

        suggested = (store.positions[well_site].inference
                     / f"{well_site}_curated.xlsx")
        chosen, _ = QFileDialog.getSaveFileName(
            panel, "Export curated summary + analysis as...",
            str(suggested), "Excel workbook (*.xlsx)")
        if not chosen:
            return
        chosen = Path(chosen)
        stem = re.sub(r"(_summary|_analysis)$", "", chosen.stem)

        written = []
        for ws in targets:
            tag = f"_{ws}" if len(targets) > 1 else ""
            summary_out = chosen.with_name(f"{stem}{tag}_summary.xlsx")
            analysis_out = chosen.with_name(f"{stem}{tag}_analysis.xlsx")
            status.setText(f"writing {summary_out.name} + {analysis_out.name}...")
            panel.repaint()
            try:
                n_kept, n_drop, n_rows = export_position(
                    store, ws, summary_out, analysis_out)
            except Exception as exc:
                QMessageBox.critical(panel, "Export failed", f"{ws}: {exc}")
                status.setText(f"<b>export failed for {ws}: {exc}</b>")
                return
            written.append(f"{ws}: {n_kept} particles kept ({n_drop} excluded), "
                           f"{n_rows} tracking rows")
        QMessageBox.information(
            panel, "Export complete",
            f"Wrote {2 * len(written)} files to {chosen.parent}\n\n"
            + "\n".join(written))
        status.setText("exported - " + "; ".join(written))

    def step_particle(delta):
        i = part_box.currentIndex() + delta
        if 0 <= i < part_box.count():
            part_box.setCurrentIndex(i)

    # --- wiring ---------------------------------------------------------
    inf_btn.clicked.connect(lambda: pick_folder(inf_box, "Inference folder"))
    img_btn.clicked.connect(lambda: pick_folder(img_box, "Image folder"))
    inf_box.editingFinished.connect(rescan)
    img_box.editingFinished.connect(rescan)
    same_box.toggled.connect(on_same_toggled)
    well_box.currentIndexChanged.connect(lambda _: refill_sites())
    site_box.currentIndexChanged.connect(lambda _: refill_particles())
    part_box.currentIndexChanged.connect(lambda _: load_particle())
    fluor_box.currentIndexChanged.connect(lambda _: load_particle())
    prev_btn.clicked.connect(lambda: step_particle(-1))
    next_btn.clicked.connect(lambda: step_particle(1))
    exclude_btn.clicked.connect(on_exclude)
    export_btn.clicked.connect(on_export)
    viewer.dims.events.current_step.connect(move_cursor)
    canvas.mpl_connect("button_press_event", on_plot_click)

    viewer.window.add_dock_widget(panel, name="Particles", area="right")
    rescan()
    return viewer


def main() -> None:
    import napari

    root = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_ROOT
    build(root)
    napari.run()


if __name__ == "__main__":
    main()
