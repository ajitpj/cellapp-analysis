"""napari browser for the particles listed in a *_summary.xlsx workbook.

    conda run -n img-env python particle_browser.py [top_level_folder]

The top-level folder holds the image stacks (`*_<well>_s<n>_<channel>.tif`) and
one `*_inference` folder per position; each inference folder holds the tracking
table (`*_analysis*.xlsx`) and the per-particle summary (`*_summary*.xlsx`).
Stacks and inference folders are paired on the well_site key (`A12_s2`) rather
than the file stem, because the two can carry different dates.

Selecting a particle pulls a 100x100 ROI that follows its tracked centroid
through every channel, and plots `semantic`, one or all of the fluorescence
columns, and `dead_proba` on a shared 0-1 axis. Columns a legacy workbook does
not carry are left off the plot rather than drawn as zeros. Particles can be
excluded and annotated, and the surviving ones written back out as a new
summary/analysis pair carrying the annotations.

Coordinates: the analysis table's `x`/`y` are centroid row/column on the
half-resolution segmentation grid, so they are doubled for the raw stacks; the
`frame` column indexes the stacks directly (no offset).
"""

from __future__ import annotations

import hashlib
import json
import queue
import re
import sys
import threading
from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

ROI = 100
COORD_SCALE = 2  # analysis x/y are half-resolution; raw stacks are full
PHASE_OFFSET = 0  # `frame` indexes the raw stacks directly

# ROI stacks are held by total size rather than by count, since a track can be
# anything from a handful of frames to the whole movie. 1.5 GB is roughly 80
# two-channel 450-frame particles - enough to re-visit a whole position without
# touching the share again, and small enough beside a workstation's RAM.
ROI_CACHE_BYTES = 1_500_000_000
ROI_CACHE_MIN = 4  # never evict below this, however large the entries are

# Below this many particles a bulk pass is not worth it. Reading a whole plane
# costs ~17 crops locally and ~2-3 over SMB, where a crop already drags in
# ~400 KB of pages across 100 separate rows; past a handful of particles one
# sequential pass beats per-particle reads by two orders of magnitude.
BULK_MIN_PARTICLES = 4

# Home, so the browser opens somewhere that exists on any machine; scan() only
# looks one level deep, so this costs a couple of directory listings and finds
# nothing until the user points the folder boxes at a dataset.
DEFAULT_ROOT = Path.home()

# Machine-local and fully regenerable, so it lives outside the repo.
_STATE = Path.home() / ".cache" / "particle_browser"
STATE_FILE = _STATE / "exclusions.json"
TABLE_CACHE = _STATE / "tables"
ROI_CACHE_DIR = _STATE / "roi"

# `_A12_s2_` in a stack name, an inference folder name or a workbook name.
WELL_SITE = re.compile(r"_([A-H]\d{1,2}_s\d{1,2})(?=[_.])")
# Trailing channel token of a stack name. Underscores are excluded so that
# derived maps (`..._GFP_background_map.tif`) are not mistaken for channels.
CHANNEL = re.compile(r"_[A-H]\d{1,2}_s\d{1,2}_([A-Za-z0-9 ]+)\.tif$")

SEMANTIC_COLOR = "tab:orange"
DEAD_COLOR = "tab:red"
# Cycled when several fluorescence traces are drawn at once. Orange and red are
# reserved for semantic and dead_proba.
FLUOR_COLORS = ("tab:green", "tab:blue", "tab:purple", "tab:olive",
                "tab:cyan", "tab:brown", "tab:pink")
# Preference order for the fluorescence trace; the first one present wins.
FLUOR_PREF = ("GFP", "GFP_bkg_corr", "GFP_int_corr")
# Shown next to the selected particle when present. Both the current names and
# the ones they replaced are listed, so workbooks written before the rename
# still display.
INFO_COLS = ("track_length", "n_peaks", "mitotic_start_frame",
             "frames_to_death", "sem_frames_in_mitosis",
             "corrected_frames_in_mitosis", "dead_cell_score", "fate_label",
             "death_frame",
             "mito_start", "mitosis", "time_to_death", "n_sem_mitotic")
ANNOTATION_COL = "user_annotation"

# Grouping. The pipeline's own per-particle call is `fate_label`
# (mitotic_survived / dead_in_mitosis / dead_no_mitosis / dead_post_mitosis),
# so it leads the list; the user's own notes are offered as a second axis.
NOTE_GROUP = "user annotation"
NO_GROUP = "(none)"
NO_NOTE = "(no note)"
BLANK_GROUP = "(unlabelled)"
MAX_GROUPS = 40  # a column with more distinct values than this is not a grouping


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
        self.notes: dict = defaultdict(dict)

    # -- discovery -----------------------------------------------------
    def rescan(self, inference_root: Path, image_root: Path) -> None:
        self.positions = scan(Path(inference_root), Path(image_root))
        self.inference_root, self.image_root = Path(inference_root), Path(image_root)
        self._summaries.clear()
        self._analyses.clear()
        self.load_state()

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
            df = _cached_table(pos.summary_path, sheet_name="Summary")
            self._summaries[well_site] = df
            self._seed_notes(well_site, df)
        return self._summaries[well_site]

    def _seed_notes(self, well_site: str, summary: pd.DataFrame) -> None:
        """Adopt annotations already in the workbook, so re-curating an
        exported file picks up where the last pass left off. Notes held in the
        local state file win, since they are the more recent edit."""
        if ANNOTATION_COL not in summary.columns:
            return
        held = self.notes[well_site]
        for particle, note in zip(summary["particle"], summary[ANNOTATION_COL]):
            text = "" if pd.isna(note) else str(note).strip()
            if text and int(particle) not in held:
                held[int(particle)] = text

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

    # -- grouping ------------------------------------------------------
    def group_columns(self, well_site: str) -> list:
        """Summary columns that sort the particles into a few named groups.

        Categorical columns only: a float per particle is a measurement, not a
        group. `fate_label` leads when present since it is the pipeline's own
        verdict on each particle. Legacy summaries carry no such column, which
        is why the caller must cope with an empty list.
        """
        s = self.summary(well_site)
        cols = []
        for c in s.columns:
            if c in ("particle", ANNOTATION_COL):
                continue
            if not (s[c].dtype == object or pd.api.types.is_bool_dtype(s[c])):
                continue
            if 1 <= s[c].nunique(dropna=True) <= MAX_GROUPS:
                cols.append(c)
        cols.sort(key=lambda c: (c != "fate_label", str(c)))
        return cols

    def group_of(self, well_site: str, particle: int, group_by: str) -> str:
        """Which group a particle falls in, as displayed."""
        if group_by == NOTE_GROUP:
            return self.note(well_site, particle) or NO_NOTE
        row = self.summary_row(well_site, particle)
        value = None if row is None else row.get(group_by)
        return BLANK_GROUP if value is None or pd.isna(value) else str(value)

    def groups(self, well_site: str, group_by: str) -> list:
        """(group, count) for one position, largest group first."""
        counts: dict = {}
        for p in self.particles(well_site):
            g = self.group_of(well_site, p, group_by)
            counts[g] = counts.get(g, 0) + 1
        return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))

    def particles_in_group(self, well_site: str, group_by: str, group) -> list:
        """Particles of one group, in summary order. `group` None means all."""
        particles = self.particles(well_site)
        if not group_by or group_by == NO_GROUP or group is None:
            return particles
        return [p for p in particles
                if self.group_of(well_site, p, group_by) == group]

    def track(self, well_site: str, particle: int) -> pd.DataFrame:
        a = self.analysis(well_site)
        return a[a["particle"] == particle].sort_values("frame")

    # -- curation state ------------------------------------------------
    def _state_key(self) -> str:
        return str(self.inference_root)

    def load_state(self) -> None:
        """Read this root's exclusions and annotations from the state file."""
        self.excluded, self.notes = defaultdict(set), defaultdict(dict)
        if not STATE_FILE.exists():
            return
        try:
            blob = json.loads(STATE_FILE.read_text())
        except (json.JSONDecodeError, OSError):
            return
        mine = blob.get(self._state_key(), {})
        # Files written before annotations existed map well_site straight to a
        # list of excluded particles.
        if mine and all(isinstance(v, list) for v in mine.values()):
            mine = {"excluded": mine, "notes": {}}
        for ws, particles in mine.get("excluded", {}).items():
            self.excluded[ws] = {int(p) for p in particles}
        for ws, notes in mine.get("notes", {}).items():
            self.notes[ws] = {int(p): str(t) for p, t in notes.items() if t}

    def save_state(self) -> None:
        blob = {}
        if STATE_FILE.exists():
            try:
                blob = json.loads(STATE_FILE.read_text())
            except (json.JSONDecodeError, OSError):
                blob = {}
        blob[self._state_key()] = {
            "excluded": {ws: sorted(int(p) for p in ps)
                         for ws, ps in self.excluded.items() if ps},
            "notes": {ws: {str(p): t for p, t in sorted(notes.items()) if t}
                      for ws, notes in self.notes.items() if any(notes.values())},
        }
        STATE_FILE.parent.mkdir(parents=True, exist_ok=True)
        STATE_FILE.write_text(json.dumps(blob, indent=1))

    def is_excluded(self, well_site: str, particle: int) -> bool:
        return int(particle) in self.excluded[well_site]

    def toggle(self, well_site: str, particle: int) -> bool:
        s = self.excluded[well_site]
        particle = int(particle)
        if particle in s:
            s.discard(particle)
        else:
            s.add(particle)
        self.save_state()
        return particle in s

    def note(self, well_site: str, particle: int) -> str:
        return self.notes[well_site].get(int(particle), "")

    def set_note(self, well_site: str, particle: int, text: str) -> None:
        """In-memory only; callers persist with save_state at a natural pause
        so that a note is not written to disk on every keystroke."""
        text = text.strip()
        if text:
            self.notes[well_site][int(particle)] = text
        else:
            self.notes[well_site].pop(int(particle), None)

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


def roi_movies_bulk(paths: dict, tracks: dict, size: int = ROI,
                    on_ready=None, should_stop=None) -> None:
    """One sequential pass over the stacks, scattering crops to many particles.

    The cost of a pass is dominated by reading planes, and that is independent
    of how many particles are extracted from them: a 100x100 crop already drags
    in ~400 KB of pages (100 rows, 16 KB pages) and costs one round trip per
    frame, so per-particle reads pay nearly plane price and pay it again for
    every particle. Reading each plane once and cutting every particle's crop
    out of it collapses N per-particle passes into one.

    All channels are advanced together, frame by frame, so a particle is
    complete - every channel, every frame - the moment its last frame is read,
    and can be handed over while the rest of the pass continues. Doing a whole
    channel at a time instead would leave the caller with nothing until the
    final channel finished.

    on_ready(particle, movies, frames) is called on the calling thread as each
    particle completes; should_stop() is polled per frame so a position change
    can abandon the pass promptly.
    """
    half = size // 2
    chans = channel_order(paths)
    maps, files, shapes = {}, {}, {}
    try:
        for ch in chans:
            try:
                maps[ch] = tifffile.memmap(paths[ch], mode="r")
            except (ValueError, MemoryError, OSError):
                maps[ch] = None
            files[ch] = tifffile.TiffFile(paths[ch])
            shapes[ch] = files[ch].series[0].shape[-3:]

        # Plan the pass: which particles want a crop from which plane.
        n_planes = min(s[0] for s in shapes.values())
        out, kept, per_frame = {}, {}, defaultdict(list)
        for pid, track in tracks.items():
            cols = track[["frame", "x", "y"]]
            idx = cols["frame"] + PHASE_OFFSET
            rows = cols[(idx >= 0) & (idx < n_planes)]
            dtype = files[chans[0]].series[0].dtype
            out[pid] = {ch: np.zeros((len(rows), size, size), dtype=dtype)
                        for ch in chans}
            kept[pid] = rows["frame"].to_numpy()
            for i, rec in enumerate(rows.itertuples(index=False)):
                per_frame[int(rec.frame) + PHASE_OFFSET].append(
                    (pid, i, int(round(rec.x * COORD_SCALE)),
                     int(round(rec.y * COORD_SCALE))))
        left = {pid: len(kept[pid]) for pid in out}

        for f in sorted(per_frame):
            if should_stop is not None and should_stop():
                return
            planes = {}
            for ch in chans:
                m = maps[ch]
                planes[ch] = (m[f] if m is not None
                              else files[ch].pages[f].asarray())
            for pid, i, cr, cc in per_frame[f]:
                r0, c0 = cr - half, cc - half
                for ch in chans:
                    height, width = shapes[ch][1:]
                    rs, re_ = max(r0, 0), min(r0 + size, height)
                    cs, ce = max(c0, 0), min(c0 + size, width)
                    if rs < re_ and cs < ce:
                        out[pid][ch][i, rs - r0:re_ - r0, cs - c0:ce - c0] = \
                            planes[ch][rs:re_, cs:ce]
                left[pid] -= 1
                if left[pid] == 0 and on_ready is not None:
                    on_ready(pid, out.pop(pid), kept[pid])
    finally:
        for ch in list(maps):
            maps[ch] = None
        for tf in files.values():
            tf.close()


def _stack_signature(paths: dict) -> str:
    """Identify a position's stacks, so a re-export invalidates its ROIs."""
    parts = []
    for ch in sorted(paths):
        st = paths[ch].stat()
        parts.append(f"{ch}:{paths[ch].name}:{st.st_size}:{int(st.st_mtime)}")
    return hashlib.sha1("|".join(parts).encode()).hexdigest()[:16]


def _roi_disk_path(pos: Position, particle: int, size: int) -> Path:
    return (ROI_CACHE_DIR / pos.well_site /
            f"{int(particle)}_{_stack_signature(pos.stacks)}_{size}.npz")


def roi_disk_load(pos: Position, particle: int, size: int = ROI):
    """Read a cached ROI set, or None. Never raises: the cache is disposable."""
    path = _roi_disk_path(pos, particle, size)
    if not path.exists():
        return None
    try:
        with np.load(path) as z:
            names = [str(n) for n in z["channels"]]
            return {ch: z[f"a{i}"] for i, ch in enumerate(names)}, z["frames"]
    except Exception:
        return None


def roi_disk_save(pos: Position, particle: int, movies: dict, frames,
                  size: int = ROI) -> None:
    """Persist one particle's ROI set. Channel names go in their own array -
    they contain spaces ('Texas Red'), which npz keyword names cannot."""
    path = _roi_disk_path(pos, particle, size)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        names = list(movies)
        tmp = path.with_suffix(".tmp.npz")
        # Plain unicode, not object dtype: an object array would need
        # allow_pickle on the way back in, which np.load refuses by default.
        np.savez(tmp, frames=np.asarray(frames),
                 channels=np.array(names, dtype=str),
                 **{f"a{i}": movies[ch] for i, ch in enumerate(names)})
        tmp.replace(path)
    except Exception:
        pass  # a full or read-only cache dir must not break browsing


def channel_order(stacks) -> list:
    """Phase first so it sits at the bottom of the napari layer stack.

    Phase is drawn opaque; if it were added last it would hide the additive
    fluorescence layers entirely.
    """
    return sorted(stacks, key=lambda ch: (ch.lower() not in ("phs", "phase"), ch))


class RoiCache:
    """LRU of per-particle ROI stacks, keyed by (well_site, particle).

    Bounded by total bytes rather than by entry count: a track can be 20 frames
    or the whole 450-frame movie, so a count-based bound either wastes memory
    or evicts far too eagerly. Re-reading one particle costs ~14 s per channel
    over SMB and users revisit constantly, so the budget is deliberately
    generous -- it is the difference between an interactive browser and a
    slideshow.
    """

    def __init__(self, budget: int = ROI_CACHE_BYTES, minsize: int = ROI_CACHE_MIN):
        self.budget = budget
        self.minsize = minsize
        self._d: OrderedDict = OrderedDict()
        self._bytes = 0

    @staticmethod
    def _size(movies: dict) -> int:
        return sum(int(m.nbytes) for m in movies.values())

    def load(self, pos: Position, particle: int, track: pd.DataFrame):
        key = (pos.well_site, int(particle))
        if key in self._d:
            self._d.move_to_end(key)
            return self._d[key]
        hit = roi_disk_load(pos, particle)
        if hit is not None:
            movies, frames = hit
        else:
            movies, frames = {}, np.array([], dtype=int)
            for ch in channel_order(pos.stacks):
                movies[ch], frames = roi_movie(pos.stacks[ch], track)
            roi_disk_save(pos, particle, movies, frames)
        self._insert(key, movies, frames, viewed=True)
        return movies, frames

    def _insert(self, key, movies, frames, viewed: bool) -> None:
        if key in self._d:
            return
        self._d[key] = (movies, frames)
        self._bytes += self._size(movies)
        # A particle the user actually opened outranks one fetched on spec, so
        # speculative entries go in at the evict-first end of the LRU.
        self._d.move_to_end(key, last=viewed)
        while len(self._d) > self.minsize and self._bytes > self.budget:
            _, (old, _) = self._d.popitem(last=False)
            self._bytes -= self._size(old)

    def put(self, well_site: str, particle: int, movies: dict, frames) -> None:
        """Add a speculatively fetched entry, from the GUI thread only."""
        self._insert((well_site, int(particle)), movies, frames, viewed=False)

    def has(self, pos: Position, particle: int) -> bool:
        return (pos.well_site, int(particle)) in self._d

    def ready(self, pos: Position, particle: int) -> bool:
        """In memory, or on disk and so effectively instant."""
        return (self.has(pos, particle)
                or _roi_disk_path(pos, particle, ROI).exists())

    def stats(self) -> str:
        return (f"{len(self._d)} ROI cached, "
                f"{self._bytes / 1e6:.0f}/{self.budget / 1e6:.0f} MB")

    def clear(self) -> None:
        self._d.clear()
        self._bytes = 0


class RoiPrefetcher:
    """Fills the ROI cache for a group of particles, one pass, off the GUI thread.

    The worker only reads and writes the disk cache; finished particles go on a
    queue and the GUI thread does the in-memory insert, so RoiCache stays
    single-threaded and needs no lock. Every job carries a generation, bumped
    whenever the position or group changes, and results from a superseded
    generation are dropped on arrival rather than cancelled mid-read.
    """

    def __init__(self, cache: RoiCache):
        self.cache = cache
        self._q: queue.Queue = queue.Queue()
        self._gen = 0
        self._thread: threading.Thread | None = None
        self.done = 0
        self.total = 0

    def cancel(self) -> None:
        self._gen += 1
        self.done = self.total = 0

    def start(self, pos: Position, tracks: dict) -> None:
        """Queue a pass over `tracks` ({particle: dataframe}) for `pos`."""
        self.cancel()
        tracks = {p: t for p, t in tracks.items()
                  if not self.cache.ready(pos, p) and not t.empty}
        if not pos.stacks or len(tracks) < BULK_MIN_PARTICLES:
            return
        gen = self._gen
        self.total = len(tracks)

        def ready(pid, movies, frames):
            roi_disk_save(pos, pid, movies, frames)
            self._q.put((gen, pos.well_site, pid, movies, frames))

        def run():
            try:
                roi_movies_bulk(pos.stacks, tracks, on_ready=ready,
                                should_stop=lambda: gen != self._gen)
            except Exception as exc:            # a bad stack must not kill the UI
                self._q.put((gen, None, None, None, exc))

        self._thread = threading.Thread(target=run, daemon=True)
        self._thread.start()

    def drain(self, limit: int = 8) -> bool:
        """Move finished particles into the cache. GUI thread only.

        Returns True if anything was taken, so the caller can refresh a status
        line without polling the queue itself.
        """
        took = False
        for _ in range(limit):
            try:
                gen, well_site, pid, movies, frames = self._q.get_nowait()
            except queue.Empty:
                break
            if gen != self._gen or well_site is None:
                continue                        # superseded, or a worker error
            self.cache.put(well_site, pid, movies, frames)
            self.done += 1
            took = True
        return took

    @property
    def running(self) -> bool:
        return self.total > 0 and self.done < self.total


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
    """Copy a workbook, swapping in replacement sheets and appending extras.

    An extra whose name is already in the source overwrites it rather than
    being dropped, so re-curating a previous export refreshes its `excluded`
    sheet instead of carrying the stale one forward.
    """
    sheets = pd.read_excel(src, sheet_name=None)
    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        for name, df in sheets.items():
            if name in extra:
                extra[name].to_excel(writer, sheet_name=name, index=False)
            else:
                _sheet(writer, name, replace.get(name, df))
        for name, df in extra.items():
            if name not in sheets:
                df.to_excel(writer, sheet_name=name, index=False)


def _annotate(df: pd.DataFrame, notes: dict) -> pd.DataFrame:
    """Attach the user annotation of each row's particle as a column."""
    out = df.copy()
    out[ANNOTATION_COL] = (out["particle"].map(lambda p: notes.get(int(p), ""))
                           .fillna(""))
    return out


def export_position(store: Store, well_site: str, summary_out: Path,
                    analysis_out: Path) -> tuple:
    """Write the kept particles of one position to a new workbook pair."""
    pos = store.positions[well_site]
    keep = set(store.kept(well_site))
    dropped = sorted(store.excluded[well_site])
    notes = store.notes[well_site]

    summary = store.summary(well_site)
    kept_summary = _annotate(summary[summary["particle"].isin(keep)], notes)
    analysis = store.analysis(well_site)
    kept_analysis = analysis[analysis["particle"].isin(keep)]
    # Only widen the per-frame table when there is something to carry; an empty
    # column across 50k rows is pure noise for downstream readers.
    if any(notes.get(int(p)) for p in keep):
        kept_analysis = _annotate(kept_analysis, notes)

    note = pd.DataFrame({"excluded_particle": dropped,
                         ANNOTATION_COL: [notes.get(int(p), "") for p in dropped]})
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
    from qtpy.QtCore import QTimer
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
    prefetch = RoiPrefetcher(roi_cache)
    state = {"loading": False, "frames": None, "particle": None,
             "note_owner": None}

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
    well_box, site_box = QComboBox(), QComboBox()
    group_by_box, group_box, part_box = QComboBox(), QComboBox(), QComboBox()
    group_by_box.setToolTip("Sort this position's particles into groups - the "
                            "pipeline's fate_label, or your own annotations")
    group_box.setToolTip("Review only the particles in one group")
    for label, box in (("Well", well_box), ("Site", site_box),
                       ("Group by", group_by_box), ("Group", group_box),
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
    fluor_all = QCheckBox("plot all")
    fluor_all.setToolTip("Overlay every fluorescence column found in the "
                         "analysis table instead of just the selected one")
    frow = QHBoxLayout()
    frow.addWidget(QLabel("Fluorescence"))
    frow.addWidget(fluor_box, 1)
    frow.addWidget(fluor_all)
    layout.addLayout(frow)

    annot_box = QLineEdit()
    annot_box.setPlaceholderText("free-text note for this particle")
    arow = QHBoxLayout()
    arow.addWidget(QLabel("Annotation"))
    arow.addWidget(annot_box, 1)
    layout.addLayout(arow)

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
    # Its own line: load_particle rewrites `status` on every particle, and the
    # background pass reports on its own schedule.
    cache_label = QLabel("")
    cache_label.setWordWrap(True)
    layout.addWidget(cache_label)
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
        note = store.note(well_site, particle)
        if note:
            bits.append(note if len(note) <= 30 else note[:29] + "...")
        text = "  -  ".join(bits)
        marks = ("[excluded] " if store.is_excluded(well_site, particle) else
                 "") + ("[note] " if note else "")
        return f"{marks} {text}" if marks else text

    def refresh_part_label():
        ws, particle = current()
        if particle is None:
            return
        part_box.setItemText(part_box.currentIndex(), part_label(ws, particle))

    def on_note_edited(text):
        """Keep the note with the particle it was typed against.

        The edit is held in memory on every keystroke and only written to disk
        at a pause, so nothing is lost if the user switches particles mid-word.
        """
        owner = state.get("note_owner")
        if owner is None:
            return
        store.set_note(owner[0], owner[1], text)
        refresh_part_label()

    def flush_note():
        """Persist the note being edited, if any."""
        if state.get("note_owner") is not None:
            store.set_note(*state["note_owner"], annot_box.text())
            store.save_state()

    def commit_note():
        """Persist an edited note, and when the list is grouped by the notes
        themselves, rebuild the groups around the particle in view."""
        owner = state.get("note_owner")
        if owner is None:
            return
        flush_note()
        if group_by_box.currentData() != NOTE_GROUP:
            return
        well_site, particle = owner
        if (group_box.currentData() is not None
                and store.group_of(well_site, particle, NOTE_GROUP)
                != group_box.currentData()):
            # The particle has just left the group under review; widen to all
            # rather than dropping it out from under the user.
            group_box.blockSignals(True)
            group_box.setCurrentIndex(0)
            group_box.blockSignals(False)
        # Deferred so that the click which moved focus out of the box is
        # handled before the lists it may be landing on are rebuilt.
        QTimer.singleShot(0, lambda: refill_groups(keep=particle))

    def show_traces(well_site, particle, track, fluor_cols):
        """Overlay the rescaled traces; each keeps its native range in the key.

        Legacy workbooks predate the dead classifier, so `dead_proba` and
        `death_frame` are simply absent (or all-NaN, which means the model was
        never asked about those frames). Either way they are left off the plot
        rather than drawn as a flat line at zero.
        """
        ax.clear()
        frames = track["frame"].to_numpy()
        drawn = []

        specs = [("semantic", SEMANTIC_COLOR, False)]
        specs += [(col, FLUOR_COLORS[i % len(FLUOR_COLORS)], True)
                  for i, col in enumerate(fluor_cols)]
        specs += [("dead_proba", DEAD_COLOR, False)]
        for col, color, robust in specs:
            if col is None or col not in track.columns:
                continue
            values = track[col].to_numpy(dtype=float)
            if not np.isfinite(values).any():
                continue
            scaled, (lo, hi) = rescale(values, robust=robust)
            label = f"{col} [{lo:.3g}, {hi:.3g}]"
            # A correction factor that varies by <1% of its own magnitude is
            # noise once stretched over the full axis; say so rather than
            # letting it look like signal next to the real traces.
            mid = np.nanmedian(values)
            flat = (robust and np.isfinite(mid) and mid != 0
                    and (hi - lo) / abs(mid) < 0.01)
            ax.plot(frames, scaled, lw=1.2, color=color,
                    alpha=0.5 if flat else 1.0, ls=":" if flat else "-",
                    label=label + (" ~flat" if flat else ""))
            drawn.append(col)

        row = store.summary_row(well_site, particle)
        if row is not None:
            starts = [c for c in ("mitotic_start_frame", "mito_start")
                      if c in row.index]
            for col, color in ((starts[0] if starts else None, "tab:blue"),
                               ("death_frame", "black")):
                if col is not None and col in row.index and pd.notna(row[col]):
                    ax.axvline(float(row[col]), color=color, ls="--", lw=1.0)
                    ax.annotate(col, (float(row[col]), 1.03),
                                fontsize=6.5, color=color, ha="center")

        cursor = ax.axvline(frames[0] if len(frames) else 0,
                            color="0.4", lw=0.8, alpha=0.8)
        state["cursor"] = cursor
        # Headroom above 1.0 keeps the legend and the event labels clear of the
        # traces, which use the full 0-1 range. A tall legend needs more.
        ncol = 2 if len(drawn) > 3 else max(len(drawn), 1)
        rows = int(np.ceil(len(drawn) / ncol))
        ax.set_ylim(-0.05, 1.12 + 0.14 * rows)
        ax.set_xlabel("frame")
        ax.set_ylabel("rescaled 0-1")
        ax.set_title(f"{well_site}  particle {particle}", fontsize=9)
        if drawn:
            ax.legend(fontsize=5.5, loc="upper center", ncol=ncol,
                      framealpha=0.7, borderpad=0.3, columnspacing=0.9,
                      handlelength=1.2)
        canvas.draw_idle()
        return drawn

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
        flush_note()  # the box still holds the note of the previous particle
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
        fluor_all.setEnabled(len(cols) > 1)
        if len(cols) < 2:
            fluor_all.setChecked(False)
        fluor_box.setEnabled(not fluor_all.isChecked())
        if fluor_all.isChecked():
            fluor_cols = cols
        else:
            fluor_cols = [fluor_box.currentText() or (cols[0] if cols else None)]

        annot_box.blockSignals(True)
        annot_box.setText(store.note(well_site, particle))
        annot_box.blockSignals(False)
        state["note_owner"] = (well_site, particle)

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

        drawn = show_traces(well_site, particle, track, fluor_cols)
        show_info(well_site, particle, track, frames)
        move_cursor()
        exclude_btn.setText("Include this particle"
                            if store.is_excluded(well_site, particle)
                            else "Exclude this particle")
        missing = [c for c in ("dead_proba",) if c not in drawn]
        n_excl = len(store.excluded[well_site])
        n_total = len(store.particles(well_site))
        group = group_box.currentData()
        shown = (f"{part_box.count()} of {n_total} particles "
                 f"in {group_by_box.currentData()} = {group}"
                 if group is not None else f"{n_total} particles")
        status.setText(f"{well_site}: {shown}, {n_excl} excluded. Channels: "
                       f"{', '.join(sorted(pos.stacks)) or 'none'}."
                       + (f"<br>Not plotted (absent from this workbook): "
                          f"{', '.join(missing)}" if missing else ""))

    def refill_particles(keep=None):
        well_site = site_box.currentData()
        part_box.blockSignals(True)
        part_box.clear()
        if well_site:
            try:
                particles = store.particles_in_group(
                    well_site, group_by_box.currentData(), group_box.currentData())
            except Exception as exc:  # unreadable workbook
                particles = []
                status.setText(f"<b>could not read summary for {well_site}: {exc}</b>")
            for p in particles:
                part_box.addItem(part_label(well_site, p), p)
            # Stay on the particle already under review when the group list is
            # rebuilt around it, rather than jumping back to the top.
            if keep in particles:
                part_box.setCurrentIndex(particles.index(keep))
        part_box.blockSignals(False)
        load_particle()
        start_prefetch()

    def start_prefetch():
        """Warm the whole listed group in one pass over the stacks.

        Kicked off after the particle list is rebuilt, i.e. on every position,
        group-by or group change. The particle already on screen is loaded
        first by load_particle, and anything the pass has not reached yet still
        works through the ordinary on-demand read.
        """
        well_site = site_box.currentData()
        if not well_site:
            prefetch.cancel()
            return
        pos = store.positions[well_site]
        particles = [part_box.itemData(i) for i in range(part_box.count())]
        tracks = {p: store.track(well_site, p) for p in particles
                  if p is not None}
        prefetch.start(pos, tracks)

    def drain_prefetch():
        if prefetch.drain() and prefetch.running:
            cache_label.setText(f"{roi_cache.stats()} - caching group "
                                f"{prefetch.done}/{prefetch.total}")
        elif not prefetch.running:
            cache_label.setText(roi_cache.stats())

    def refill_groups(keep=None):
        """List this position's groups, with counts, for the chosen column."""
        well_site = site_box.currentData()
        group_by = group_by_box.currentData()
        wanted = group_box.currentData()
        group_box.blockSignals(True)
        group_box.clear()
        if well_site and group_by:
            try:
                groups = store.groups(well_site, group_by)
            except Exception:
                groups = []
            group_box.addItem(f"all ({sum(n for _, n in groups)})", None)
            for name, n in groups:
                group_box.addItem(f"{name}  ({n})", name)
            # Hold the same group across positions when it exists there too.
            hit = group_box.findData(wanted)
            group_box.setCurrentIndex(max(hit, 0))
        group_box.setEnabled(bool(group_by))
        group_box.blockSignals(False)
        refill_particles(keep)

    def refill_group_bys():
        """Offer the categorical summary columns, plus the user's own notes."""
        well_site = site_box.currentData()
        # An empty box means nothing has been chosen yet, which is not the same
        # as having chosen "(none)" - both read back as None.
        wanted = group_by_box.currentData() if group_by_box.count() else "fate_label"
        group_by_box.blockSignals(True)
        group_by_box.clear()
        group_by_box.addItem(NO_GROUP, None)
        if well_site:
            try:
                for col in store.group_columns(well_site):
                    group_by_box.addItem(col, col)
            except Exception:
                pass
            group_by_box.addItem(NOTE_GROUP, NOTE_GROUP)
        # Falls back to "(none)" at index 0 when the wanted column is not in
        # this position's summary, as in the pre-classifier workbooks.
        group_by_box.setCurrentIndex(max(group_by_box.findData(wanted), 0))
        group_by_box.blockSignals(False)
        refill_groups()

    def refill_sites():
        well = well_box.currentData()
        site_box.blockSignals(True)
        site_box.clear()
        for ws in store.sites(well) if well else []:
            pos = store.positions[ws]
            text = pos.site if pos.stacks else f"{pos.site}  (no images)"
            site_box.addItem(text, ws)
        site_box.blockSignals(False)
        refill_group_bys()

    def rescan():
        flush_note()  # rescanning reloads curation state from disk
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
            for box in (site_box, group_by_box, group_box, part_box):
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
        flush_note()
        # A position counts as curated if anything was excluded or annotated.
        touched = ({ws for ws, ps in store.excluded.items() if ps}
                   | {ws for ws, ns in store.notes.items() if any(ns.values())})
        targets = [well_site]
        others = sorted(touched - {well_site})
        if others:
            answer = QMessageBox.question(
                panel, "Export scope",
                f"{len(others)} other position(s) also have exclusions or "
                f"annotations ({', '.join(others)}).\n\n"
                "Yes: export every curated position.\n"
                f"No: export only {well_site}.",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if answer == QMessageBox.Yes:
                targets = sorted(touched | {well_site})

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
            n_notes = sum(1 for p in store.kept(ws) if store.note(ws, p))
            written.append(f"{ws}: {n_kept} particles kept ({n_drop} excluded), "
                           f"{n_rows} tracking rows, {n_notes} annotated")
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
    site_box.currentIndexChanged.connect(lambda _: refill_group_bys())
    group_by_box.currentIndexChanged.connect(lambda _: refill_groups())
    group_box.currentIndexChanged.connect(lambda _: refill_particles())
    part_box.currentIndexChanged.connect(lambda _: load_particle())
    fluor_box.currentIndexChanged.connect(lambda _: load_particle())
    fluor_all.toggled.connect(lambda _: load_particle())
    annot_box.textEdited.connect(on_note_edited)
    annot_box.editingFinished.connect(commit_note)
    prev_btn.clicked.connect(lambda: step_particle(-1))
    next_btn.clicked.connect(lambda: step_particle(1))
    exclude_btn.clicked.connect(on_exclude)
    export_btn.clicked.connect(on_export)
    viewer.dims.events.current_step.connect(move_cursor)
    canvas.mpl_connect("button_press_event", on_plot_click)

    # Hands finished particles from the background pass to the cache. Polling
    # rather than a cross-thread signal keeps every cache mutation on the GUI
    # thread, so RoiCache needs no locking.
    drain_timer = QTimer(panel)
    drain_timer.timeout.connect(drain_prefetch)
    drain_timer.start(250)

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
