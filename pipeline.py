#!/usr/bin/env python
"""One driver for the whole cellaap workflow: platemap -> inference -> analysis.

The two halves of the workflow need different conda environments (inference
needs cell_AAP/detectron2 and a GPU, analysis needs img-env and only CPUs), so
they cannot run in one process. What they CAN share is a description of the
plate, and that is what this script is built around: a single `platemap.csv`
the user fills in once, listing what is in each well. Everything else - which
model to run, which cell-type parameters to analyze with, which fluorescence
channels to measure - is derived from it, so the two stages can never disagree
about what a well contains.

Typical use, from the module directory:

    python pipeline.py init    --root /scratch/.../20251025      # write template
    $EDITOR /scratch/.../20251025/platemap.csv                   # fill it in
    python pipeline.py check   --root /scratch/.../20251025      # dry run
    python pipeline.py submit  --root /scratch/.../20251025 --sbatch
    python pipeline.py status  --root /scratch/.../20251025

`submit` writes two SLURM array jobs - one GPU array for inference, one CPU
array for analysis - with the analysis array held on `aftercorr` against the
inference array, so each position is analyzed as soon as its own inference
lands rather than waiting for the whole plate. When the platemap declares
blank-media wells, a third small job builds the fluorescence correction maps
from them first, and the analysis array waits on that too. Those wells hold no
cells, so they are never segmented; declare none and the analysis simply runs
uncorrected.

Reliability, i.e. the things the old batch_* scripts got bitten by:

* One position per array task. A position that fails takes nothing else down.
* Inference writes to `<dir>.partial` and renames on success, so a job killed
  by a walltime limit cannot leave a half-written inference folder that later
  runs mistake for a finished one.
* Every finished position drops a JSON marker under `pipeline/state/`
  recording the parameters it ran with. Re-running skips finished work, and
  changing the platemap marks the affected positions stale so they re-run.
* The platemap is frozen (copied) into the job directory at submit time, so
  editing it while jobs are queued does not change what those jobs do.
* Everything is validated - wells mapped, phase files readable, channels
  present, model names real - before a single job is submitted.

This file deliberately imports nothing heavy at module scope: it has to run
under both conda environments, so `inference` and `cellaap_analysis` are
imported inside the subcommands that need them.
"""

from __future__ import annotations

import argparse
import csv
import getpass
import json
import os
import re
import resource
import shlex
import shutil
import socket
import subprocess
import sys
import time
import traceback
from dataclasses import dataclass, field as dc_field
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULT_CELLTYPE = "HeLa"

# celltype -> (inference model, analysis_pars cell_type). The focal models are
# the current generation and are what the recent runs used; override per well
# with the `model` column if a plain model is wanted.
CELLTYPES = {
    "hela": ("HeLa_focal", "hela"),
    "u2os": ("U2OS_focal", "u2os"),
    "rpe1": ("RPE1_focal", "rpe1"),
    "ht1080": ("HT1080_focal", "ht1080"),
}

# Kept in sync with inference.get_model(); duplicated here so `check` can
# validate a platemap without importing detectron2.
VALID_MODELS = [
    "HeLa", "HeLa_focal", "HT1080_focal", "HT1080", "RPE1_focal", "RPE1",
    "U2OS_focal", "U2OS", "general",
]

# Channel names as they appear in file names, per cellaap_analysis.files().
VALID_CHANNELS = ["GFP", "Texas Red", "Cy5"]

# What a well holds. Blank means cells. The two blank-media roles name the
# medium rather than the correction they produce, because that is what the
# person at the microscope knows: fluorobrite is the background well, DMEM the
# excitation-intensity well.
ROLES = {
    "": "sample", "sample": "sample", "cells": "sample",
    "background": "background", "fluorobrite": "background", "bkg": "background",
    "intensity": "intensity", "dmem": "intensity",
}
MAP_ROLES = ("background", "intensity")

# DMEM autofluoresces; FluoroBrite is formulated not to. So the DMEM well is
# always the brighter of the two, and a pair of blank wells whose brightness
# runs the other way has been mislabeled. Below this much relative difference
# the two wells are too close to call, and the labels are left alone rather
# than swapped on noise.
ROLE_CHECK_MARGIN = 0.05
ROLE_CHECK_FRAMES = 3

DEFAULT_CONFLUENCY = 1800
DEFAULT_CONF_THRESHOLD = 0.25

# The position stub cellaap_analysis keys everything off: well + site.
STUB_RE = re.compile(r"(?P<well>[A-H](?:0[1-9]|1[0-2]|[1-9]))_s(?P<site>\d{1,2})")
WELL_RE = re.compile(r"^(?P<row>[A-Ha-h])(?P<col>\d{1,2})$")

PIPELINE_DIRNAME = "pipeline"
PLATEMAP_NAME = "platemap.csv"

PLATEMAP_COLUMNS = [
    "well_ids", "role", "celltype", "transfection", "drug", "channels",
    "model", "confluency_est", "conf_threshold", "skip", "notes",
]

PLATEMAP_HEADER = """\
# cellaap platemap - fill this in before running the pipeline.
#
# One row per experimental condition. Lines starting with # are ignored.
#
# well_ids       wells this row describes: A01, or "A01,A02 B01", or a range
#                A01-A04. A single position can be named directly as A01_s3
#                when one site in a well needs different treatment.
# role           blank for wells with cells in them. For the blank-media wells
#                used to correct the fluorescence, put:
#                  background  (or fluorobrite) - the fluorobrite well
#                  intensity   (or dmem)        - the DMEM well
#                These wells are NOT segmented; they are used only to build the
#                correction maps, one per channel, before the analysis runs.
#                Leave them out entirely and the analysis proceeds with no
#                background or intensity correction.
# celltype       one of: {celltypes}. Sets BOTH the inference model and the
#                analysis parameters. Unlisted wells default to {default}.
#                Ignored on background/intensity rows.
# transfection   free text, e.g. pEN2 or siBUB1. Use none/-- if not applicable.
# drug           free text, e.g. DMSO or "STLC 5uM".
# channels       fluorescence channels to measure, e.g. "Texas Red" or
#                "GFP;Texas Red". Leave BLANK to measure every channel found
#                next to the phase file. Valid: {channels}. On a
#                background/intensity row this is which channels that blank
#                well provides a map for.
# model          optional. Overrides the model implied by celltype.
#                One of: {models}.
# confluency_est optional, (0, 2000]. Default {confluency}.
# conf_threshold optional, (0, 1). Default {threshold}.
# skip           put y/yes/true here to leave the well out of the run.
# notes          free text, ignored by the pipeline.
#
# celltype/transfection/drug are also the grouping keys used by
# cellaap_aggregate.load_experiment(), so this same file drives the
# downstream compilation of summaries.
""".format(
    celltypes=", ".join(sorted(CELLTYPES)),
    default=DEFAULT_CELLTYPE,
    channels=", ".join(VALID_CHANNELS),
    models=", ".join(VALID_MODELS),
    confluency=DEFAULT_CONFLUENCY,
    threshold=DEFAULT_CONF_THRESHOLD,
)

def default_mail_user() -> str:
    """<login>@umich.edu for whoever is submitting.

    Resolved here rather than written as "$USER@umich.edu" into the script:
    SBATCH directives are not shell-expanded, so the literal string would
    reach SLURM as an invalid address.
    """
    try:
        user = os.environ.get("USER") or getpass.getuser()
    except Exception:
        return ""
    return f"{user}@umich.edu" if user else ""


# SLURM defaults. The walltimes and memory figures are per array task - one
# position for the two stages, one folder for the maps job - and are sized to
# measured runs rather than to the worst case. A task that exceeds its walltime
# is killed by SLURM, so raise the relevant flag and resubmit; finished
# positions are skipped, so a resubmit only re-runs what died.
SLURM_DEFAULTS = {
    "account": "ajitj99",
    "mail_user": default_mail_user(),
    "gpu_partition": "gpu",
    "cpu_partition": "standard",
    "infer_env": "cellaap-env",
    "analysis_env": "img-env",
    "infer_time": "0-00:40:00",
    "analysis_time": "0-01:00:00",
    "infer_mem": "12g",
    "analysis_mem": "25g",
    "analysis_cpus": 1,
    "maps_time": "0-00:30:00",
    "maps_mem": "20g",
    "infer_concurrent": 4,
    "analysis_concurrent": 12,
}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str) -> "NoReturn":  # type: ignore[valid-type]
    print(f"ERROR: {msg}", file=sys.stderr, flush=True)
    raise SystemExit(2)


def peak_rss_bytes() -> int:
    """Peak resident set size of this process, in bytes.

    Under SLURM one array task is one position in its own process, so this is
    that position's peak. Run locally over a whole plate it is the high-water
    mark of the process so far, which over-reports every position after the
    heaviest one - `report` says so rather than pretending otherwise.

    ru_maxrss is kilobytes on Linux and bytes on macOS.
    """
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(raw) if sys.platform == "darwin" else int(raw) * 1024


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def normalize_well(token: str) -> str:
    """A01 / a1 / A1 -> A01. Returns '' if the token is not a well id."""
    m = WELL_RE.match(token.strip())
    if not m:
        return ""
    return f"{m.group('row').upper()}{int(m.group('col')):02d}"


def expand_well_token(token: str) -> list[str]:
    """Expand one well_ids token: a well, a position (A01_s3), or a range."""
    token = token.strip()
    if not token:
        return []

    # An explicit position, e.g. A01_s3 - kept verbatim as a position key.
    m = STUB_RE.fullmatch(token.replace(" ", ""))
    if m:
        return [f"{normalize_well(m.group('well'))}_s{int(m.group('site'))}"]

    if "-" in token:
        lo, _, hi = token.partition("-")
        lo, hi = normalize_well(lo), normalize_well(hi)
        if lo and hi and lo[0] == hi[0]:
            return [f"{lo[0]}{c:02d}" for c in range(int(lo[1:]), int(hi[1:]) + 1)]

    well = normalize_well(token)
    return [well] if well else [token.strip()]


def split_list(value: str) -> list[str]:
    """Split a user-typed list on comma/semicolon. Channel names have spaces
    in them ('Texas Red'), so whitespace cannot be a separator here."""
    return [p.strip() for p in re.split(r"[,;]+", value or "") if p.strip()]


def as_bool(value: str) -> bool:
    return str(value).strip().lower() in {"y", "yes", "true", "1", "skip"}


# ---------------------------------------------------------------------------
# Positions and platemap
# ---------------------------------------------------------------------------

@dataclass
class Position:
    """One microscope position: a `*phs.tif` stack and its parsed stub."""
    stem: str          # file name without .tif, e.g. 20251009_HeLa_G03_s8_phs
    phs: Path
    well: str          # G03
    site: int          # 8

    @property
    def stub(self) -> str:
        return f"{self.well}_s{self.site}"


@dataclass
class PlateRow:
    wells: list[str]
    role: str
    celltype: str
    transfection: str
    drug: str
    channels: list[str]
    model: str
    confluency_est: int | None
    conf_threshold: float | None
    skip: bool
    notes: str
    lineno: int


@dataclass
class MapWell:
    """A blank-media position: fluorobrite (background) or DMEM (intensity)."""
    pos: Position
    role: str
    channels: list[str]
    lineno: int


@dataclass
class MapPair:
    """One correction map to build: this channel, from this stack.

    The output name is not ours to choose - create_correction_maps derives it
    from the source file stem, and cellaap_analysis._load_maps finds it again
    by looking for the channel name and the words background/intensity in the
    file name.
    """
    channel: str
    role: str
    source: Path

    @property
    def output(self) -> Path:
        return self.source.parent / f"{self.source.stem}_{self.role}_map.tif"

    @property
    def key(self) -> str:
        return f"{self.channel}_{self.role}"


@dataclass
class Task:
    """A position plus everything both stages need to know about it."""
    pos: Position
    celltype: str
    transfection: str
    drug: str
    model: str
    confluency_est: int
    conf_threshold: float
    analysis_cell_type: str
    channels: list[str]
    mapped: bool           # False when the well fell through to the default
    root: Path
    # signal_correction knobs, carried on the task so that changing one marks
    # finished positions stale rather than silently leaving them uncorrected
    # under the old settings.
    correction_params: dict = dc_field(default_factory=dict)

    @property
    def stem(self) -> str:
        return self.pos.stem

    @property
    def inference_dir(self) -> Path:
        """Must match what batch_inference_APJ.py produced, so that folders
        made by the old scripts are picked up unchanged."""
        name = (f"{self.pos.stem}_{self.model}_{self.confluency_est}_"
                f"{round(self.conf_threshold, 2)}_inference")
        return self.root / name

    def params(self, stage: str) -> dict:
        """The parameters whose change should invalidate a finished stage."""
        if stage == "inference":
            return {"model": self.model, "confluency_est": self.confluency_est,
                    "conf_threshold": round(self.conf_threshold, 2)}
        return {"cell_type": self.analysis_cell_type,
                "channels": sorted(self.channels),
                "inference_dir": self.inference_dir.name,
                # Changing how the signal is corrected has to re-run the
                # analysis, the same as changing the cell type would.
                "correction": dict(self.correction_params)}


def discover_positions(root: Path, pattern: str = "*phs.tif") -> list[Position]:
    positions = []
    for phs in sorted(root.glob(pattern)):
        m = STUB_RE.search(phs.name)
        if not m:
            log(f"  ! {phs.name}: no <well>_s<site> in the name, skipped")
            continue
        positions.append(Position(stem=phs.name[:-len(phs.suffix)], phs=phs,
                                  well=normalize_well(m.group("well")),
                                  site=int(m.group("site"))))
    return positions


def platemap_path(root: Path) -> Path:
    return root / PLATEMAP_NAME


def write_template(root: Path, positions: list[Position], force: bool) -> Path:
    path = platemap_path(root)
    if path.exists() and not force:
        die(f"{path} already exists. Edit it, or pass --force to overwrite.")

    wells: dict[str, list[Position]] = {}
    for p in positions:
        wells.setdefault(p.well, []).append(p)

    with open(path, "w", newline="") as fh:
        fh.write(PLATEMAP_HEADER)
        writer = csv.DictWriter(fh, fieldnames=PLATEMAP_COLUMNS)
        writer.writeheader()
        for well in sorted(wells):
            sites = ",".join(f"s{p.site}" for p in sorted(wells[well],
                                                          key=lambda x: x.site))
            writer.writerow({
                "well_ids": well,
                "role": "",
                "celltype": DEFAULT_CELLTYPE,
                "transfection": "",
                "drug": "",
                "channels": "",
                "model": "",
                "confluency_est": "",
                "conf_threshold": "",
                "skip": "",
                "notes": f"{len(wells[well])} position(s): {sites}",
            })
    return path


def read_platemap(path: Path) -> list[PlateRow]:
    if not path.exists():
        die(f"No platemap at {path}. Run `pipeline.py init --root {path.parent}` first.")

    rows: list[PlateRow] = []
    with open(path, newline="") as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
    reader = csv.DictReader(lines)
    if reader.fieldnames is None:
        die(f"{path} has no header row.")
    missing = {"well_ids", "celltype"} - set(reader.fieldnames)
    if missing:
        die(f"{path} is missing required column(s): {', '.join(sorted(missing))}")

    for lineno, raw in enumerate(reader, start=2):
        raw = {k: (v or "").strip() for k, v in raw.items() if k}
        if not raw.get("well_ids"):
            continue

        wells: list[str] = []
        for token in re.split(r"[,;\s]+", raw["well_ids"]):
            wells.extend(expand_well_token(token))

        def number(key, cast, default=None):
            value = raw.get(key, "")
            if not value:
                return default
            try:
                return cast(value)
            except ValueError:
                die(f"{path} line {lineno}: {key}={value!r} is not a number")

        role = ROLES.get(raw.get("role", "").lower())
        if role is None:
            die(f"{path} line {lineno}: role={raw.get('role')!r} is not one of "
                f"{', '.join(sorted(set(ROLES) - {''}))}")

        rows.append(PlateRow(
            wells=wells,
            role=role,
            celltype=raw.get("celltype") or DEFAULT_CELLTYPE,
            transfection=raw.get("transfection", ""),
            drug=raw.get("drug", ""),
            channels=split_list(raw.get("channels", "")),
            model=raw.get("model", ""),
            confluency_est=number("confluency_est", int),
            conf_threshold=number("conf_threshold", float),
            skip=as_bool(raw.get("skip", "")),
            notes=raw.get("notes", ""),
            lineno=lineno,
        ))
    return rows


def detect_channels(root: Path, pos: Position) -> list[str]:
    """Channels with a stack sitting next to the phase file for this position.

    cellaap_analysis.files() finds channel stacks by the same rule - a file
    carrying this position's stub and a known channel name - so anything this
    returns is measurable.
    """
    found = []
    for path in sorted(root.glob("*.tif")) + sorted(root.glob("*.tiff")):
        if f"{pos.stub}_" not in path.name:
            continue
        m = re.search(r"GFP|Texas Red|Cy5", path.name)
        if m and m.group() not in found:
            found.append(m.group())
    return found


def channel_stack(root: Path, pos: Position, channel: str) -> Path | None:
    """The stack for one channel of one position, by the same naming rule
    cellaap_analysis.files() uses."""
    for path in sorted(root.glob("*.tif")) + sorted(root.glob("*.tiff")):
        if f"{pos.stub}_" in path.name and channel in path.name:
            return path
    return None


def resolve_correction_maps(root: Path, map_wells: list[MapWell]
                            ) -> tuple[list[MapPair], list[str]]:
    """Which correction maps the platemap asks for, and from which stack.

    One map per channel per role - a second fluorobrite well for the same
    channel would produce a second file that cellaap_analysis._load_maps picks
    up in whatever order os.walk returns, so extras are refused here rather
    than left to chance.
    """
    warnings: list[str] = []
    pairs: dict[str, MapPair] = {}
    for well in sorted(map_wells, key=lambda w: (w.pos.well, w.pos.site)):
        for channel in well.channels:
            source = channel_stack(root, well.pos, channel)
            if source is None:
                warnings.append(f"{well.pos.stub} is marked {well.role} for "
                                f"{channel}, but no {channel} stack was found "
                                f"for it")
                continue
            key = f"{channel}_{well.role}"
            if key in pairs:
                warnings.append(f"more than one {well.role} well provides "
                                f"{channel}; using {pairs[key].source.name} and "
                                f"ignoring {source.name}")
                continue
            pairs[key] = MapPair(channel=channel, role=well.role, source=source)
    return list(pairs.values()), warnings


@dataclass
class RoleCheck:
    """Whether the two blank wells of one channel are labeled the right way."""
    channel: str
    intensity_mean: float
    background_mean: float
    verdict: str                 # ok | swapped | inconclusive | not-checked
    intensity_source: str        # after any swap
    background_source: str

    @property
    def ratio(self) -> float:
        return (self.intensity_mean / self.background_mean
                if self.background_mean else float("nan"))

    def as_dict(self) -> dict:
        return {"channel": self.channel, "verdict": self.verdict,
                "dmem_mean": round(self.intensity_mean, 1),
                "fluorobrite_mean": round(self.background_mean, 1),
                "ratio": round(self.ratio, 3),
                "intensity_from": self.intensity_source,
                "background_from": self.background_source}


def sample_mean(path: Path, frames: int = ROLE_CHECK_FRAMES) -> float | None:
    """Mean intensity of a few frames of a stack.

    A few pages, not the whole stack: these are gigabyte files on shared
    storage, and this only has to separate DMEM from FluoroBrite, which differ
    severalfold. Returns None if the file cannot be read.
    """
    try:
        import numpy as np
        import tifffile
        with tifffile.TiffFile(path) as fh:
            n = len(fh.pages)
            if n == 0:
                return None
            idx = sorted({int(round(i * (n - 1) / max(frames - 1, 1)))
                          for i in range(min(frames, n))})
            values = [float(np.asarray(fh.pages[i].asarray()).mean()) for i in idx]
        return sum(values) / len(values)
    except Exception:
        return None


def verify_map_roles(pairs: list[MapPair]) -> tuple[list[MapPair], list[RoleCheck]]:
    """Check each channel's blank wells against the brightness rule, and swap
    the two sources when they have been labeled the wrong way round.

    Only wells the platemap marks `dmem` build intensity maps and only wells it
    marks `fluorobrite` build background maps - this makes sure the label
    matches what is actually in the well. The pair is swapped, not dropped,
    because a mislabeled pair still contains both media; it is the names that
    are wrong.
    """
    checks: list[RoleCheck] = []
    by_channel: dict[str, dict[str, MapPair]] = {}
    for pair in pairs:
        by_channel.setdefault(pair.channel, {})[pair.role] = pair

    for channel, roles in sorted(by_channel.items()):
        intensity, background = roles.get("intensity"), roles.get("background")
        if not (intensity and background):
            continue

        i_mean = sample_mean(intensity.source)
        b_mean = sample_mean(background.source)
        if i_mean is None or b_mean is None or b_mean <= 0:
            checks.append(RoleCheck(channel, i_mean or float("nan"),
                                    b_mean or float("nan"), "not-checked",
                                    intensity.source.name, background.source.name))
            continue

        if i_mean >= b_mean * (1 + ROLE_CHECK_MARGIN):
            verdict = "ok"
        elif b_mean >= i_mean * (1 + ROLE_CHECK_MARGIN):
            verdict = "swapped"
            intensity.source, background.source = background.source, intensity.source
            i_mean, b_mean = b_mean, i_mean
        else:
            verdict = "inconclusive"

        checks.append(RoleCheck(channel, i_mean, b_mean, verdict,
                                intensity.source.name, background.source.name))
    return pairs, checks


def report_role_checks(checks: list[RoleCheck]) -> list[str]:
    """The lines to show the user. Loud for a swap - it means the platemap is
    wrong, and the platemap is what everything else is read from."""
    lines = []
    for c in checks:
        if c.verdict == "ok":
            lines.append(f"  role check {c.channel}: dmem {c.intensity_mean:.0f} vs "
                         f"fluorobrite {c.background_mean:.0f} ({c.ratio:.1f}x) - "
                         f"labels look right")
        elif c.verdict == "swapped":
            lines.append(
                f"  ! role check {c.channel}: the well labeled dmem was DIMMER "
                f"than the well labeled fluorobrite. DMEM autofluoresces and "
                f"FluoroBrite does not, so these labels are swapped in the "
                f"platemap. SWAPPING THEM: intensity map from "
                f"{c.intensity_source}, background map from {c.background_source} "
                f"(means now {c.intensity_mean:.0f} vs {c.background_mean:.0f}). "
                f"Fix the platemap to make this permanent.")
        elif c.verdict == "inconclusive":
            lines.append(
                f"  ! role check {c.channel}: the two blank wells are within "
                f"{ROLE_CHECK_MARGIN:.0%} of each other "
                f"({c.intensity_mean:.0f} vs {c.background_mean:.0f}), too close "
                f"to tell apart. Labels left as written - check they are right.")
        else:
            lines.append(f"  ! role check {c.channel}: could not read one of the "
                         f"blank stacks; labels left as written")
    return lines


def find_map_files(root: Path, channel: str, role: str) -> list[Path]:
    """Correction-map files already in the folder for this channel and role.

    Matches the way _load_maps finds them - the channel name and the word
    background/intensity somewhere in the file name - so that this sees
    exactly what the analysis will see, including maps made by hand.
    """
    out = []
    for path in root.rglob("*.tif"):
        if pipeline_dir(root) in path.parents:
            continue
        if channel in path.name and role in path.name and "_map" in path.name:
            out.append(path)
    return sorted(out)


def check_map_files(root: Path, pairs: list[MapPair]) -> list[str]:
    """Two ways the correction maps already on disk can break every analysis.

    Both are properties of cellaap_analysis._load_maps(), which runs at the
    start of every analysis: it walks the whole root folder, treats any file
    whose name contains background/intensity as a map, and reads the channel
    out of the same name.
    """
    problems = []

    # Two maps for one channel and role: _load_maps keeps whichever os.walk
    # yields last, so which one is applied is not defined.
    for pair in pairs:
        existing = find_map_files(root, pair.channel, pair.role)
        if len(existing) > 1:
            problems.append(
                f"{len(existing)} {pair.role} maps for {pair.channel} in the "
                f"folder ({', '.join(p.name for p in existing)}); the analysis "
                f"would pick one at random - delete all but the right one")

    # A file that looks like a map but names no channel: _load_maps calls
    # .group() on a failed match and dies before any analysis starts.
    for path in root.rglob("*"):
        if not path.is_file() or pipeline_dir(root) in path.parents:
            continue
        if re.search(r"background|intensity", path.name) and \
                not re.search(r"GFP|Texas Red|Cy5|phs", path.name):
            problems.append(
                f"{path.name} has 'background' or 'intensity' in its name but "
                f"no channel name; cellaap_analysis reads it as a correction "
                f"map and every analysis in this folder will crash on it")

    return problems


def quarantine_stale_maps(root: Path, pair: MapPair) -> list[Path]:
    """Move aside correction maps for this channel and role that came from a
    different stack than the one we are about to use.

    After a swap the folder would otherwise hold two intensity maps for one
    channel - the wrong one from the earlier run and the right one from this
    one - and _load_maps picks between them by directory order. Moved, not
    deleted: they are somebody's data, and the move is recorded in the log.
    """
    moved = []
    for path in find_map_files(root, pair.channel, pair.role):
        if path.name == pair.output.name:
            continue
        destination = pipeline_dir(root) / "superseded" / path.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(path), str(destination))
        moved.append(destination)
    return moved


def correction_params(args) -> dict:
    """The signal_correction knobs, from whichever subcommand is running.

    Every subcommand carries the same three flags with the same defaults, so
    this returns the same dict everywhere and a finished position only reports
    stale when the operator actually changed one.
    """
    return {"dilation": getattr(args, "correction_dilation", 121),
            "n_frames": getattr(args, "correction_frames", 24),
            "block": getattr(args, "correction_block", 64),
            "model": getattr(args, "correction_model", "auto")}


def build_tasks(root: Path, positions: list[Position], rows: list[PlateRow],
                autodetect_channels: bool = True, correction: dict | None = None
                ) -> tuple[list[Task], list[MapWell], list[str]]:
    """Join positions to platemap rows. Returns (tasks, map_wells, warnings).

    Precedence is exact position (G03_s8) over well (G03) over the HeLa
    default, so a single odd site can be overridden without splitting the well.

    Blank-media wells come back separately: they hold no cells, so segmenting
    them would be meaningless work, and neither stage should ever see them as
    a position to process.
    """
    warnings: list[str] = []
    by_position: dict[str, PlateRow] = {}
    by_well: dict[str, PlateRow] = {}
    for row in rows:
        for well in row.wells:
            target = by_position if "_s" in well else by_well
            if well in target:
                warnings.append(f"{well} appears twice in the platemap "
                                f"(lines {target[well].lineno} and {row.lineno}); "
                                f"line {row.lineno} wins")
            target[well] = row

    seen_wells = set()
    tasks: list[Task] = []
    map_wells: list[MapWell] = []
    for pos in positions:
        seen_wells.add(pos.well)
        row = by_position.get(pos.stub) or by_well.get(pos.well)
        mapped = row is not None
        if row is None:
            row = PlateRow(wells=[pos.well], role="sample",
                           celltype=DEFAULT_CELLTYPE,
                           transfection="", drug="", channels=[], model="",
                           confluency_est=None, conf_threshold=None,
                           skip=False, notes="", lineno=-1)
            warnings.append(f"well {pos.well} is not in the platemap; "
                            f"falling back to {DEFAULT_CELLTYPE}")
        if row.skip:
            continue

        if row.role in MAP_ROLES:
            channels = list(row.channels)
            bad = [c for c in channels if c not in VALID_CHANNELS]
            if bad:
                die(f"platemap line {row.lineno}: unknown channel(s) "
                    f"{', '.join(bad)}; valid names are {', '.join(VALID_CHANNELS)}")
            if not channels and autodetect_channels:
                channels = detect_channels(root, pos)
            if not channels:
                warnings.append(f"{pos.stub} is marked {row.role} but has no "
                                f"fluorescence stack; no map can be built from it")
            map_wells.append(MapWell(pos=pos, role=row.role, channels=channels,
                                     lineno=row.lineno))
            continue

        celltype = row.celltype.strip()
        key = celltype.lower()
        if key not in CELLTYPES:
            die(f"platemap line {row.lineno}: celltype {celltype!r} is not one of "
                f"{', '.join(sorted(CELLTYPES))}")
        default_model, analysis_cell_type = CELLTYPES[key]

        model = row.model or default_model
        if model not in VALID_MODELS:
            die(f"platemap line {row.lineno}: model {model!r} is not one of "
                f"{', '.join(VALID_MODELS)}")

        channels = list(row.channels)
        if channels:
            bad = [c for c in channels if c not in VALID_CHANNELS]
            if bad:
                die(f"platemap line {row.lineno}: unknown channel(s) "
                    f"{', '.join(bad)}; valid names are {', '.join(VALID_CHANNELS)}")
        elif autodetect_channels:
            channels = detect_channels(root, pos)

        confluency = row.confluency_est or DEFAULT_CONFLUENCY
        if not 0 < confluency <= 2000:
            die(f"platemap line {row.lineno}: confluency_est must be in (0, 2000]")
        threshold = row.conf_threshold or DEFAULT_CONF_THRESHOLD
        if not 0 < threshold < 1:
            die(f"platemap line {row.lineno}: conf_threshold must be in (0, 1)")

        tasks.append(Task(
            pos=pos, celltype=celltype, transfection=row.transfection,
            drug=row.drug, model=model, confluency_est=confluency,
            conf_threshold=threshold, analysis_cell_type=analysis_cell_type,
            channels=channels, mapped=mapped, root=root,
            correction_params=dict(correction or {}),
        ))

    # A row covering a whole plate row (B01-B12) when only three wells were
    # imaged is normal; report the unused wells once rather than per well.
    unused = [w for row in rows for w in row.wells
              if w.split("_s")[0] not in seen_wells]
    if unused:
        warnings.append(f"{len(unused)} well(s) in the platemap have no phase "
                        f"stack: {', '.join(sorted(unused))}")

    return tasks, map_wells, warnings


# ---------------------------------------------------------------------------
# State: what has finished, and with which parameters
# ---------------------------------------------------------------------------

def pipeline_dir(root: Path) -> Path:
    return root / PIPELINE_DIRNAME


def state_file(root: Path, stage: str, stem: str) -> Path:
    return pipeline_dir(root) / "state" / stage / f"{stem}.json"


def read_state(root: Path, stage: str, stem: str) -> dict | None:
    path = state_file(root, stage, stem)
    if not path.exists():
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError):
        return None


def write_state(root: Path, stage: str, task: Task, status: str,
                started: float, error: str = "", outputs: list[str] | None = None) -> None:
    path = state_file(root, stage, task.stem)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": stage,
        "stem": task.stem,
        "well": task.pos.well,
        "site": task.pos.site,
        "status": status,
        "started": datetime.fromtimestamp(started).isoformat(timespec="seconds"),
        "finished": now(),
        "duration_s": round(time.time() - started, 1),
        "params": task.params(stage),
        "host": socket.gethostname(),
        "slurm_job": os.environ.get("SLURM_JOB_ID", ""),
        "peak_rss_bytes": peak_rss_bytes(),
        # True when this process handled one position only, which is what makes
        # peak_rss_bytes attributable to it. False for a local sequential run.
        "rss_is_exclusive": bool(os.environ.get("SLURM_ARRAY_TASK_ID")),
        "outputs": outputs or [],
        "error": error,
    }
    tmp = path.with_suffix(".json.tmp")
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)
    os.replace(tmp, path)


def inference_outputs_present(task: Task) -> bool:
    d = task.inference_dir
    if not d.is_dir():
        return False
    have = {kind: [p for p in d.glob("*.tif") if kind in p.name and p.stat().st_size > 0]
            for kind in ("semantic", "instance")}
    return all(have.values())


def analysis_outputs_present(task: Task) -> bool:
    d = task.inference_dir
    return d.is_dir() and any(p.stat().st_size > 0 for p in d.glob("*_summary.xlsx"))


def outputs_mtime(task: Task, stage: str) -> float | None:
    """When this stage's outputs were last written, or None if there are none."""
    d = task.inference_dir
    if not d.is_dir():
        return None
    if stage == "inference":
        files = [p for p in d.glob("*.tif")
                 if ("semantic" in p.name or "instance" in p.name)
                 and p.stat().st_size > 0]
    else:
        files = [p for p in d.glob("*_summary.xlsx") if p.stat().st_size > 0]
    return max((p.stat().st_mtime for p in files), default=None)


# Populated once per root per process by killed_index(). stage_status consults
# it for every position that looks adoptable, and re-walking the job
# directories each time would turn `check` into a filesystem crawl.
_KILLED_INDEX: dict = {}


def killed_index(root: Path) -> dict:
    """{(stage, stem): when SLURM last killed it} for this folder.

    Only kills that destroy work in progress count - a walltime or memory kill,
    a lost node - not a `scancel`, which is usually the user replacing one
    submit with another and says nothing about the outputs.
    """
    key = str(root)
    if key in _KILLED_INDEX:
        return _KILLED_INDEX[key]

    index: dict = {}
    for jobdir in job_dirs(root):
        for r in scan_job_dir(root, jobdir, with_logs=False):
            if (r["kind"] not in ("timeout", "oom", "node_fail", "preempted")
                    or not r["position"]):
                continue
            when = None
            m = CANCEL_AT_RE.search(r["evidence"])
            if m:
                try:
                    when = datetime.fromisoformat(m.group(1))
                except ValueError:
                    when = None
            if when is None:
                try:
                    when = datetime.fromtimestamp(
                        (pipeline_dir(root) / r["out"]).stat().st_mtime)
                except OSError:
                    continue
            k = (r["stage"], r["position"])
            if k not in index or when > index[k]:
                index[k] = when
    _KILLED_INDEX[key] = index
    return index


def stage_status(root: Path, task: Task, stage: str) -> str:
    """One of: done, stale, failed, pending.

    `stale` means the stage finished, but under parameters the platemap no
    longer asks for - it will be re-run. A finished folder with no marker (an
    old batch_* run) is adopted as done rather than redone, unless SLURM killed
    the position after those outputs were written - see below.
    """
    present = (inference_outputs_present(task) if stage == "inference"
               else analysis_outputs_present(task))
    state = read_state(root, stage, task.stem)

    if state is None:
        if not present:
            return "pending"
        # Adoption is for folders produced before markers existed. It becomes a
        # trap when a position was killed mid-run: the kill writes no marker, so
        # an older run's outputs are still sitting there and the position looks
        # finished. That is backwards - the killed work is exactly what needs
        # re-running - so outputs older than the kill do not count as done.
        killed_at = killed_index(root).get((stage, task.stem))
        written = outputs_mtime(task, stage)
        if killed_at and written and written < killed_at.timestamp():
            return "pending"
        return "done"
    if state.get("status") == "done":
        if not present:
            return "pending"      # marker without outputs: something removed them
        if state.get("params") != task.params(stage):
            return "stale"
        return "done"
    return "failed"


# The correction-map ledger is one file for all channels, and its name must
# contain neither "background" nor "intensity": cellaap_analysis._load_maps
# walks the whole root folder looking for those words in a file name and tries
# to imread whatever it finds, so a state file named after a role would be
# read as an image and crash every analysis in the plate.
# ---------------------------------------------------------------------------
# Flat fields for signal_correction
#
# The correction this pipeline applies is signal_correction's: the background
# is measured per position from that position's own frames, and the flat field
# - which is a property of the optics, not of the well - is built once for the
# plate from the *difference* of the two blank wells. The difference is what
# cancels the camera offset, which neither blank alone can be separated from;
# see SIGNAL_CORRECTION_DESIGN.md section 4.
#
# The old `*_map.tif` correction maps are no longer built. `pipeline.py maps`
# still exists for reproducing an old run by hand, but nothing submits it.
# ---------------------------------------------------------------------------

def surface_dir(root: Path) -> Path:
    """Where the per-position background surfaces are kept for inspection."""
    return pipeline_dir(root) / "state" / "surfaces"


def surface_path(root: Path, stem: str, channel: str) -> Path:
    """`<stem>_<channel>_bkg.tif` - deliberately NOT `..._background.tif`.

    cellaap_analysis._load_maps walks the whole root, this directory included,
    and treats any file whose name matches `background|intensity` as a
    plate-wide correction map. A per-position surface caught that way would be
    applied to every position on the plate. signal_correction.save_background_stack
    refuses such a name too; this is the belt to its braces.
    """
    safe = re.sub(r"[^A-Za-z0-9]+", "_", channel).strip("_")
    return surface_dir(root) / f"{stem}_{safe}_bkg.tif"


def flatfield_dir(root: Path) -> Path:
    return pipeline_dir(root) / "state" / "flatfield"


def flatfield_path(root: Path, channel: str) -> Path:
    """Where this channel's flat field is cached.

    The file name must contain neither "background" nor "intensity":
    cellaap_analysis._load_maps walks the whole root looking for those words
    and tries to imread whatever it finds. Channel names are safe on that
    count; the space in "Texas Red" is not safe as a file name, so it goes.
    """
    safe = re.sub(r"[^A-Za-z0-9]+", "_", channel).strip("_")
    return flatfield_dir(root) / f"{safe}.npz"


def flatfield_sources(root: Path, map_wells: list[MapWell]
                      ) -> tuple[dict[str, tuple[Path, Path]], list[str]]:
    """Per channel, the (bright, dim) blank stacks a flat field needs.

    Reuses the platemap's existing blank-well roles: the DMEM well (role
    `intensity`) is the bright one and the FluoroBrite well (role
    `background`) the dim one. `verify_map_roles` has already swapped them if
    the measured means say the platemap had them the wrong way round, so this
    reads whatever survived that check rather than trusting the labels.

    A channel with only one blank declared gets no flat field, and its
    positions are corrected for background only - which is the larger of the
    two errors anyway.
    """
    pairs, warnings = resolve_correction_maps(root, map_wells)
    pairs, checks = verify_map_roles(pairs)

    by_channel: dict[str, dict[str, Path]] = {}
    for pair in pairs:
        by_channel.setdefault(pair.channel, {})[pair.role] = pair.source

    sources, notes = {}, list(warnings)
    for channel, roles in sorted(by_channel.items()):
        bright, dim = roles.get("intensity"), roles.get("background")
        if bright and dim:
            sources[channel] = (bright, dim)
        else:
            have = "dmem" if bright else "fluorobrite"
            notes.append(f"{channel}: only the {have} blank well is declared, so "
                         f"no flat field can be built - positions will be "
                         f"corrected for background only")
    return sources, notes


def flatfield_status(root: Path, channel: str) -> str:
    return "built" if flatfield_path(root, channel).exists() else "pending"


def build_flatfield(root: Path, channel: str, bright: Path, dim: Path,
                    frames: int = 12) -> dict:
    """Build and cache one channel's flat field. Returns its provenance."""
    import numpy as np
    import tifffile
    import signal_correction as sc

    def sample(path: Path):
        with tifffile.TiffFile(path) as fh:
            n = fh.series[0].shape[0] if fh.series[0].ndim == 3 else 1
            step = max(1, n // frames)
            return tifffile.imread(path, key=range(0, n, step))

    flat = sc.flatfield_from_blank_pair(sample(bright), sample(dim))
    meta = {
        "channel": channel,
        "bright_source": bright.name,
        "dim_source": dim.name,
        "centre_edge": round(sc.centre_edge_ratio(flat), 4),
        "min": round(float(flat.min()), 4),
        "max": round(float(flat.max()), 4),
        "built": now(),
    }
    path = flatfield_path(root, channel)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, flatfield=flat, meta=json.dumps(meta))
    return meta


def load_flatfield(root: Path, channel: str):
    """The cached flat field and its provenance, or (None, None)."""
    import numpy as np

    path = flatfield_path(root, channel)
    if not path.exists():
        return None, None
    with np.load(path, allow_pickle=False) as z:
        return z["flatfield"], json.loads(str(z["meta"]))


def cmd_flatfield(args) -> int:
    """Build the plate's flat fields from the blank wells named in the platemap."""
    root = Path(args.root).resolve()
    mapfile = Path(args.map).resolve() if getattr(args, "map", None) else platemap_path(root)
    rows = read_platemap(mapfile)
    positions = discover_positions(root, args.pattern)
    _, map_wells, _ = build_tasks(root, positions, rows)

    sources, notes = flatfield_sources(root, map_wells)
    for note in dict.fromkeys(notes):
        log(f"  ! {note}")
    if not sources:
        log("No channel has both blank wells declared in the platemap, so no "
            "flat field is built. Every position will be corrected for "
            "background only.")
        return 0

    todo = [c for c in sources
            if args.force or flatfield_status(root, c) == "pending"]
    if not todo:
        log("Every flat field is already built.")
        for channel in sources:
            _, meta = load_flatfield(root, channel)
            log(f"  {channel:<10} centre/edge {meta['centre_edge']} "
                f"from {meta['bright_source'][:40]}")
        return 0

    tee = tee_to(root, "flatfield", "flatfield")
    saved, sys.stdout = sys.stdout, tee
    failures = 0
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        for channel in todo:
            bright, dim = sources[channel]
            started = time.time()
            try:
                log(f"[{now()}] {channel}: bright {bright.name}")
                log(f"           dim    {dim.name}")
                meta = build_flatfield(root, channel, bright, dim)
                log(f"  centre/edge {meta['centre_edge']}, range "
                    f"{meta['min']}-{meta['max']}, {time.time()-started:.0f}s")
            except Exception as exc:
                failures += 1
                log(f"  ! {channel}: {exc}")
                log(traceback.format_exc())
    finally:
        sys.stdout = saved
        tee.flush()

    built = [c for c in sources if flatfield_status(root, c) == "built"]
    ratios = {}
    for channel in built:
        _, meta = load_flatfield(root, channel)
        ratios[channel] = meta["centre_edge"]
        log(f"{channel:<10} centre/edge {meta['centre_edge']}")
    # The flat field is a property of the optics, so every channel should
    # measure the same vignette. Two channels disagreeing means one of the
    # blank pairs is not what the platemap says it is.
    if len(ratios) > 1:
        spread = max(ratios.values()) - min(ratios.values())
        if spread > 0.05:
            log(f"  ! the channels disagree on the vignette by {spread:.3f} "
                f"({ratios}); a flat field is optics, so they should agree to "
                f"~0.01. Check that both blank wells are what the platemap says.")
        else:
            log(f"  channels agree on the vignette to {spread:.3f} - good")
    if built:
        log("")
        log("Positions analyzed before this was built keep their old numbers; "
            "re-run them with  pipeline.py submit --force-analysis")
    return 1 if failures else 0


def maps_ledger_path(root: Path) -> Path:
    return pipeline_dir(root) / "state" / "maps" / "corrections.json"


def read_maps_ledger(root: Path) -> dict:
    path = maps_ledger_path(root)
    if not path.exists():
        return {}
    try:
        with open(path) as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError):
        return {}


def record_map(root: Path, pair: MapPair, started: float,
               check: RoleCheck | None = None) -> None:
    ledger = read_maps_ledger(root)
    ledger[pair.key] = {
        "channel": pair.channel,
        "role": pair.role,
        "source": pair.source.name,
        "source_well": (STUB_RE.search(pair.source.name).group()
                        if STUB_RE.search(pair.source.name) else ""),
        "output": pair.output.name,
        "finished": now(),
        "duration_s": round(time.time() - started, 1),
        "slurm_job": os.environ.get("SLURM_JOB_ID", ""),
        # Provenance for the summary sheet: which well this really came from,
        # and whether the platemap had the two blank wells the right way round.
        "role_check": check.as_dict() if check else {"verdict": "not-checked"},
    }
    path = maps_ledger_path(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    with open(tmp, "w") as fh:
        json.dump(ledger, fh, indent=2, sort_keys=True)
    os.replace(tmp, path)


def filesystem_of(path: Path) -> str:
    """Mount point and type a path sits on, from /proc/mounts.

    Worth a line in the log: these stacks are gigabytes, and whether the root
    folder is on a parallel scratch filesystem or on shared NFS is the
    difference between a map job that takes a minute and one that takes an
    hour. Empty string off Linux.
    """
    try:
        with open("/proc/mounts") as fh:
            mounts = [(parts[1], parts[2]) for parts in
                      (line.split() for line in fh) if len(parts) >= 3]
    except OSError:
        return ""
    target = str(path.resolve())
    matching = [m for m in mounts if target == m[0] or target.startswith(m[0].rstrip("/") + "/")]
    if not matching:
        return ""
    point, kind = max(matching, key=lambda m: len(m[0]))
    return f"{point} ({kind})"


def map_status(pair: MapPair) -> str:
    """done or pending. The map file itself is the completion record - it is
    named after its source stack, so asking for a map from a different well
    asks for a different file and this correctly reports pending."""
    return ("done" if pair.output.exists() and pair.output.stat().st_size > 0
            else "pending")


# ---------------------------------------------------------------------------
# Logging that works the same whether launched by SLURM or by hand
# ---------------------------------------------------------------------------

class Tee:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.fh = open(path, "a", buffering=1)
        self.stdout = sys.stdout

    def write(self, data):
        self.stdout.write(data)
        self.fh.write(data)
        return len(data)

    def flush(self):
        self.stdout.flush()
        self.fh.flush()


def tee_to(root: Path, stage: str, stem: str):
    return Tee(pipeline_dir(root) / "logs" / stage / f"{stem}.log")


# ---------------------------------------------------------------------------
# Stage 1: inference (needs the cellaap/detectron2 environment and a GPU)
# ---------------------------------------------------------------------------

def run_inference_one(task: Task, container, force: bool = False) -> str:
    """Segment one position. Returns 'done' or 'skipped'.

    Writes into `<inference_dir>.partial` and renames only after every stack
    is on disk, so an interrupted run leaves nothing that looks finished.
    """
    import numpy as np
    import tifffile
    import inference as inf  # type: ignore

    out = task.inference_dir
    if inference_outputs_present(task) and not force:
        log(f"  inference already present at {out.name}; skipping")
        return "skipped"

    partial = out.with_name(out.name + ".partial")
    if partial.exists():
        log(f"  clearing stale partial output {partial.name}")
        shutil.rmtree(partial)
    partial.mkdir(parents=True)

    with tifffile.TiffFile(task.pos.phs) as phs:
        n_pages = len(phs.pages)
    interval = [0, n_pages - 1]
    log(f"  {task.pos.phs.name}: {n_pages} frames, model {task.model}, "
        f"confluency {task.confluency_est}, threshold {task.conf_threshold}")

    result = inf.run_inference(container, task.pos.phs, interval)

    written = []
    for kind, array, dtype in (
        ("semantic", result["semantic_movie"], np.uint8),
        ("instance", result["instance_movie"], np.uint16),
        ("scores", result["scores_movie"], np.uint16),
    ):
        path = partial / f"{task.pos.stem}_{kind}.tif"
        tifffile.imwrite(path, array.astype(dtype))
        written.append(path.name)

    frames_out = int(np.asarray(result["semantic_movie"]).shape[0])
    if frames_out != n_pages:
        shutil.rmtree(partial)
        raise RuntimeError(f"inference returned {frames_out} frames for a "
                           f"{n_pages}-frame stack; output discarded")

    if out.exists():
        if not force:
            shutil.rmtree(partial)
            raise RuntimeError(f"{out} appeared while this job was running; "
                               f"not overwriting")
        shutil.rmtree(out)
    partial.rename(out)
    log(f"  wrote {out.name}: {', '.join(written)}")
    return "done"


def cmd_infer(args) -> int:
    root = Path(args.root).resolve()
    tasks = load_tasks_for_run(root, args)
    tasks = select_tasks(root, tasks, "inference", args)
    if not tasks:
        log("Nothing to infer.")
        return 0

    import inference as inf  # type: ignore

    # One model load serves every position in this process; on an array job
    # that is one position, running locally it is the whole plate.
    containers: dict[tuple, object] = {}
    failures = 0
    for i, task in enumerate(tasks, start=1):
        tee = tee_to(root, "inference", task.stem)
        saved, sys.stdout = sys.stdout, tee
        started = time.time()
        try:
            log(f"[{now()}] inference {i}/{len(tasks)}: {task.stem} "
                f"({task.celltype}, well {task.pos.well})")
            key = (task.model, task.confluency_est, task.conf_threshold)
            if key not in containers:
                log(f"  configuring {task.model} ...")
                containers[key] = inf.configure(task.model, task.confluency_est,
                                                task.conf_threshold)
            status = run_inference_one(task, containers[key], force=args.force)
            write_state(root, "inference", task, "done", started,
                        outputs=[task.inference_dir.name])
            log(f"[{now()}] {status} in {time.time() - started:.0f}s")
        except Exception:
            failures += 1
            err = traceback.format_exc()
            log(err)
            write_state(root, "inference", task, "failed", started, error=err)
        finally:
            sys.stdout = saved
            tee.flush()

    if failures:
        log(f"{failures}/{len(tasks)} position(s) failed inference.")
    return 1 if failures else 0


# ---------------------------------------------------------------------------
# Stage 2: analysis (img-env, CPU)
# ---------------------------------------------------------------------------

def write_corrections_sheet(session, task: Task, root: Path,
                            corrections: list[dict] | None = None) -> None:
    """Append a `corrections` sheet to this position's summary workbook.

    The summary is the file that outlives the run, so the provenance of the
    correction applied to its numbers belongs in it. Two things are recorded,
    because two corrections can be present at once:

    * what `signal_correction` did for this position - the background it
      measured, how far it drifted, how close to cells it had to work, and
      which blank pair the flat field came from. This is the correction behind
      the `<ch>_corrected` columns, and the one to read.
    * whether the legacy `*_map.tif` files were also on disk and loaded by
      cellaap_analysis, which is what the `<ch>_bkg_corr` / `<ch>_int_corr`
      columns carry. The pipeline no longer builds those, but a folder from an
      earlier run still has them, and a reader months later needs to know that
      the two sets of numbers are not the same correction.
    """
    import pandas as pd

    rows = []
    for record in (corrections or []):
        rows.append({
            "correction": "signal_correction (per position)",
            "channel": record.get("channel", ""),
            "applied": "yes",
            "flatfield_from": record.get("flatfield", "") or "none - background only",
            "flatfield_centre_edge": record.get("flatfield_centre_edge", ""),
            "background_counts": record.get("background_mean_range", ""),
            "background_drift_percent": round(
                float(record.get("background_drift_percent", 0)), 1),
            "dilation_used_px": record.get("dilation_used", ""),
            "unusable_block_fraction": round(
                float(record.get("unusable_block_fraction", 0)), 3),
            "frames_sampled": record.get("n_frames_sampled", ""),
            "background_model": record.get("background_model", ""),
            "model_switch_reason": record.get("background_model_reason", ""),
            "surface_file": record.get("surface_file", ""),
        })
    if not rows:
        rows.append({"correction": "signal_correction (per position)",
                     "channel": "", "applied": "no",
                     "flatfield_from": "", "flatfield_centre_edge": "",
                     "background_counts": "", "background_drift_percent": "",
                     "dilation_used_px": "", "unusable_block_fraction": "",
                     "frames_sampled": "", "background_model": "",
                     "model_switch_reason": "", "surface_file": ""})

    legacy = {"background": session.background_map_present,
              "intensity": session.intensity_map_present}
    for role, present in legacy.items():
        rows.append({
            "correction": f"legacy {role} map (not built by this pipeline)",
            "channel": "", "applied": "yes" if present else "no",
            "flatfield_from": "", "flatfield_centre_edge": "",
            "background_counts": "", "background_drift_percent": "",
            "dilation_used_px": "", "unusable_block_fraction": "",
            "frames_sampled": "", "background_model": "",
            "model_switch_reason": "", "surface_file": "",
        })
    frame = pd.DataFrame(rows)

    note = ("`<ch>_corrected` is signal_correction's number: background measured "
            "from this position's own frames, illumination from the difference "
            "of the two blank wells. `<ch>` is raw. `<ch>_bkg_corr` and "
            "`<ch>_int_corr` come from the legacy *_map.tif files if any were "
            "still in the folder, and are a DIFFERENT correction - do not mix "
            "them with `<ch>_corrected`. See SIGNAL_CORRECTION_README.md. "
            "`surface_file` is the background this position actually "
            "subtracted, saved at grid resolution; read it with "
            "signal_correction.read_background_stack().")

    for summary in sorted(task.inference_dir.glob("*_summary*.xlsx")):
        try:
            # Replace on a re-run so the sheet never accumulates stale rows,
            # then append the note through openpyxl - a second to_excel call
            # on the same sheet would replace the table rather than follow it.
            with pd.ExcelWriter(summary, mode="a", engine="openpyxl",
                                if_sheet_exists="replace") as writer:
                frame.to_excel(writer, sheet_name="corrections", index=False)
            from openpyxl import load_workbook
            book = load_workbook(summary)
            sheet = book["corrections"]
            sheet.append([])
            sheet.append([note])
            book.save(summary)
        except Exception as exc:
            # The summary itself is the deliverable; failing to annotate it is
            # not a reason to fail the position.
            log(f"  ! could not add the corrections sheet to {summary.name} ({exc})")


def apply_signal_correction(session, task: Task, args) -> list[dict]:
    """Add `<ch>_corrected` to this position's per-frame table.

    Runs after every channel has been measured and before summarize_data, so
    the corrected signal is in `cell_data` - where the per-frame dynamics can
    be recovered from it - and gets averaged into the summary along with
    everything else.

    The background is estimated from this position's own frames. The flat
    field is the plate's, built earlier by `pipeline.py flatfield`; when the
    platemap declares no blank pair for a channel there is none, and that
    channel is corrected for background only.
    """
    import numpy as np
    import tifffile
    import signal_correction as sc

    records = []
    labels = session.stacks.get("instance")
    if labels is None:
        log("  ! no instance stack in memory; signal correction skipped")
        return records

    for channel in task.channels:
        if channel not in session.paths or channel not in session.tracked.columns:
            continue
        flat, flat_meta = load_flatfield(task.root, channel)

        with tifffile.TiffFile(session.paths[channel]) as fh:
            stack = fh.series[0].asarray(out="memmap")
            correction = sc.estimate_position_correction(
                stack, stem=task.stem, channel=channel, labels=labels,
                flatfield=flat, n_frames=args.correction_frames,
                block=args.correction_block, dilation=args.correction_dilation,
                background_model=args.correction_model)

            # Centroids are at segmentation scale; the fluorescence frame is
            # bigger by whatever the inference downsampling was. Derive the
            # factor from the two shapes rather than assuming the 2x that
            # cellaap currently uses - a different binning would otherwise
            # sample the correction from the wrong part of the field, silently.
            frame_shape = stack.shape[-2:]
        scale_y = frame_shape[0] / labels.shape[-2]
        scale_x = frame_shape[1] / labels.shape[-1]

        session.tracked[f"{channel}_corrected"] = correction.correct_measurements(
            session.tracked[channel].to_numpy(dtype=float),
            session.tracked["x"].to_numpy(dtype=float) * scale_x,
            session.tracked["y"].to_numpy(dtype=float) * scale_y,
            session.tracked["frame"].to_numpy(dtype=int))

        # Keep the surface that was subtracted, at grid resolution: ~490 kB
        # per position-channel, against ~1.1 GB for the stack it came from.
        # It is redundant with the correction object, which reconstructs it
        # exactly, but a TIFF opens in Fiji and an npz does not.
        try:
            written = sc.save_background_stack(
                correction, surface_path(task.root, task.stem, channel))
        except Exception as exc:
            written = None
            log(f"  ! could not save the {channel} background surface ({exc})")

        d = correction.diagnostics
        measured = session.tracked[f"{channel}_corrected"].notna().sum()
        fell_back = d["background_model"] != d.get("background_model_requested")
        log(f"  corrected {channel}: background "
            f"{d['background_mean_range'][0]:.1f}-{d['background_mean_range'][1]:.1f} "
            f"counts (drift {d['background_drift_percent']:.0f}%), "
            f"dilation {d['dilation_used']}, "
            f"{100*d['unusable_block_fraction']:.0f}% blocks unusable, "
            f"flat field {'yes' if flat is not None else 'NO - background only'}, "
            f"model {d['background_model']}, {measured} rows")
        if written:
            log(f"    surface saved to {written.parent.name}/{written.name} "
                f"({written.stat().st_size/1024:.0f} kB)")
        if d.get("background_model_reason"):
            log(f"    ! switched to the {d['background_model']} background model: "
                f"{d['background_model_reason']}")
        records.append({
            "channel": channel,
            "scale": round(scale_x, 3),
            "flatfield": flat_meta["bright_source"] if flat_meta else "",
            "flatfield_dim_source": flat_meta["dim_source"] if flat_meta else "",
            "flatfield_centre_edge": flat_meta["centre_edge"] if flat_meta else "",
            "surface_file": str(written.relative_to(task.root)) if written else "",
            **{k: (json.dumps(v) if isinstance(v, list) else v)
               for k, v in d.items()},
        })

    if records:
        session.write_analysis_file()
    return records


def run_analysis_one(session, task: Task, semantic_gap: int | None,
                     args=None) -> None:
    """Track, measure, correct and summarize one inference folder."""
    import numpy as np

    if not inference_outputs_present(task):
        raise RuntimeError(f"no inference output at {task.inference_dir}")

    session.files(task.inference_dir, cell_type=task.analysis_cell_type)
    if semantic_gap:
        session.defaults.semantic_gap_closing = semantic_gap
        session.defaults.semantic_footprint = np.ones(semantic_gap)

    # The session object is reused across positions when this runs locally,
    # and summarize_data/measure_signal would otherwise see the previous
    # position's frames.
    for attr in ("summaryDF", "tracked"):
        if hasattr(session, attr):
            delattr(session, attr)

    session.track_centroids(session.defaults.track_mode, save_flag=True)
    for channel in task.channels:
        if channel not in session.paths:
            log(f"  ! no {channel} stack for this position; not measured")
            continue
        log(f"  measuring {channel} ...")
        session.measure_signal(channel, True, -1)

    # Correction before the summary, deliberately. summarize_data averages
    # whatever per-channel columns it finds over each track's window, so
    # `<ch>_corrected` only reaches the summary if it is on the table first -
    # and analysis_outputs_present keys on the summary existing, so a position
    # whose correction failed reports pending rather than finished-but-wrong.
    corrections = apply_signal_correction(session, task, args) if args else []

    session.summarize_data(True)
    write_corrections_sheet(session, task, task.root, corrections)


def cmd_analyze(args) -> int:
    root = Path(args.root).resolve()
    tasks = load_tasks_for_run(root, args)
    tasks = select_tasks(root, tasks, "analysis", args)
    if not tasks:
        log("Nothing to analyze.")
        return 0

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import cellaap_analysis  # type: ignore

    session = cellaap_analysis.analysis(root, plotting_only=False)
    # Say it once, up front: which corrections these numbers carry is the kind
    # of thing that has to be recoverable from the log months later.
    built = [c for c in ("GFP", "Texas Red", "Cy5")
             if flatfield_status(root, c) == "built"]
    log(f"signal_correction: background measured per position; flat field "
        f"{'for ' + ', '.join(built) if built else 'NOT built - background only'}")
    if session.background_map_present or session.intensity_map_present:
        log(f"  legacy *_map.tif files are also on disk and loaded "
            f"(background {'yes' if session.background_map_present else 'no'}, "
            f"intensity {'yes' if session.intensity_map_present else 'no'}); "
            f"they fill <ch>_bkg_corr/<ch>_int_corr, NOT <ch>_corrected")
    failures = 0
    for i, task in enumerate(tasks, start=1):
        tee = tee_to(root, "analysis", task.stem)
        saved, sys.stdout = sys.stdout, tee
        started = time.time()
        try:
            log(f"[{now()}] analysis {i}/{len(tasks)}: {task.stem} "
                f"(cell_type={task.analysis_cell_type}, "
                f"channels={task.channels or 'none'})")
            run_analysis_one(session, task, args.semantic_gap, args)
            outputs = sorted(p.name for p in task.inference_dir.glob("*.xlsx"))
            write_state(root, "analysis", task, "done", started, outputs=outputs)
            log(f"[{now()}] done in {time.time() - started:.0f}s")
        except Exception:
            failures += 1
            err = traceback.format_exc()
            log(err)
            write_state(root, "analysis", task, "failed", started, error=err)
        finally:
            sys.stdout = saved
            tee.flush()

    if failures:
        log(f"{failures}/{len(tasks)} position(s) failed analysis.")
    return 1 if failures else 0


# ---------------------------------------------------------------------------
# Task selection shared by both stages
# ---------------------------------------------------------------------------

def load_tasks_for_run(root: Path, args) -> list[Task]:
    """Tasks for a run, from the frozen task list when one was given.

    `--tasks` is what makes an array job reproducible: the file fixes both the
    set of positions and their order, so index N means the same position in
    the inference array and in the analysis array that depends on it.
    """
    mapfile = Path(args.map).resolve() if getattr(args, "map", None) else platemap_path(root)
    rows = read_platemap(mapfile)
    positions = discover_positions(root, args.pattern)
    tasks, map_wells, warnings = build_tasks(
        root, positions, rows, correction=correction_params(args))

    if getattr(args, "tasks", None):
        wanted = [ln.strip() for ln in Path(args.tasks).read_text().splitlines()
                  if ln.strip() and not ln.startswith("#")]
        by_stem = {t.stem: t for t in tasks}
        missing = [s for s in wanted if s not in by_stem]
        if missing:
            die(f"task list {args.tasks} names positions that the platemap no "
                f"longer covers: {', '.join(missing[:5])}")
        tasks = [by_stem[s] for s in wanted]
    else:
        for w in warnings:
            log(f"  ! {w}")
    return tasks


def select_tasks(root: Path, tasks: list[Task], stage: str, args) -> list[Task]:
    """Apply --index / --stem filters, then drop finished work."""
    if getattr(args, "index", None) is not None:
        if not 0 <= args.index < len(tasks):
            die(f"--index {args.index} is out of range (0..{len(tasks) - 1})")
        tasks = [tasks[args.index]]
    if getattr(args, "stem", None):
        wanted = set(args.stem)
        tasks = [t for t in tasks if t.stem in wanted or t.pos.well in wanted
                 or t.pos.stub in wanted]
        if not tasks:
            die("no position matched --stem")

    if getattr(args, "force", False):
        return tasks

    keep = []
    for task in tasks:
        status = stage_status(root, task, stage)
        if status == "done":
            log(f"  {task.stem}: {stage} already done, skipping")
            continue
        if status == "stale":
            log(f"  {task.stem}: {stage} was run with different parameters, re-running")
        keep.append(task)
    return keep


# ---------------------------------------------------------------------------
# init / check / status
# ---------------------------------------------------------------------------

def cmd_init(args) -> int:
    root = Path(args.root).resolve()
    if not root.is_dir():
        die(f"{root} is not a directory")
    positions = discover_positions(root, args.pattern)
    if not positions:
        die(f"no files matching {args.pattern} in {root}")

    path = write_template(root, positions, args.force)
    wells = sorted({p.well for p in positions})
    log(f"{len(positions)} position(s) in {len(wells)} well(s): {', '.join(wells)}")
    log(f"Wrote {path}")
    log("")
    log(f"Fill in celltype/transfection/drug (celltype defaults to {DEFAULT_CELLTYPE}),")
    log(f"then run:  python pipeline.py check --root {root}")
    return 0


def cmd_check(args) -> int:
    root = Path(args.root).resolve()
    rows = read_platemap(platemap_path(root))
    positions = discover_positions(root, args.pattern)
    if not positions:
        die(f"no files matching {args.pattern} in {root}")
    tasks, map_wells, warnings = build_tasks(
        root, positions, rows, correction=correction_params(args))

    problems: list[str] = []
    for task in tasks:
        if not task.pos.phs.exists():
            problems.append(f"{task.stem}: phase stack missing")
        if not task.channels:
            warnings.append(f"{task.stem}: no fluorescence channel found; "
                            f"only mitotic timing will be reported")

    # A readable phase stack with a plausible frame count is worth the few
    # seconds it costs here: the alternative is finding out on the GPU node.
    try:
        import tifffile
        for task in tasks:
            try:
                with tifffile.TiffFile(task.pos.phs) as fh:
                    n = len(fh.pages)
                if n < 2:
                    problems.append(f"{task.stem}: phase stack has {n} page(s)")
            except Exception as exc:
                problems.append(f"{task.stem}: cannot read phase stack ({exc})")
    except ImportError:
        warnings.append("tifffile not importable here; phase stacks not verified")

    pairs, map_warnings = resolve_correction_maps(root, map_wells)
    warnings.extend(map_warnings)
    pairs, role_checks = verify_map_roles(pairs)
    ff_sources, ff_notes = flatfield_sources(root, map_wells)
    warnings.extend(ff_notes)
    # Legacy *_map.tif files are no longer built, but cellaap_analysis still
    # loads any it finds, so a folder from an earlier run needs saying so.
    stale_maps = [p.name for role in MAP_ROLES for channel in
                  {q.channel for q in pairs}
                  for p in find_map_files(root, channel, role)]
    if stale_maps:
        warnings.append(
            f"{len(stale_maps)} legacy *_map.tif file(s) are still in the "
            f"folder and will be loaded into <ch>_bkg_corr/<ch>_int_corr. "
            f"They are a different correction from <ch>_corrected; the "
            f"pipeline no longer builds them")

    log(f"Root:      {root}")
    log(f"Platemap:  {platemap_path(root)}")
    log(f"Positions: {len(tasks)} to run "
        f"({len(positions) - len(tasks) - len(map_wells)} skipped by the "
        f"platemap, {len(map_wells)} blank-media)")
    log("")

    header = f"{'position':<46} {'well':<5} {'celltype':<9} {'model':<14} {'channels':<20} {'infer':<8} {'analysis'}"
    log(header)
    log("-" * len(header))
    counts = {"inference": {}, "analysis": {}}
    for task in tasks:
        s_inf = stage_status(root, task, "inference")
        s_ana = stage_status(root, task, "analysis")
        counts["inference"][s_inf] = counts["inference"].get(s_inf, 0) + 1
        counts["analysis"][s_ana] = counts["analysis"].get(s_ana, 0) + 1
        flag = "" if task.mapped else "*"
        log(f"{task.stem[:45]:<46} {task.pos.well:<5} {task.celltype + flag:<9} "
            f"{task.model:<14} {';'.join(task.channels)[:19]:<20} {s_inf:<8} {s_ana}")

    log("")
    log(f"inference: {fmt_counts(counts['inference'])}")
    log(f"analysis:  {fmt_counts(counts['analysis'])}")
    if any(not t.mapped for t in tasks):
        log(f"* well not in the platemap; defaulted to {DEFAULT_CELLTYPE}")

    log("")
    log("Signal correction:")
    log(f"  background   measured per position from its own frames "
        f"(dilation {args.correction_dilation} px, "
        f"{args.correction_frames} frames) - always applied")
    if map_wells:
        for line in report_role_checks(role_checks):
            log(line)
        if any(c.verdict == "swapped" for c in role_checks):
            log("  the sources below are AFTER that swap")
        if ff_sources:
            for channel, (bright, dim) in sorted(ff_sources.items()):
                state = flatfield_status(root, channel)
                log(f"  flat field   {channel:<10} {state:<8} from "
                    f"{bright.name[:34]:<35} - {dim.name[:34]}")
        else:
            log("  flat field   no channel has BOTH blank wells declared, so "
                "illumination is not corrected")
    else:
        log("  flat field   no blank wells in the platemap; positions get the "
            "background correction only")

    pending_inf = counts["inference"].get("pending", 0) + counts["inference"].get("stale", 0)
    pending_ana = counts["analysis"].get("pending", 0) + counts["analysis"].get("stale", 0)
    log("")
    if not pending_inf and pending_ana:
        log("Every position already has an inference folder: submit will run "
            "the analysis only, with no GPU job and no wait on the GPU queue.")
        log("")
    log(f"Rough wall time at ~35 min/inference and ~45 min/analysis: "
        f"{pending_inf * 35 / 60:.1f} GPU-hours, {pending_ana * 45 / 60:.1f} CPU-hours "
        f"(array jobs run these in parallel).")

    if warnings:
        log("")
        log("Warnings:")
        for w in dict.fromkeys(warnings):
            log(f"  ! {w}")
    if problems:
        log("")
        log("Problems that will stop the run:")
        for p in problems:
            log(f"  x {p}")
        return 1
    log("")
    log("Platemap is consistent with the folder. Ready to submit.")
    return 0


def fmt_counts(counts: dict) -> str:
    order = ["done", "pending", "stale", "failed"]
    parts = [f"{k} {counts[k]}" for k in order if counts.get(k)]
    return ", ".join(parts) if parts else "nothing to do"


def cmd_status(args) -> int:
    root = Path(args.root).resolve()
    rows = read_platemap(platemap_path(root))
    positions = discover_positions(root, args.pattern)
    tasks, map_wells, _ = build_tasks(
        root, positions, rows, correction=correction_params(args))

    records = []
    for task in tasks:
        rec = {
            "position": task.stem, "well": task.pos.well, "site": task.pos.site,
            "celltype": task.celltype, "transfection": task.transfection,
            "drug": task.drug, "model": task.model,
            "channels": ";".join(task.channels),
            "inference": stage_status(root, task, "inference"),
            "analysis": stage_status(root, task, "analysis"),
            "inference_dir": task.inference_dir.name,
        }
        for stage in ("inference", "analysis"):
            st = read_state(root, stage, task.stem) or {}
            rec[f"{stage}_duration_s"] = st.get("duration_s", "")
            rec[f"{stage}_finished"] = st.get("finished", "")
        records.append(rec)

    counts = {"inference": {}, "analysis": {}}
    for rec in records:
        for stage in ("inference", "analysis"):
            counts[stage][rec[stage]] = counts[stage].get(rec[stage], 0) + 1

    ff_sources, _ = flatfield_sources(root, map_wells)
    ff_counts = {}
    for channel in ff_sources:
        s = flatfield_status(root, channel)
        ff_counts[s] = ff_counts.get(s, 0) + 1
    if ff_sources:
        log(f"  flat field: {fmt_counts(ff_counts)} "
            f"({', '.join(sorted(ff_sources))})")
    else:
        log(f"  flat field: no blank pair declared; background correction only")

    failed = [(r, s) for r in records for s in ("inference", "analysis")
              if r[s] == "failed"]
    if failed:
        log("")
        log("Failures (last line of the traceback; full log under pipeline/logs/):")
        for rec, stage in failed:
            st = read_state(root, stage, rec["position"]) or {}
            last = (st.get("error", "").strip().splitlines() or ["?"])[-1]
            log(f"  {stage:<9} {rec['position']}")
            log(f"            {last[:150]}")

    if args.csv:
        out = Path(args.csv) if args.csv != "-" else pipeline_dir(root) / "status.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(records[0]) if records else ["position"])
            writer.writeheader()
            writer.writerows(records)
        log("")
        log(f"Wrote {out}")
    return 0


# ---------------------------------------------------------------------------
# Report: what a run actually did, read back out of the logs
#
# `status` answers "what is finished?" by looking at state markers and outputs.
# That misses the one failure the markers cannot record: a task killed by SLURM
# for exceeding its walltime is SIGKILLed, so it never writes a marker, and the
# position reads `pending` - indistinguishable from one that was never
# submitted. The only record is the job's .out file. This command reads those,
# joins them back to positions through the frozen task lists, and says plainly
# which stacks were cut off and what to do about it.
# ---------------------------------------------------------------------------

# How SLURM ends a task, in the words it uses in the .out file. Order matters:
# a walltime kill is reported as a cancellation, so TIME LIMIT has to be tested
# before the generic CANCELLED or every timeout would be filed as a manual
# scancel.
SLURM_SIGNATURES = (
    ("timeout",   re.compile(r"DUE TO TIME LIMIT", re.I)),
    ("oom",       re.compile(r"oom[-_ ]kill|Out Of Memory|Exceeded job memory", re.I)),
    ("preempted", re.compile(r"DUE TO PREEMPTION", re.I)),
    ("node_fail", re.compile(r"NODE_FAIL|DUE TO NODE FAILURE", re.I)),
    ("cancelled", re.compile(r"CANCELLED AT .*\*\*\*", re.I)),
    ("env_error", re.compile(r"could not activate .*conda environment", re.I)),
    ("traceback", re.compile(r"^Traceback \(most recent call last\)", re.M)),
)

# What each kind means for the user, and whether the work is recoverable by
# simply asking for more of something.
SIGNATURE_HELP = {
    "timeout":   ("killed at the walltime limit", "raise the time limit"),
    "oom":       ("killed for exceeding its memory request", "raise the memory"),
    "preempted": ("preempted off the node", "just resubmit"),
    "node_fail": ("lost to a node failure", "just resubmit"),
    "cancelled": ("cancelled (scancel, or a failed dependency)", "resubmit if unintended"),
    "env_error": ("died before python started - conda activate failed", "fix the environment"),
    "traceback": ("raised a python exception", "read the traceback"),
}

CANCEL_AT_RE = re.compile(r"CANCELLED AT (\S+)")
ARRAY_OUT_RE = re.compile(r"^(?P<name>.+)_(?P<jobid>\d+)_(?P<index>\d+)\.out$")
PLAIN_OUT_RE = re.compile(r"^(?P<name>.+)_(?P<jobid>\d+)\.out$")
LOG_STAMP_RE = re.compile(r"^\[(\d{4}-\d{2}-\d{2}T[\d:]+)\]", re.M)
# The line every stage writes when it picks a position up. Position logs
# are APPENDED across attempts, so a log can hold several of these and only
# the last one before a given cancellation belongs to that run.
LOG_START_RE = re.compile(
    r"^\[(\d{4}-\d{2}-\d{2}T[\d:]+)\] (?:inference|analysis) \d+/\d+:", re.M)
SBATCH_DIRECTIVE_RE = re.compile(r"^#SBATCH --(?P<key>[a-z-]+)=(?P<value>.+)$", re.M)

# job-name suffix -> (stage, frozen task list). submit builds these names, so
# they are the link from a .out file back to a position when jobids.json is
# absent (scripts written but submitted by hand).
JOBNAME_SUFFIXES = (
    ("_ana0", "analysis", "tasks_ready.txt"),
    ("_inf", "inference", "tasks.txt"),
    ("_ana", "analysis", "tasks.txt"),
    ("_flat", "flatfield", ""),
)


def parse_walltime(value: str) -> int | None:
    """SLURM walltime to seconds. Accepts D-HH:MM:SS, HH:MM:SS, MM:SS, MM."""
    value = (value or "").strip()
    if not value:
        return None
    days = 0
    if "-" in value:
        head, _, value = value.partition("-")
        try:
            days = int(head)
        except ValueError:
            return None
    parts = value.split(":")
    try:
        nums = [int(p) for p in parts]
    except ValueError:
        return None
    if len(parts) == 3:
        h, m, s = nums
    elif len(parts) == 2:
        h, m, s = 0, nums[0], nums[1]
    elif len(parts) == 1:
        h, m, s = 0, nums[0], 0
    else:
        return None
    return days * 86400 + h * 3600 + m * 60 + s


def fmt_duration(seconds) -> str:
    if seconds is None:
        return "?"
    seconds = int(round(float(seconds)))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m"
    if m:
        return f"{m}m{s:02d}s"
    return f"{s}s"


def job_dirs(root: Path) -> list[Path]:
    """Submit directories, oldest first. The name is a timestamp, so sorting
    by name is sorting by time."""
    base = pipeline_dir(root) / "jobs"
    return sorted((d for d in base.glob("*") if d.is_dir()), key=lambda d: d.name)


def sbatch_request(path: Path) -> dict:
    """The resources a generated sbatch script asked for."""
    if not path.exists():
        return {}
    text = path.read_text(errors="replace")
    got = {m.group("key"): m.group("value").strip()
           for m in SBATCH_DIRECTIVE_RE.finditer(text)}
    return {
        "time": got.get("time", ""),
        "mem": got.get("mem") or got.get("mem-per-gpu", ""),
        "cpus": got.get("cpus-per-task", ""),
    }


def read_tail(path: Path, limit: int = 200_000) -> str:
    """The end of a log. Capped because a chatty position can write megabytes
    and only the end carries the outcome."""
    try:
        size = path.stat().st_size
        with open(path, "rb") as fh:
            if size > limit:
                fh.seek(size - limit)
            return fh.read().decode("utf-8", errors="replace")
    except OSError:
        return ""


def classify_slurm_out(path: Path) -> tuple[str, str]:
    """(kind, evidence line) for one SLURM .out file, or ("", "")."""
    text = read_tail(path)
    if not text:
        return "", ""
    for kind, pattern in SLURM_SIGNATURES:
        m = pattern.search(text)
        if not m:
            continue
        # Quote the line the match landed on - that is the evidence, and it is
        # what the user would have gone looking for by hand.
        start = text.rfind("\n", 0, m.start()) + 1
        end = text.find("\n", m.start())
        line = text[start:end if end != -1 else len(text)].strip()
        return kind, line
    return "", ""


def position_log_span(root: Path, stage: str, stem: str) -> tuple[str, str, bool]:
    """(first stamp, last stamp, looks finished) from a per-position log.

    Every stage opens its log with `[<iso>] <stage> i/N: <stem>` and writes a
    closing `[<iso>] ... in Ns` only on success, so a log with one stamp and no
    closing line is a position that started and never came back.
    """
    path = pipeline_dir(root) / "logs" / stage / f"{stem}.log"
    if not path.exists():
        return [], "", False
    text = read_tail(path)
    stamps = LOG_STAMP_RE.findall(text)
    if not stamps:
        return [], "", False
    opens = list(LOG_START_RE.finditer(text))
    # Only the stretch after the final opening line describes the latest
    # attempt; searching the whole file would let a run that succeeded this
    # morning mark this evening's kill as finished.
    tail = text[opens[-1].end():] if opens else text
    finished = bool(re.search(r"\] .* in \d+s", tail)) or "Traceback" in tail
    return [m.group(1) for m in opens], stamps[-1], finished


def scan_job_dir(root: Path, jobdir: Path, with_logs: bool = True) -> list[dict]:
    """One record per SLURM .out file in a submit directory.

    `with_logs=False` skips reading the per-position logs and the sbatch
    scripts, which is most of the I/O. killed_index only needs the verdict, and
    it runs inside stage_status, so it takes the cheap path.
    """
    logdir = jobdir / "slurm"
    if not logdir.is_dir():
        return []

    try:
        jobids = json.loads((jobdir / "jobids.json").read_text())
    except (OSError, json.JSONDecodeError):
        jobids = {}
    # job id -> stage, when submit actually submitted. Authoritative when
    # present; the job-name suffix is the fallback.
    by_jobid = {str(v): k for k, v in jobids.items()
                if k != "submitted" and v}

    tasks_cache: dict[str, list[str]] = {}

    def task_list(filename: str) -> list[str]:
        if filename not in tasks_cache:
            path = jobdir / filename
            # splitlines(), not split(): a stem can contain a space
            # ("20250213_HT1080 pPS18_A08_s1_phs"), and splitting on whitespace
            # shatters it into two entries and misaligns every array index
            # after it against the wrong position.
            tasks_cache[filename] = ([ln.strip() for ln in
                                      path.read_text().splitlines() if ln.strip()]
                                     if path.exists() else [])
        return tasks_cache[filename]

    records = []
    for out in sorted(logdir.glob("*.out")):
        m = ARRAY_OUT_RE.match(out.name) or PLAIN_OUT_RE.match(out.name)
        if not m:
            continue
        name = m.group("name")
        jobid = m.group("jobid")
        index = int(m.groupdict().get("index") or -1) if "index" in m.groupdict() else -1

        stage, tasks_file = "", ""
        for suffix, st, tf in JOBNAME_SUFFIXES:
            if name.endswith(suffix):
                stage, tasks_file = st, tf
                break
        # jobids.json distinguishes the two analysis arrays, which share a
        # suffix family; prefer it when we have it.
        mapped = by_jobid.get(jobid, "")
        if mapped == "analysis_ready":
            stage, tasks_file = "analysis", "tasks_ready.txt"
        elif mapped in ("inference", "analysis", "flatfield"):
            stage = "inference" if mapped == "inference" else (
                "analysis" if mapped == "analysis" else "flatfield")
            if mapped != "flatfield" and not tasks_file:
                tasks_file = "tasks.txt"

        stem = ""
        if tasks_file and index >= 0:
            stems = task_list(tasks_file)
            if index < len(stems):
                stem = stems[index]

        kind, evidence = classify_slurm_out(out)
        sbatch_name = {"inference": "infer.sbatch", "analysis": "analyze.sbatch",
                       "flatfield": "flatfield.sbatch"}.get(stage, "")
        if stage == "analysis" and tasks_file == "tasks_ready.txt":
            sbatch_name = "analyze_ready.sbatch"

        starts, last, finished = ([], "", False)
        if with_logs and stem and stage in ("inference", "analysis"):
            starts, last, finished = position_log_span(root, stage, stem)

        records.append({
            "job": jobdir.name, "jobid": jobid, "index": index,
            "stage": stage or "?", "position": stem,
            "kind": kind, "evidence": evidence,
            "out": str(out.relative_to(pipeline_dir(root))),
            "request": (sbatch_request(jobdir / sbatch_name)
                        if with_logs and sbatch_name else {}),
            "log_starts": starts,
            "log_started": starts[-1] if starts else "",
            "log_last": last, "log_finished": finished,
        })
    return records


def marker_state(root: Path, stage: str, stem: str) -> str:
    """What the state marker says, ignoring the platemap. `report` must work on
    a folder whose platemap has since been edited or broken, so it cannot go
    through stage_status()."""
    st = read_state(root, stage, stem)
    if st is None:
        return "no marker"
    return st.get("status", "?")


def stage_durations(root: Path, stage: str) -> list[tuple[str, float]]:
    """(position, seconds) for every position whose marker says done."""
    d = pipeline_dir(root) / "state" / stage
    out = []
    for path in sorted(d.glob("*.json")) if d.is_dir() else []:
        try:
            st = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if st.get("status") == "done" and st.get("duration_s") is not None:
            try:
                out.append((st.get("stem", path.stem), float(st["duration_s"])))
            except (TypeError, ValueError):
                continue
    return out


def killed_runtime(record: dict) -> str:
    """How long a walltime-killed task actually ran.

    Measured, not assumed: the position log's opening stamp against the time
    in SLURM's `CANCELLED AT ...` line. Worth the trouble because it separates
    a task that genuinely used its whole limit - raise the limit - from one
    that sat idle on a wedged filesystem and got reaped, where a bigger limit
    changes nothing. Falls back to the requested walltime, which is the upper
    bound.
    """
    asked = parse_walltime(record.get("request", {}).get("time", ""))
    m = CANCEL_AT_RE.search(record.get("evidence", ""))
    starts = record.get("log_starts") or []
    if m and starts:
        try:
            killed_at = datetime.fromisoformat(m.group(1))
            # The run this cancellation ended is the last one that began before
            # it. Taking the first opening line instead reports the age of the
            # log rather than the length of the run - on a position retried
            # three times that reads as nine hours under a two-hour limit.
            began = [datetime.fromisoformat(s) for s in starts]
            began = [b for b in began if b <= killed_at]
            if began:
                measured = (killed_at - began[-1]).total_seconds()
                # A task cannot outlive its own limit. If the arithmetic says
                # otherwise the log and the .out belong to different runs, so
                # report the limit and flag it rather than print a nonsense
                # number with a straight face.
                if asked is None or measured <= asked * 1.1:
                    return fmt_duration(measured)
        except ValueError:
            pass
    return f"{fmt_duration(asked)}?" if asked else "?"


def fmt_bytes(n) -> str:
    """Bytes as GB/MB, or '-' when we simply do not know."""
    if n in (None, "", 0):
        return "-"
    try:
        n = float(n)
    except (TypeError, ValueError):
        return "-"
    for unit, size in (("GB", 1 << 30), ("MB", 1 << 20), ("KB", 1 << 10)):
        if n >= size:
            return f"{n / size:.1f} {unit}"
    return f"{int(n)} B"


def parse_sacct_mem(value: str) -> int | None:
    """sacct MaxRSS ('4194304K', '3.91G', '512M') to bytes."""
    value = (value or "").strip()
    if not value or value in ("", "0"):
        return None
    mult = {"K": 1 << 10, "M": 1 << 20, "G": 1 << 30, "T": 1 << 40}.get(value[-1].upper())
    try:
        return int(float(value[:-1]) * mult) if mult else int(float(value))
    except ValueError:
        return None


def sacct_memory(jobids: set) -> dict:
    """{'<jobid>_<index>': peak bytes} from sacct, or {} if it is not available.

    Only source of per-position memory for runs that happened before markers
    started recording it. sacct reports MaxRSS against the `.batch`/`.extern`
    steps rather than the job row, so the steps are folded back onto their
    parent by taking the largest.
    """
    if not jobids or not shutil.which("sacct"):
        return {}
    try:
        out = subprocess.run(
            ["sacct", "-j", ",".join(sorted(jobids)), "-P", "-n",
             "--format=JobID,MaxRSS"],
            capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return {}
    if out.returncode != 0:
        return {}

    peaks: dict = {}
    for line in out.stdout.splitlines():
        parts = line.split("|")
        if len(parts) < 2:
            continue
        jobid, maxrss = parts[0].strip(), parts[1].strip()
        base = jobid.split(".")[0]          # 123_4.batch -> 123_4
        mem = parse_sacct_mem(maxrss)
        if mem is None:
            continue
        peaks[base] = max(peaks.get(base, 0), mem)
    return peaks


def adopted_inference(root: Path, stem: str) -> bool:
    """Segmentation outputs on disk for this stem, under any parameters.

    Deliberately not `inference_outputs_present`, which needs a Task and so a
    platemap; `report` has to work on a folder whose platemap has since been
    edited, moved or broken.
    """
    for d in root.glob(f"{stem}_*_inference"):
        if d.is_dir() and any(
                p.stat().st_size > 0 for kind in ("semantic", "instance")
                for p in d.glob(f"*{kind}*.tif")):
            return True
    return False


def position_table(root: Path, records: list[dict], use_sacct: bool = False) -> list[dict]:
    """One row per position per stage pair, joining markers to SLURM verdicts.

    The position list is the union of every state marker and every entry in the
    scanned task lists, so a position that was submitted and killed before it
    could write anything still appears - which is the whole point.
    """
    # SLURM verdict per (stage, position), and the job id that carried it.
    verdict: dict = {}
    for r in records:
        if r["position"] and r["stage"] in ("inference", "analysis"):
            key = (r["stage"], r["position"])
            if r["kind"] or key not in verdict:
                verdict[key] = r

    positions: set = set()
    for r in records:
        if r["position"]:
            positions.add(r["position"])
    for stage in ("inference", "analysis"):
        d = pipeline_dir(root) / "state" / stage
        for path in sorted(d.glob("*.json")) if d.is_dir() else []:
            positions.add(path.stem)

    sacct = {}
    if use_sacct:
        ids = {f"{r['jobid']}_{r['index']}" if r["index"] >= 0 else r["jobid"]
               for r in records if r["jobid"]}
        sacct = sacct_memory(ids)

    rows = []
    for stem in sorted(positions):
        row = {"position": stem}
        for stage in ("inference", "analysis"):
            st = read_state(root, stage, stem) or {}
            rec = verdict.get((stage, stem))
            kind = rec["kind"] if rec else ""

            # Marker and SLURM verdict can both exist and disagree - a position
            # killed on one submit and finished on the next carries both. Neither
            # wins by rank; the later one wins, which needs their timestamps.
            kill_at = marker_at = None
            if rec:
                m = CANCEL_AT_RE.search(rec.get("evidence", ""))
                if m:
                    try:
                        kill_at = datetime.fromisoformat(m.group(1))
                    except ValueError:
                        pass
            if st.get("finished"):
                try:
                    marker_at = datetime.fromisoformat(st["finished"])
                except ValueError:
                    pass
            marker_wins = bool(marker_at and (kill_at is None or marker_at >= kill_at))

            if st.get("status") == "done" and marker_wins:
                state, secs = "ok", st.get("duration_s")
            elif kind == "timeout":
                state, secs = "TIMEOUT", None
            elif kind == "oom":
                state, secs = "OOM", None
            elif kind in ("cancelled", "preempted", "node_fail", "env_error"):
                state, secs = kind.upper(), None
            elif st.get("status") == "done":
                state, secs = "ok", st.get("duration_s")
            elif st.get("status") == "failed":
                state, secs = "failed", st.get("duration_s")
            elif rec and rec.get("log_started") and not rec.get("log_finished"):
                state, secs = "running?", None
            elif stage == "inference" and adopted_inference(root, stem):
                # Outputs on disk with no marker: a folder segmented by the old
                # batch scripts, or before markers existed. `status` adopts these
                # as done rather than redoing them, so saying `pending` here
                # would contradict it.
                state, secs = "adopted", None
            else:
                state, secs = "pending", None

            if secs is None and state == "TIMEOUT" and rec:
                measured = killed_runtime(rec)
                secs = parse_walltime(rec["request"].get("time", "")) \
                    if measured.endswith("?") else None
                row[f"{stage}_time"] = measured
            else:
                row[f"{stage}_time"] = fmt_duration(secs) if secs is not None else "-"

            mem = st.get("peak_rss_bytes")
            source = "marker"
            if not mem and rec:
                key = (f"{rec['jobid']}_{rec['index']}" if rec["index"] >= 0
                       else rec["jobid"])
                mem = sacct.get(key)
                source = "sacct" if mem else ""
            row[f"{stage}_state"] = state
            row[f"{stage}_mem"] = fmt_bytes(mem)
            row[f"{stage}_mem_bytes"] = mem or ""
            row[f"{stage}_mem_source"] = source if mem else ""
            row[f"{stage}_rss_exclusive"] = st.get("rss_is_exclusive", "")
        rows.append(row)
    return rows


def write_report_txt(root: Path, rows: list[dict], path: Path,
                     use_sacct: bool) -> None:
    """The position-wise table, as a plain fixed-width text file."""
    width = max([len(r["position"]) for r in rows] + [len("position")])
    width = min(width, 60)

    ok = {s: sum(1 for r in rows if r[f"{s}_state"] == "ok")
          for s in ("inference", "analysis")}
    lost = {s: sum(1 for r in rows if r[f"{s}_state"] in
                   ("TIMEOUT", "OOM", "CANCELLED", "PREEMPTED", "NODE_FAIL",
                    "ENV_ERROR"))
            for s in ("inference", "analysis")}

    lines = [
        "Position-wise run report",
        f"root:      {root}",
        f"generated: {now()}",
        f"positions: {len(rows)}",
        "",
        f"{'':<{width}}   {'inference':<26}  {'analysis':<26}",
        f"{'position':<{width}}   {'state':<9}{'time':<9}{'memory':<8}  "
        f"{'state':<9}{'time':<9}{'memory':<8}",
        "-" * (width + 3 + 26 + 2 + 26),
    ]
    for r in rows:
        name = r["position"]
        if len(name) > width:
            name = name[:width - 1] + "…"
        lines.append(
            f"{name:<{width}}   "
            f"{r['inference_state']:<9}{r['inference_time']:<9}{r['inference_mem']:<8}  "
            f"{r['analysis_state']:<9}{r['analysis_time']:<9}{r['analysis_mem']:<8}")

    lines += [
        "-" * (width + 3 + 26 + 2 + 26),
        f"{'ok':<{width}}   {ok['inference']:<26}  {ok['analysis']:<26}",
        f"{'lost to slurm':<{width}}   {lost['inference']:<26}  {lost['analysis']:<26}",
        "",
        "states",
        "  ok         finished and recorded a state marker",
        "  adopted     segmentation already on disk with no marker - `status`",
        "             treats these as done rather than redoing them",
        "  failed     ran to a python exception; traceback in pipeline/logs/",
        "  TIMEOUT    killed at the walltime limit - no marker was written, so",
        "             `status` calls this position `pending`",
        "  OOM        killed for exceeding its memory request",
        "  CANCELLED  scancel, or a dependency that could not be satisfied",
        "  running?   log opened, never closed, no marker and no slurm verdict",
        "  pending    no evidence it has run",
        "",
        "time",
        "  measured from the state marker, or for a TIMEOUT from the position",
        "  log's first timestamp against slurm's cancellation time",
        "",
        "memory",
        "  peak resident set size. Recorded in the state marker from the run",
        "  itself; runs from before that was added show '-' unless --sacct is",
        "  given and sacct still has the job. A '-' means not measured, never",
        "  zero. Under slurm one array task is one position, so the figure is",
        "  that position's; a plate run locally in one process reports the",
        "  high-water mark to that point instead.",
    ]
    if not use_sacct:
        lines += ["", "  (--sacct was not used; it can recover memory for older runs)"]

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def cmd_report(args) -> int:
    root = Path(args.root).resolve()
    pdir = pipeline_dir(root)
    if not pdir.is_dir():
        die(f"no {PIPELINE_DIRNAME}/ under {root}; nothing has been run here")

    dirs = job_dirs(root)
    if args.job:
        dirs = [d for d in dirs if d.name == args.job]
        if not dirs:
            die(f"no job directory named {args.job} under {pdir / 'jobs'}")
    elif not args.all_jobs:
        dirs = dirs[-1:]

    log(f"Run report for {root}")
    log(f"Generated {now()}")
    if dirs:
        scope = "all submits" if args.all_jobs else "most recent submit"
        log(f"Reading {len(dirs)} job director{'y' if len(dirs) == 1 else 'ies'} "
            f"({scope}): {', '.join(d.name for d in dirs)}")
    else:
        log("No submit directories - reporting from the per-position logs and "
            "state markers only.")

    records = [r for d in dirs for r in scan_job_dir(root, d)]

    # ---- the headline: work lost to the walltime -------------------------
    # One position retried three times is three timeout records, and scanning
    # every submit finds all of them. Collapse to one row per position, keeping
    # the largest limit that has been tried - that is what the next one has to
    # beat.
    attempts: dict = {}
    for r in (r for r in records if r["kind"] == "timeout"):
        key = (r["stage"], r["position"])
        prev = attempts.get(key)
        if prev is None or ((parse_walltime(r["request"].get("time", "")) or 0)
                            > (parse_walltime(prev["request"].get("time", "")) or 0)):
            attempts[key] = r
    # A position killed on Monday and finished on Tuesday is history, not work
    # to do. Only the ones with no `done` marker still need anything.
    killed = [r for (stage, stem), r in sorted(attempts.items())
              if marker_state(root, stage, stem) != "done"]
    recovered = [r for (stage, stem), r in sorted(attempts.items())
                 if marker_state(root, stage, stem) == "done"]

    log("")
    log("=" * 78)
    if killed:
        log(f"CUT OFF BY THE WALLTIME LIMIT - {len(killed)} position(s) still unfinished")
        log("=" * 78)
        log("These were killed mid-run, so they wrote no state marker and")
        log("`status` reports them as `pending`, not `failed`. Nothing is")
        log("corrupt; the work simply did not finish.")
        if recovered:
            log("")
            log(f"({len(recovered)} more hit the limit on an earlier attempt and have")
            log(" since completed - listed at the end of this section, nothing to do.)")
        log("")
        log(f"  {'stage':<10} {'position':<44} {'largest limit tried':<20} ran for")
        log(f"  {'-' * 10} {'-' * 44} {'-' * 20} {'-' * 9}")
        for r in sorted(killed, key=lambda r: (r["stage"], r["position"])):
            log(f"  {r['stage']:<10} {(r['position'] or '(unmapped)'):<44} "
                f"{r['request'].get('time', '?'):<20} {killed_runtime(r)}")
        if any(killed_runtime(r).endswith("?") for r in killed):
            log("")
            log("  (a `?` in `ran for` means the position log had no opening")
            log("   timestamp to measure against)")
        log("")
        log("  Full logs:")
        for r in sorted(killed, key=lambda r: r["out"])[:args.max_logs]:
            log(f"    {PIPELINE_DIRNAME}/{r['out']}")
        if len(killed) > args.max_logs:
            log(f"    ... and {len(killed) - args.max_logs} more")

        # What to actually run. The limit that was too small is the one to
        # raise, and only the stages that were actually killed need it.
        by_stage = {}
        for r in killed:
            asked = parse_walltime(r["request"].get("time", "")) or 0
            by_stage[r["stage"]] = max(by_stage.get(r["stage"], 0), asked)
        flags = []
        for stage, asked in sorted(by_stage.items()):
            flag = {"inference": "--infer-time", "analysis": "--analysis-time",
                    "flatfield": "--maps-time"}.get(stage)
            if flag:
                doubled = max(asked * 2, 3600)
                h, rem = divmod(int(doubled), 3600)
                flags.append(f"{flag} {h//24}-{h%24:02d}:{rem//60:02d}:00")
        # STUB_RE is the parser the rest of the pipeline uses for this; picking
        # the stub apart with string splits breaks on the first stem whose
        # naming differs.
        found = {m.group(0) for m in
                 (STUB_RE.search(r["position"]) for r in killed) if m}
        stubs = " ".join(sorted(found))
        log("")
        log("  To finish them, raise the limit and resubmit just these:")
        log("")
        log(f"    python pipeline.py submit --root {sh(root)} \\")
        log(f"        --stem {stubs} \\")
        log(f"        {' '.join(flags)} --sbatch")
        log("")
        log("  Without --stem the whole plate is considered, which re-runs")
        log("  anything the platemap has since made stale as well - correct,")
        log("  but a great deal more work than these positions need.")
        if recovered:
            log("")
            log("  Hit the limit earlier but have since completed - no action:")
            for r in recovered:
                log(f"    {r['stage']:<10} {r['position']}")
    else:
        log("CUT OFF BY THE WALLTIME LIMIT - none outstanding")
        log("=" * 78)
        if recovered:
            log(f"{len(recovered)} position(s) hit the limit on an earlier attempt and")
            log("have since completed; nothing is outstanding:")
            for r in recovered:
                log(f"  {r['stage']:<10} {r['position']}")
        elif dirs:
            log("No task in the scanned job directories hit its time limit.")
        else:
            log("There are no submit directories to scan, so nothing here can")
            log("prove a walltime kill either way. A plate run by hand cannot")
            log("be cut off by SLURM, but one submitted from another checkout")
            log("keeps its .out files there, not here.")

    # ---- everything else SLURM ended for us ------------------------------
    other = [r for r in records if r["kind"] and r["kind"] != "timeout"]
    if other:
        groups: dict[str, list[dict]] = {}
        for r in other:
            groups.setdefault(r["kind"], []).append(r)
        log("")
        log("=" * 78)
        log(f"OTHER TASKS SLURM OR PYTHON ENDED - {len(other)} task(s)")
        log("=" * 78)
        for kind, rs in sorted(groups.items()):
            what, remedy = SIGNATURE_HELP.get(kind, (kind, ""))
            log(f"  {kind} - {what} ({remedy}): {len(rs)}")
            for r in sorted(rs, key=lambda r: (r["stage"], r["position"]))[:args.max_logs]:
                log(f"    {r['stage']:<10} {r['position'] or '(unmapped)'}")
                if r["evidence"]:
                    log(f"      {r['evidence'][:150]}")
                log(f"      {PIPELINE_DIRNAME}/{r['out']}")
            if len(rs) > args.max_logs:
                log(f"    ... and {len(rs) - args.max_logs} more")

    # ---- started and never came back, with no SLURM verdict --------------
    # A task whose .out has no signature but whose position log opens and never
    # closes is the same lost work by a different route: the node went away, or
    # the run is still going right now.
    stalled = [r for r in records
               if not r["kind"] and r["position"] and r["log_started"]
               and not r["log_finished"]
               and marker_state(root, r["stage"], r["position"]) == "no marker"]
    if stalled:
        log("")
        log("=" * 78)
        log(f"STARTED BUT NEVER FINISHED - {len(stalled)} task(s)")
        log("=" * 78)
        log("The position log opens and stops, with no marker and no SLURM")
        log("verdict. Either still running, or the job vanished without")
        log("writing one. Check `squeue` before resubmitting.")
        for r in sorted(stalled, key=lambda r: (r["stage"], r["position"])):
            log(f"  {r['stage']:<10} {r['position']:<40} last log line "
                f"{r['log_last'] or '?'}")

    # ---- python-level failures, from the markers -------------------------
    failed = []
    for stage in ("inference", "analysis"):
        d = pdir / "state" / stage
        for path in sorted(d.glob("*.json")) if d.is_dir() else []:
            try:
                st = json.loads(path.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            if st.get("status") == "failed":
                failed.append((stage, st))
    if failed:
        log("")
        log("=" * 78)
        log(f"FAILED WITH AN EXCEPTION - {len(failed)} position(s)")
        log("=" * 78)
        log("These ran to a python error and recorded it, so they are `failed`")
        log("in `status` and will be retried by the next submit.")
        for stage, st in failed:
            last = (st.get("error", "").strip().splitlines() or ["?"])[-1]
            log(f"  {stage:<10} {st.get('stem', '?')}")
            log(f"    {last[:150]}")
            log(f"    {PIPELINE_DIRNAME}/logs/{stage}/{st.get('stem', '?')}.log")

    # ---- what did finish, and how close to the limit ---------------------
    log("")
    log("=" * 78)
    log("COMPLETED, AND HOW MUCH HEADROOM IS LEFT")
    log("=" * 78)
    all_dirs = job_dirs(root)
    for stage, candidates in (("inference", ("infer.sbatch",)),
                              ("analysis", ("analyze.sbatch",
                                            "analyze_ready.sbatch"))):
        runs = stage_durations(root, stage)
        if not runs:
            log(f"  {stage:<10} nothing recorded as done")
            continue
        slowest_pos, slowest = max(runs, key=lambda kv: kv[1])
        median = sorted(v for _, v in runs)[len(runs) // 2]
        line = (f"  {stage:<10} {len(runs)} done, median {fmt_duration(median)}, "
                f"slowest {fmt_duration(slowest)} ({slowest_pos})")
        # The newest submit need not contain this stage's script - a plate that
        # was already segmented gets an analysis array and no infer.sbatch - so
        # walk back until one turns up. Falling through to another stage's
        # script would compare against the wrong limit entirely.
        limit = None
        for d in reversed(all_dirs):
            for name in candidates:
                limit = parse_walltime(sbatch_request(d / name).get("time", ""))
                if limit:
                    break
            if limit:
                break
        if limit:
            used = 100 * slowest / limit
            line += f" - {used:.0f}% of the {fmt_duration(limit)} limit"
        log(line)
        # A completed run that used most of its limit is the next timeout.
        if limit and slowest > 0.8 * limit:
            log(f"             ! the slowest position used {100*slowest/limit:.0f}% "
                f"of its walltime; the next plate will lose positions here")

    log("")
    log("=" * 78)
    log("HOW THE PLATE STANDS NOW")
    log("=" * 78)
    for stage in ("inference", "analysis"):
        d = pdir / "state" / stage
        counts: dict[str, int] = {}
        for path in sorted(d.glob("*.json")) if d.is_dir() else []:
            try:
                st = json.loads(path.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            s = st.get("status", "?")
            counts[s] = counts.get(s, 0) + 1
        log(f"  {stage:<10} {fmt_counts(counts) if counts else 'no markers yet'}")
    if killed:
        log("")
        log(f"  Remember: the {len(killed)} walltime-killed task(s) above are NOT")
        log("  in these counts. They have no marker, so `status` calls them")
        log("  `pending` - the same word it uses for work never submitted.")

    if args.txt is not None:
        out = Path(args.txt) if args.txt != "-" else pdir / "report.txt"
        rows = position_table(root, records, use_sacct=args.sacct)
        write_report_txt(root, rows, out, use_sacct=args.sacct)
        log("")
        log(f"Wrote {out}  ({len(rows)} position(s))")

    if args.csv:
        out = Path(args.csv) if args.csv != "-" else pdir / "report.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        fields = ["job", "jobid", "index", "stage", "position", "kind",
                  "evidence", "out", "log_started", "log_last", "log_finished"]
        with open(out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            for r in records:
                writer.writerow(r)
        log("")
        log(f"Wrote {out}")

    return 1 if killed or other else 0


# ---------------------------------------------------------------------------
# Correction maps
# ---------------------------------------------------------------------------

def manual_map_pairs(root: Path, args) -> list[MapPair]:
    """Correction maps named on the command line rather than in the platemap.

    Blank wells are sometimes acquired in a separate run, so the stacks are not
    always positions of this plate and cannot be described by a well id.
    """
    pairs = []
    for role, values in (("intensity", args.intensity), ("background", args.background)):
        for value in values or []:
            path = Path(value)
            if not path.is_absolute():
                path = root / path
            if not path.exists():
                die(f"{path} does not exist")
            # create_correction_maps takes the channel from the file stem, so
            # the stack must be named ..._<Channel>.tif - check it here rather
            # than let it write a map called e.g. "s8_intensity_map".
            channel = path.stem.split("_")[-1]
            if channel not in VALID_CHANNELS:
                die(f"{path.name}: the name must end in _<channel>.tif with "
                    f"channel one of {', '.join(VALID_CHANNELS)} (got {channel!r})")
            pairs.append(MapPair(channel=channel, role=role, source=path))
    return pairs


def cmd_maps(args) -> int:
    """Build the correction maps the platemap's blank-media wells describe.

    Fluorobrite wells give the background map, DMEM wells the excitation
    intensity map, one of each per channel. Once written they sit in the root
    folder and cellaap_analysis picks them up on its own, so this runs once per
    experiment, before the analysis. With no blank wells declared there is
    nothing to do and the analysis simply runs uncorrected.
    """
    root = Path(args.root).resolve()

    pairs = manual_map_pairs(root, args)
    if pairs:
        log(f"{len(pairs)} correction map(s) named on the command line; "
            f"the platemap's blank wells are ignored for this run")
    else:
        mapfile = Path(args.map).resolve() if args.map else platemap_path(root)
        rows = read_platemap(mapfile)
        positions = discover_positions(root, args.pattern)
        _, map_wells, _ = build_tasks(root, positions, rows)
        if not map_wells:
            log("No background/intensity wells in the platemap, and no stacks "
                "given on the command line. Nothing to do - the analysis will "
                "run without corrections.")
            return 0
        pairs, warnings = resolve_correction_maps(root, map_wells)
        for w in warnings:
            log(f"  ! {w}")

    # Confirm the labels against the wells themselves before anything is built:
    # only a dmem well may produce an intensity map, only a fluorobrite well a
    # background map, and the brightness of the two says which is which.
    pairs, checks = verify_map_roles(pairs)
    for line in report_role_checks(checks):
        log(line)
    check_by_channel = {c.channel: c for c in checks}
    swapped = {c.channel for c in checks if c.verdict == "swapped"}

    # A swapped channel has to rebuild both of its maps: the ones on disk were
    # built from the other well.
    todo = [p for p in pairs if args.force or p.channel in swapped
            or map_status(p) == "pending"]
    for pair in pairs:
        if pair not in todo:
            log(f"  {pair.channel} {pair.role} map already built "
                f"({pair.output.name}); skipping")
    if not todo:
        log("All correction maps are already built.")
        return 0

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import cellaap_analysis  # type: ignore

    tee = tee_to(root, "maps", "corrections")
    saved, sys.stdout = sys.stdout, tee
    started = time.time()
    written: list[MapPair] = []
    missing: list[MapPair] = list(todo)
    crashed = False
    try:
        # What this job is about to cost, in the log, before it costs it. A
        # background map is built by median-filtering every frame of the blank
        # stack, and gen_background_correction_map holds the result as int64 -
        # four times the size of the int16 it eventually writes. That buffer,
        # not the filtering, is what makes this job large.
        for pair in todo:
            gb = pair.source.stat().st_size / 2**30
            note = f", ~{gb * 4:.1f} GB working buffer" if pair.role == "background" else ""
            log(f"  {pair.role:<10} {pair.channel:<10} from {pair.source.name} "
                f"({gb:.2f} GB{note})")
        fs = filesystem_of(root)
        if fs:
            log(f"  root folder is on {fs}")
        for pair in todo:
            for moved in quarantine_stale_maps(root, pair):
                log(f"  ! moved {moved.name} to pipeline/superseded/ - it was "
                    f"built from a different well than this run uses")
        log(f"[{now()}] building {len(todo)} correction map(s)...")

        # One call for the whole batch. create_correction_maps() re-reads every
        # map in the folder at the end of each call, so calling it once per map
        # made a plate with several channels re-read gigabytes it had just
        # written, over and over.
        session = cellaap_analysis.analysis(root, plotting_only=True)
        session.create_correction_maps({p.key: p.source for p in todo})

        written = [p for p in todo if p.output.exists()]
        for pair in written:
            record_map(root, pair, started, check_by_channel.get(pair.channel))
            log(f"  wrote {pair.output.name} "
                f"({pair.output.stat().st_size / 2**30:.2f} GB)")
        missing = [p for p in todo if p not in written]
        for pair in missing:
            log(f"  ! {pair.output.name} was not written")

        elapsed = time.time() - started
        # Read of each source, write of each map, and the one read-back
        # create_correction_maps does at the end. A floor, not the exact
        # figure, but enough to tell compute-bound from I/O-bound.
        moved = (sum(p.source.stat().st_size for p in todo)
                 + 2 * sum(p.output.stat().st_size for p in written))
        log(f"[{now()}] {len(written)}/{len(todo)} map(s) in {elapsed / 60:.1f} min "
            f"(>= {moved / 2**30:.1f} GB moved, "
            f"{moved / 2**20 / max(elapsed, 1):.0f} MB/s effective)")
        if elapsed > 15 * 60:
            # Median filtering a 150-frame 2048x2048 stack is about a minute of
            # CPU. Much more than that is not the filter - it is the machine
            # swapping on the int64 buffer, or the filesystem. `seff` on this
            # job id separates the two: low CPU efficiency means waiting, not
            # computing.
            log(f"  note: that is far longer than the filtering itself costs "
                f"(~1 min per 150-frame 2048x2048 stack). Run "
                f"`seff {os.environ.get('SLURM_JOB_ID', '<jobid>')}` - low CPU "
                f"efficiency means the time went to memory pressure or file "
                f"I/O, not computation. See PIPELINE_README.md.")
    except Exception:
        crashed = True
        log(traceback.format_exc())
    finally:
        sys.stdout = saved
        tee.flush()

    if crashed or missing:
        log(f"{len(missing)} correction map(s) were not built; see "
            f"{pipeline_dir(root) / 'logs' / 'maps' / 'corrections.log'}")
        return 1
    log(f"{len(written)} correction map(s) written to {root}. Analyses started "
        f"from now on will apply them automatically.")
    return 0


# ---------------------------------------------------------------------------
# SLURM submission
# ---------------------------------------------------------------------------

SBATCH_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={account}
#SBATCH --partition={partition}
{extra_directives}#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --time={walltime}
{array_directive}#SBATCH --output={logdir}/{output_pattern}
#SBATCH --mail-user={mail_user}
#SBATCH --mail-type=FAIL,END

# Generated by pipeline.py on {stamp}. The task list and the platemap in this
# directory are frozen copies, so editing the originals does not change a run
# that is already queued.

# No `set -u` here. Conda's activate.d hooks read variables that are not set
# yet - cellaap-env's MKL hook opens with $MKL_INTERFACE_LAYER - and under -u
# that is a fatal error, so the job dies during `conda activate` and never
# reaches python. The submit_*.sh scripts this replaces set no shell options
# at all; the guards below cover what -u was there for.
set -o pipefail
eval "$(conda shell.bash hook)"
conda activate {env} || {{ echo "could not activate the {env} conda environment"; exit 1; }}
cd {module_dir} || exit 1

python pipeline.py {stage} \\
    --root {root} \\
    --map {jobdir}/platemap.csv{stage_args}
"""

# Every stage the template can launch. SBATCH_TEMPLATE passes --root and --map
# to all of them, so a stage listed here that does not accept --map dies on
# argparse the moment SLURM starts it - which is a long way from where the
# mistake was made. build_parser asserts the flag exists on each of these.
SBATCH_STAGES = ("infer", "analyze", "flatfield", "maps")

# One array task = one position.
ARRAY_ARGS = (' \\\n'
              '    --tasks {jobdir}/{tasks} \\\n'
              '    --index "$SLURM_ARRAY_TASK_ID"')


def sh(path) -> str:
    """A path as a shell word. Data folders on turbo/scratch do pick up
    spaces, and an unquoted one would silently truncate the sbatch script."""
    return shlex.quote(str(path))


def submit_one(script: Path, dependency: str = "") -> tuple[str, str]:
    """sbatch one script. Returns (jobid, error); exactly one is non-empty."""
    cmd = ["sbatch", "--parsable"]
    if dependency:
        cmd.append(f"--dependency={dependency}")
    cmd.append(str(script))
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        return "", out.stderr.strip() or f"sbatch exited {out.returncode}"
    return out.stdout.strip().split(";")[0], ""


def cmd_submit(args) -> int:
    root = Path(args.root).resolve()
    module_dir = Path(__file__).resolve().parent
    rows = read_platemap(platemap_path(root))
    positions = discover_positions(root, args.pattern)
    tasks, map_wells, warnings = build_tasks(
        root, positions, rows, correction=correction_params(args))
    if not tasks:
        die("the platemap leaves nothing to run")

    unmapped = [t for t in tasks if not t.mapped]
    if unmapped and not args.allow_unmapped:
        for t in unmapped[:10]:
            log(f"  ! {t.pos.well} ({t.stem}) is not in the platemap")
        die(f"{len(unmapped)} position(s) are not in the platemap. Add them, or "
            f"pass --allow-unmapped to run them as {DEFAULT_CELLTYPE}.")
    for w in dict.fromkeys(warnings):
        log(f"  ! {w}")

    # --stem narrows the plate before anything else looks at it, so the arrays
    # carry only what was asked for. Same matching as `infer`/`analyze`: a full
    # stem, a position stub (B12_s5), or a whole well.
    if args.stem:
        wanted = set(args.stem)
        tasks = [t for t in tasks if t.stem in wanted or t.pos.well in wanted
                 or t.pos.stub in wanted]
        if not tasks:
            die(f"no position matched --stem {' '.join(args.stem)}")
        unmatched = wanted - {t.stem for t in tasks} - {t.pos.well for t in tasks} \
            - {t.pos.stub for t in tasks}
        for u in sorted(unmatched):
            log(f"  ! --stem {u} matched no position")
        log(f"--stem: {len(tasks)} position(s) selected")

    if args.all_positions or args.force_analysis:
        selected = tasks
    else:
        selected = [t for t in tasks
                    if stage_status(root, t, "inference") != "done"
                    or stage_status(root, t, "analysis") != "done"]
    if not selected:
        log("Everything is already done. Nothing submitted.")
        log("To analyze the existing inference folders again, add "
            "--force-analysis.")
        return 0

    # Positions whose inference folder is already there need no GPU at all, and
    # making their analysis wait on a GPU array that would only no-op for them
    # means waiting in the GPU queue for nothing. They go into their own
    # analysis array with no inference dependency, so re-analyzing a folder
    # that was segmented earlier starts as soon as the CPU queue allows.
    if args.analysis_only:
        blocked = [t for t in selected
                   if stage_status(root, t, "inference") != "done"]
        if blocked:
            for t in blocked[:10]:
                log(f"  ! {t.stem}: no inference folder "
                    f"({t.inference_dir.name})")
            die(f"--analysis-only, but {len(blocked)} position(s) have no "
                f"inference output. Drop the flag, or exclude them in the "
                f"platemap.")
        ready, needs_inference = selected, []
    else:
        needs_inference = [t for t in selected
                           if stage_status(root, t, "inference") != "done"]
        ready = [t for t in selected if t not in needs_inference]

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    jobdir = pipeline_dir(root) / "jobs" / stamp
    logdir = jobdir / "slurm"
    logdir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(platemap_path(root), jobdir / "platemap.csv")

    common = dict(
        account=args.account, mail_user=args.mail_user, logdir=logdir,
        module_dir=sh(module_dir), root=sh(root), jobdir=sh(jobdir), stamp=now(),
    )
    name = args.job_name or root.name
    pattern_arg = f" \\\n    --pattern '{args.pattern}'"
    analysis_args = pattern_arg
    if args.semantic_gap:
        analysis_args += f" \\\n    --semantic-gap {args.semantic_gap}"
    if args.force_analysis:
        analysis_args += " \\\n    --force"

    def write_tasks(filename: str, group: list[Task]) -> str:
        (jobdir / filename).write_text("".join(f"{t.stem}\n" for t in group))
        return ARRAY_ARGS.format(jobdir=sh(jobdir), tasks=filename)

    infer_sbatch = analyze_sbatch = ready_sbatch = None

    if needs_inference:
        # These two arrays share one task list, which is what makes
        # `--dependency=aftercorr` correct: analysis element N waits for
        # inference element N and nothing else.
        array = write_tasks("tasks.txt", needs_inference)
        last = len(needs_inference) - 1

        infer_sbatch = jobdir / "infer.sbatch"
        infer_sbatch.write_text(SBATCH_TEMPLATE.format(
            job_name=f"{name}_inf", partition=args.gpu_partition,
            extra_directives=f"#SBATCH --gpus={args.gpus}\n"
                             f"#SBATCH --mem-per-gpu={args.infer_mem}\n",
            cpus=args.infer_cpus, walltime=args.infer_time,
            output_pattern="%x_%A_%a.out",
            array_directive=f"#SBATCH --array=0-{last}%{args.infer_concurrent}\n",
            env=args.infer_env, stage="infer",
            stage_args=array + pattern_arg, **common))

        analyze_sbatch = jobdir / "analyze.sbatch"
        analyze_sbatch.write_text(SBATCH_TEMPLATE.format(
            job_name=f"{name}_ana", partition=args.cpu_partition,
            extra_directives=f"#SBATCH --mem={args.analysis_mem}\n",
            cpus=args.analysis_cpus, walltime=args.analysis_time,
            output_pattern="%x_%A_%a.out",
            array_directive=f"#SBATCH --array=0-{last}%{args.analysis_concurrent}\n",
            env=args.analysis_env, stage="analyze",
            stage_args=array + analysis_args, **common))

    if ready:
        array = write_tasks("tasks_ready.txt", ready)
        ready_sbatch = jobdir / "analyze_ready.sbatch"
        ready_sbatch.write_text(SBATCH_TEMPLATE.format(
            job_name=f"{name}_ana0", partition=args.cpu_partition,
            extra_directives=f"#SBATCH --mem={args.analysis_mem}\n",
            cpus=args.analysis_cpus, walltime=args.analysis_time,
            output_pattern="%x_%A_%a.out",
            array_directive=f"#SBATCH --array=0-{len(ready) - 1}"
                            f"%{args.analysis_concurrent}\n",
            env=args.analysis_env, stage="analyze",
            stage_args=array + analysis_args, **common))

    # The flat field is plate-wide - it is a property of the optics, not of a
    # position - so it is built once, in its own job, before any analysis array
    # is released. Building it inside the array instead would have every task
    # racing to write the same file. The per-position background needs no such
    # job: each analysis task measures its own.
    ff_sources, ff_notes = flatfield_sources(root, map_wells)
    for w in dict.fromkeys(ff_notes):
        log(f"  ! {w}")
    pending_ff = [c for c in ff_sources if flatfield_status(root, c) == "pending"]
    maps_sbatch = None
    if pending_ff:
        maps_sbatch = jobdir / "flatfield.sbatch"
        maps_sbatch.write_text(SBATCH_TEMPLATE.format(
            job_name=f"{name}_flat", partition=args.cpu_partition,
            extra_directives=f"#SBATCH --mem={args.maps_mem}\n",
            cpus=1, walltime=args.maps_time, array_directive="",
            output_pattern="%x_%j.out",
            env=args.analysis_env, stage="flatfield",
            stage_args=pattern_arg, **common))

    log("")
    log(f"{len(selected)} position(s) to run")
    log(f"  job directory : {jobdir}")
    if needs_inference:
        log(f"  inference     : {len(needs_inference)} position(s), "
            f"{args.gpus} GPU, {args.infer_time}, "
            f"{args.infer_concurrent} at a time, env {args.infer_env}")
        log(f"  analysis      : {len(needs_inference)} position(s) after their "
            f"own inference (aftercorr)")
    else:
        log(f"  inference     : nothing to do - every position already has an "
            f"inference folder, so no GPU job is submitted")
    if ready:
        log(f"  analysis      : {len(ready)} position(s) on existing inference "
            f"folders, no wait on the GPU queue")
    log(f"  analysis res. : {args.analysis_cpus} cpu, {args.analysis_mem}, "
        f"{args.analysis_time}, {args.analysis_concurrent} at a time, "
        f"env {args.analysis_env}"
        f"{', forced re-run' if args.force_analysis else ''}")
    if ff_sources:
        done_ff = len(ff_sources) - len(pending_ff)
        log(f"  flat field    : {len(pending_ff)} channel(s) to build"
            f"{f' ({done_ff} already built)' if done_ff else ''}, "
            f"{args.maps_mem}, {args.maps_time}")
    else:
        log(f"  flat field    : no blank pair declared; positions are corrected "
            f"for background only")
    log(f"  background    : measured per position inside each analysis task "
        f"(dilation {args.correction_dilation} px, "
        f"{args.correction_frames} frames)")

    if not args.sbatch:
        log("")
        log("Nothing submitted (add --sbatch to submit). To submit by hand:")
        maps_dep = ""
        if maps_sbatch:
            log(f"  mapid=$(sbatch --parsable {sh(maps_sbatch)})")
            maps_dep = ",afterok:$mapid"
        if infer_sbatch:
            log(f"  jobid=$(sbatch --parsable {sh(infer_sbatch)})")
            log(f"  sbatch --dependency=aftercorr:$jobid{maps_dep} "
                f"{sh(analyze_sbatch)}")
        if ready_sbatch:
            dep = f" --dependency=afterok:$mapid" if maps_sbatch else ""
            log(f"  sbatch{dep} {sh(ready_sbatch)}")
        return 0

    if not shutil.which("sbatch"):
        die("sbatch not found on this machine; run submit from a login node "
            "or use the printed commands there")

    # Flat field first: both analysis arrays depend on it, and a failure here
    # is the one case where submitting nothing else is the right answer.
    maps_job = ""
    if maps_sbatch:
        maps_job, err = submit_one(maps_sbatch)
        if err:
            die(f"sbatch failed for the flat-field job ({err}); "
                f"nothing else was submitted")
        log(f"submitted flat-field job {maps_job}")
    maps_dep = f"afterok:{maps_job}" if maps_job else ""

    jobids = {"flatfield": maps_job, "submitted": now()}
    failed = False

    if ready_sbatch:
        job, err = submit_one(ready_sbatch, maps_dep)
        if err:
            log(f"! the analysis array for existing inference folders was NOT "
                f"submitted: {err}")
            log(f"  sbatch {f'--dependency={maps_dep} ' if maps_dep else ''}"
                f"{sh(ready_sbatch)}")
            failed = True
        else:
            jobids["analysis_ready"] = job
            log(f"submitted analysis array {job} over {len(ready)} existing "
                f"inference folder(s)")

    if infer_sbatch:
        infer_job, err = submit_one(infer_sbatch)
        if err:
            die(f"sbatch failed for the inference array: {err}")
        jobids["inference"] = infer_job
        log(f"submitted inference array {infer_job}")

        dependency = f"aftercorr:{infer_job}" + (f",{maps_dep}" if maps_dep else "")
        job, err = submit_one(analyze_sbatch, dependency)
        if err:
            log(f"! the analysis array was NOT submitted: {err}")
            log(f"  inference {infer_job} is running; submit analysis with")
            log(f"  sbatch --dependency={dependency} {sh(analyze_sbatch)}")
            failed = True
        else:
            jobids["analysis"] = job
            log(f"submitted analysis array {job} (waits on {infer_job})")

    (jobdir / "jobids.json").write_text(json.dumps(jobids, indent=2))
    log("")
    log(f"Watch it with:  python pipeline.py status --root {root}")
    running = [j for k, j in jobids.items() if k != "submitted" and j]
    if running:
        log(f"           or:  squeue -j {','.join(running)}")
    return 1 if failed else 0


# ---------------------------------------------------------------------------
# For downstream notebooks: the platemap as a DataFrame
# ---------------------------------------------------------------------------

def read_platemap_df(path):
    """The platemap as a DataFrame, for reading the plate layout by hand.

    The comment lines at the top of the file need `comment='#'`, which is easy
    to forget; use this instead of a bare pd.read_csv.

    To compile a plate, prefer `cellaap_aggregate.load_experiment()`: it goes
    through `read_platemap`/`build_tasks` rather than the DataFrame, so it
    applies the same position-over-well precedence the run applied.
    """
    import pandas as pd
    df = pd.read_csv(path, comment="#")
    return df.fillna("")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pipeline.py",
        description="Platemap-driven cellaap inference + analysis.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Order of business: init -> edit platemap.csv -> check -> "
               "submit --sbatch -> status")
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp):
        sp.add_argument("--root", required=True,
                        help="folder holding the *phs.tif stacks")
        sp.add_argument("--pattern", default="*phs.tif",
                        help="glob for the phase stacks (default: %(default)s)")
        # These live on every subcommand on purpose. They go into the analysis
        # stage's parameter fingerprint, so if `check` and `analyze` disagreed
        # about their defaults every finished position would report stale.
        sp.add_argument("--correction-dilation", type=int, default=121,
                        metavar="PX",
                        help="how far from a cell the background is measured, "
                             "in fluorescence pixels (default: %(default)s)")
        sp.add_argument("--correction-frames", type=int, default=24,
                        metavar="N",
                        help="frames sampled per position for the background "
                             "(default: %(default)s)")
        sp.add_argument("--correction-block", type=int, default=64,
                        metavar="PX",
                        help="block size for the background grid "
                             "(default: %(default)s)")
        sp.add_argument("--correction-model", default="auto",
                        choices=("auto", "grid", "flatfield"),
                        help="background surface model. 'auto' uses the grid "
                             "and falls back to a two-parameter fit against "
                             "the flat field when a position is too crowded "
                             "to measure a grid (default: %(default)s)")
        return sp

    sp = common(sub.add_parser("init", help="write a platemap template for this folder"))
    sp.add_argument("--force", action="store_true", help="overwrite an existing platemap")
    sp.set_defaults(func=cmd_init)

    sp = common(sub.add_parser("check", help="validate the platemap and show the plan"))
    sp.set_defaults(func=cmd_check)

    sp = common(sub.add_parser(
        "flatfield",
        help="build the plate's flat fields from the blank wells (runs before "
             "the analysis array)"))
    # Every stage SBATCH_TEMPLATE launches is passed --map, so a stage that
    # does not accept it dies on argparse before it runs. This one is
    # submitted by `submit`, so it needs the flag - and it should honour it,
    # since the point of the frozen copy is that a resubmitted job sees the
    # platemap the run started with rather than one edited since.
    sp.add_argument("--map", help="platemap to use (default: <root>/platemap.csv)")
    sp.add_argument("--force", action="store_true",
                    help="rebuild even if already built")
    sp.set_defaults(func=cmd_flatfield)

    sp = common(sub.add_parser("status", help="what has finished, what failed"))
    sp.add_argument("--csv", nargs="?", const="-", default=None,
                    help="also write a per-position CSV (default: pipeline/status.csv)")
    sp.set_defaults(func=cmd_status)

    sp = common(sub.add_parser(
        "report",
        help="human-readable run report from the SLURM and position logs",
        description="Reads the .out files SLURM wrote, joins them back to "
                    "positions through the frozen task lists, and reports what "
                    "each task did - in particular which stacks were killed at "
                    "the walltime limit, which `status` cannot tell you because "
                    "a killed task never writes a state marker."))
    sp.add_argument("--job", metavar="STAMP",
                    help="one submit directory by name (default: the most recent)")
    sp.add_argument("--all-jobs", action="store_true",
                    help="scan every submit directory, not just the last")
    sp.add_argument("--max-logs", type=int, default=12, metavar="N",
                    help="log paths to list per section (default: %(default)s)")
    sp.add_argument("--csv", nargs="?", const="-", default=None,
                    help="also write a per-task CSV (default: pipeline/report.csv)")
    sp.add_argument("--txt", nargs="?", const="-", default=None,
                    metavar="PATH",
                    help="write the position-wise table - state, run time and "
                         "peak memory for each stage - as plain text "
                         "(default: pipeline/report.txt)")
    sp.add_argument("--sacct", action="store_true",
                    help="ask sacct for peak memory of runs whose markers "
                         "predate memory recording (cluster only)")
    sp.set_defaults(func=cmd_report)

    sp = common(sub.add_parser("infer", help="run segmentation (GPU env)"))
    sp.add_argument("--map", help="platemap to use (default: <root>/platemap.csv)")
    sp.add_argument("--tasks", help="frozen task list from submit")
    sp.add_argument("--index", type=int, help="run one position by array index")
    sp.add_argument("--stem", nargs="+", help="run only these positions/wells")
    sp.add_argument("--force", action="store_true", help="redo finished positions")
    sp.set_defaults(func=cmd_infer)

    sp = common(sub.add_parser("analyze", help="track, measure and summarize (img-env)"))
    sp.add_argument("--map", help="platemap to use (default: <root>/platemap.csv)")
    sp.add_argument("--tasks", help="frozen task list from submit")
    sp.add_argument("--index", type=int, help="run one position by array index")
    sp.add_argument("--stem", nargs="+", help="run only these positions/wells")
    sp.add_argument("--force", action="store_true", help="redo finished positions")
    sp.add_argument("--semantic-gap", type=int, default=None,
                    help="override analysis_pars.semantic_gap_closing (frames)")
    sp.set_defaults(func=cmd_analyze)

    sp = common(sub.add_parser(
        "maps", help="build the correction maps the platemap's blank wells describe"))
    sp.add_argument("--map", help="platemap to use (default: <root>/platemap.csv)")
    sp.add_argument("--intensity", action="append", metavar="STACK",
                    help="DMEM stack, named ..._<channel>.tif; repeatable. "
                         "Overrides the platemap, for blanks acquired elsewhere")
    sp.add_argument("--background", action="append", metavar="STACK",
                    help="fluorobrite stack, named ..._<channel>.tif; repeatable")
    sp.add_argument("--force", action="store_true", help="rebuild existing maps")
    sp.set_defaults(func=cmd_maps)

    sp = common(sub.add_parser("submit", help="write and optionally submit the SLURM arrays"))
    d = SLURM_DEFAULTS
    sp.add_argument("--sbatch", action="store_true",
                    help="actually submit; without it the scripts are only written")
    sp.add_argument("--job-name", default=None, help="default: the root folder name")
    sp.add_argument("--account", default=d["account"])
    sp.add_argument("--mail-user", default=d["mail_user"],
                    help="default: $USER@umich.edu (%(default)s here)")
    sp.add_argument("--gpu-partition", default=d["gpu_partition"])
    sp.add_argument("--cpu-partition", default=d["cpu_partition"])
    sp.add_argument("--infer-env", default=d["infer_env"])
    sp.add_argument("--analysis-env", default=d["analysis_env"])
    sp.add_argument("--infer-time", default=d["infer_time"])
    sp.add_argument("--analysis-time", default=d["analysis_time"])
    sp.add_argument("--infer-mem", default=d["infer_mem"], help="per GPU")
    sp.add_argument("--analysis-mem", default=d["analysis_mem"])
    sp.add_argument("--maps-mem", default=d["maps_mem"],
                    help="the background map is a median filter over the whole "
                         "blank stack, so this is larger than it looks")
    sp.add_argument("--maps-time", default=d["maps_time"])
    sp.add_argument("--gpus", type=int, default=1)
    sp.add_argument("--infer-cpus", type=int, default=1)
    sp.add_argument("--analysis-cpus", type=int, default=d["analysis_cpus"])
    sp.add_argument("--infer-concurrent", type=int, default=d["infer_concurrent"])
    sp.add_argument("--analysis-concurrent", type=int, default=d["analysis_concurrent"])
    sp.add_argument("--semantic-gap", type=int, default=None,
                    help="passed through to analyze")
    sp.add_argument("--analysis-only", action="store_true",
                    help="submit no GPU job at all; refuses if any position "
                         "lacks an inference folder")
    sp.add_argument("--force-analysis", action="store_true",
                    help="re-analyze every position, including ones already "
                         "analyzed (for a change the platemap cannot see, such "
                         "as an edit to analysis_pars.py)")
    sp.add_argument("--stem", nargs="+", metavar="POSITION",
                    help="submit only these positions - a full stem, a stub "
                         "like B12_s5, or a whole well. Everything else on the "
                         "plate is left alone.")
    sp.add_argument("--all-positions", action="store_true",
                    help="include positions that are already finished")
    sp.add_argument("--allow-unmapped", action="store_true",
                    help=f"run wells missing from the platemap as {DEFAULT_CELLTYPE}")
    sp.set_defaults(func=cmd_submit)

    # A stage submitted to SLURM but missing a flag the template always sends
    # fails on the compute node, minutes or hours after submit looked fine.
    # Catch it here instead: this runs on every invocation and costs nothing.
    launched = {name: parser for name, parser in
                sub.choices.items() if name in SBATCH_STAGES}
    missing = [name for name, parser in launched.items()
               if not any(a.dest == "map" for a in parser._actions)]
    if missing:
        raise AssertionError(
            f"SBATCH_TEMPLATE passes --map to every stage it launches, but "
            f"{', '.join(sorted(missing))} do(es) not accept it; those jobs "
            f"would die on argparse under SLURM")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    root = Path(args.root)
    if not root.is_dir():
        die(f"{root} is not a directory")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
