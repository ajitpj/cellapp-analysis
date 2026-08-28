'''
Augment previously generated cellaap *_analysis.xlsx files with the
mitotic/dead label from the pooled classifier
(models/dead_classifier_pooled.joblib).

For every *_inference folder under the root folder that contains an
*_analysis.xlsx, the script classifies each detection the pipeline would score
- those whose semantic label is mitotic (100 or 101), plus post_peak_frames
after each episode - from its instance-masked phase crop, and rewrites the cell
table with updated dead_flag (1 = dead) and mitotic_proba / dead_proba columns.
Unscored rows keep NaN probabilities. All other sheets are preserved.

By default both files are replaced in place. --suffix writes the augmented pair
alongside the originals; --archive instead renames the originals to *_old.xlsx
and writes the augmented pair under the plain names, so downstream tools that
look for the unsuffixed names pick up the new data with the previous run still
on disk.

The matching *_summary.xlsx is then rebuilt from those labels with
analysis.summarize_data, since a summary is derived entirely from them and a
stale one beside a re-scored analysis file is worse than none.

Usage:
  python augment_dead_label.py <root_folder>                 # all positions
  python augment_dead_label.py <root_folder> --wells A03 B03 # well filter
  python augment_dead_label.py <root_folder> --suffix _dead  # write copies
  python augment_dead_label.py <root_folder> --archive       # keep the old pair
  python augment_dead_label.py <root_folder> --phase-offset 20
  python augment_dead_label.py <root_folder> --post-peak-frames 0  # mitotic only
                                     (default augments the files in place)
'''
import argparse
import re
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import tifffile

from analysis_pars import analysis_pars
from cellaap_analysis import analysis
from dead_classifier import (MITOTIC_SEMANTIC_VALUES, classify_dead,
                             rows_to_classify)

MODEL_PATH = Path(__file__).parent / "models" / "dead_classifier_pooled.joblib"


def _archive_existing(analysis_xlsx: Path) -> None:
    """Move the current analysis/summary pair aside as *_old.xlsx.

    Called just before the augmented pair is written under the plain names, so
    the previous run stays readable next to the new one. Both workbooks have
    already been read into memory by this point, so renaming them is safe.

    An existing *_old.xlsx is never overwritten: it holds the pre-augmentation
    original, and clobbering it on a second run would replace that original
    with an already-augmented copy - losing the only untouched version.
    """
    summary_xlsx = analysis_xlsx.with_name(
        analysis_xlsx.name.replace("_analysis", "_summary"))
    for src in (analysis_xlsx, summary_xlsx):
        if not src.exists():
            continue
        dst = src.with_name(src.name.replace(".xlsx", "_old.xlsx"))
        if dst.exists():
            print(f"  {dst.name} already exists, keeping it and leaving "
                  f"{src.name} to be overwritten")
            continue
        src.rename(dst)
        print(f"  {src.name} -> {dst.name}")


def augment_file(inference_dir: Path, model, suffix: str = "",
                 phase_offset: int = 0, post_peak_frames: int = 0,
                 cell_type: str = "hela", archive: bool = False,
                 frame_interval = None) -> None:
    xlsx = sorted(inference_dir.glob("*_analysis.xlsx"))
    if not xlsx:
        print(f"{inference_dir.name}: no analysis file, skipping")
        return
    xlsx = xlsx[0]

    # phase stack lives in the parent data folder, instance stack alongside.
    # Instance files are named either *_phs_instance.tif or *_instance_movie.tif
    # depending on the pipeline version that produced them.
    stub = re.search(r"[A-H]([1-9]|[0][1-9]|[1][0-2])_s(\d{2}|\d{1})", xlsx.name)
    if stub is None:
        print(f"{xlsx.name}: could not parse well/position stub, skipping")
        return
    stub = stub.group()
    phase_hits = [p for p in inference_dir.parent.glob("*_phs.tif")
                  if f"{stub}_" in p.name]
    instance_hits = [p for p in sorted(inference_dir.glob("*_instance*.tif"))
                     if "scores" not in p.name and "semantic" not in p.name]
    if not phase_hits or not instance_hits:
        print(f"{xlsx.name}: phase or instance stack not found, skipping")
        return

    # the cell table is called 'cell_data' in newer files and 'Sheet1' in older ones
    sheets = pd.read_excel(xlsx, sheet_name=None, index_col=0)
    key = "cell_data" if "cell_data" in sheets else list(sheets)[0]
    cell_data = sheets[key]
    required = {"semantic", "frame", "x", "y", "label", "area"}
    if not required.issubset(cell_data.columns):
        print(f"{xlsx.name}: missing columns {required - set(cell_data.columns)}, skipping")
        return

    # same scored set as the pipeline: mitotic frames plus a short tail after
    # each episode, so summarize_data can see a death on mitotic exit
    n_mitotic = len(rows_to_classify(cell_data, MITOTIC_SEMANTIC_VALUES,
                                     post_peak_frames))
    if n_mitotic == 0:
        print(f"{xlsx.name}: no semantic label in {MITOTIC_SEMANTIC_VALUES}, skipping")
        return
    phase = tifffile.imread(phase_hits[0])
    instance = tifffile.imread(instance_hits[0])

    print(f"{xlsx.name}: classifying {n_mitotic} detections...")
    try:
        label_df = classify_dead(phase, instance, cell_data, model,
                                 phase_offset=phase_offset,
                                 post_peak_frames=post_peak_frames)
    except ValueError as e:      # frame-count mismatch: never guess an offset
        print(f"{xlsx.name}: SKIPPED - {e}")
        return

    # non-mitotic rows keep NaN: the model was never asked about them
    cell_data["mitotic_proba"] = np.nan
    cell_data["dead_proba"] = np.nan
    cell_data["dead_flag"] = 0
    cell_data.loc[label_df.index, "mitotic_proba"] = label_df.mitotic_proba
    cell_data.loc[label_df.index, "dead_proba"] = label_df.dead_proba
    cell_data.loc[label_df.index, "dead_flag"] = label_df.dead_flag
    sheets[key] = cell_data

    if archive:
        _archive_existing(xlsx)

    out = xlsx if not suffix else xlsx.with_name(
        xlsx.name.replace("_analysis.xlsx", f"_analysis{suffix}.xlsx"))
    with pd.ExcelWriter(out) as writer:
        for name, df in sheets.items():
            df.to_excel(writer, sheet_name=name)
    print(f"{out.name}: {int(cell_data.dead_flag.sum())}/{n_mitotic} "
          f"detections flagged dead-like")

    # The summary is derived entirely from the labels just rewritten, so a
    # stale one next to a re-scored analysis file is worse than none. Rebuild
    # it from the file that was actually written.
    summary = analysis.from_analysis_file(
        out, cell_type=cell_type, frame_interval=frame_interval
    ).summarize_data(True, suffix=suffix)
    print(f"{out.name.replace('_analysis', '_summary')}: {len(summary)} tracks "
          f"summarized")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root_folder", type=Path,
                    help="folder containing *_inference directories")
    ap.add_argument("--wells", nargs="+", default=None,
                    help="only process these wells, e.g. --wells A03 B03")
    ap.add_argument("--suffix", default="",
                    help="write augmented copies with this suffix instead of "
                         "updating the analysis files in place")
    ap.add_argument("--phase-offset", type=int, default=0,
                    help="analysis frame f maps to phase page f+OFFSET; use when "
                         "the phase stack keeps leading unsegmented frames "
                         "(e.g. 361 phase frames vs 341 segmented -> 20)")
    ap.add_argument("--post-peak-frames", type=int,
                    default=analysis_pars().post_peak_frames,
                    help="frames scored after each mitotic episode, so a death "
                         "on mitotic exit is seen (default from analysis_pars)")
    ap.add_argument("--archive", action="store_true",
                    help="rename the existing analysis and summary files to "
                         "*_old.xlsx, then write the augmented pair under the "
                         "plain names; cannot be combined with --suffix")
    ap.add_argument("--cell-type", default="hela",
                    help="analysis_pars defaults used for the rebuilt summary "
                         "(default: hela)")
    ap.add_argument("--frame-interval", type=float, default=None,
                    help="minutes between frames; read from the acquisition "
                         "metadata when omitted. Pass it when that metadata is "
                         "wrong - it sets the minimum mitotic duration, so a "
                         "wrong value changes which runs count as a mitosis")
    args = ap.parse_args()

    if not args.root_folder.is_dir():
        raise SystemExit(f"{args.root_folder} is not a directory")
    if args.archive and args.suffix:
        raise SystemExit("--archive and --suffix are mutually exclusive: "
                         "--suffix already leaves the originals untouched")

    # pass the whole bundle: it declares which feature blocks the model takes
    model = joblib.load(MODEL_PATH)
    folders = sorted(p for p in args.root_folder.iterdir()
                     if p.is_dir() and "_inference" in p.name)
    if args.wells:
        folders = [f for f in folders
                   if any(f"_{w}_" in f.name for w in args.wells)]
    print(f"{len(folders)} inference folders to process")
    for folder in folders:
        augment_file(folder, model, suffix=args.suffix,
                     phase_offset=args.phase_offset,
                     post_peak_frames=args.post_peak_frames,
                     cell_type=args.cell_type, archive=args.archive,
                     frame_interval=args.frame_interval)


if __name__ == "__main__":
    main()
