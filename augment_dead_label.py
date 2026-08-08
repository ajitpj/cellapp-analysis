'''
Augment previously generated cellaap *_analysis.xlsx files with the
mitotic/dead label from the hand-labeled classifier
(models/dead_classifier_handlabeled.joblib).

For every *_inference folder under the root folder that contains an
*_analysis.xlsx, the script classifies each detection labeled mitotic by the
semantic segmentation (semantic_smoothed == 1) from its instance-masked phase
crop and rewrites the cell table with updated dead_flag (1 = dead)
and new dead_proba columns. All other sheets are preserved.

Usage:
  python augment_dead_label.py <root_folder>                 # all positions
  python augment_dead_label.py <root_folder> --wells A03 B03 # well filter
  python augment_dead_label.py <root_folder> --suffix _dead  # write copies
  python augment_dead_label.py <root_folder> --phase-offset 20
                                     (default augments the files in place)
'''
import argparse
import re
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import tifffile

from dead_classifier import classify_dead

MODEL_PATH = Path(__file__).parent / "models" / "dead_classifier_handlabeled.joblib"


def augment_file(inference_dir: Path, model, suffix: str = "",
                 phase_offset: int = 0) -> None:
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
    required = {"semantic_smoothed", "frame", "x", "y", "label", "area"}
    if not required.issubset(cell_data.columns):
        print(f"{xlsx.name}: missing columns {required - set(cell_data.columns)}, skipping")
        return

    phase = tifffile.imread(phase_hits[0])
    instance = tifffile.imread(instance_hits[0])

    n_mitotic = int(cell_data.semantic_smoothed.sum())
    print(f"{xlsx.name}: classifying {n_mitotic} mitotic-labeled detections...")
    try:
        label_df = classify_dead(phase, instance, cell_data, model,
                                 phase_offset=phase_offset)
    except ValueError as e:      # frame-count mismatch: never guess an offset
        print(f"{xlsx.name}: SKIPPED - {e}")
        return

    cell_data["dead_flag"] = 0
    cell_data["dead_proba"] = np.nan
    cell_data.loc[label_df.index, "dead_flag"] = label_df.dead_flag
    cell_data.loc[label_df.index, "dead_proba"] = label_df.dead_proba
    sheets[key] = cell_data

    out = xlsx if not suffix else xlsx.with_name(
        xlsx.name.replace("_analysis.xlsx", f"_analysis{suffix}.xlsx"))
    with pd.ExcelWriter(out) as writer:
        for name, df in sheets.items():
            df.to_excel(writer, sheet_name=name)
    print(f"{out.name}: {int(cell_data.dead_flag.sum())}/{n_mitotic} "
          f"detections flagged dead-like")


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
    args = ap.parse_args()

    if not args.root_folder.is_dir():
        raise SystemExit(f"{args.root_folder} is not a directory")

    model = joblib.load(MODEL_PATH)["model"]
    folders = sorted(p for p in args.root_folder.iterdir()
                     if p.is_dir() and "_inference" in p.name)
    if args.wells:
        folders = [f for f in folders
                   if any(f"_{w}_" in f.name for w in args.wells)]
    print(f"{len(folders)} inference folders to process")
    for folder in folders:
        augment_file(folder, model, suffix=args.suffix,
                     phase_offset=args.phase_offset)


if __name__ == "__main__":
    main()
