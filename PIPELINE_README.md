# `pipeline.py` — running a whole plate

One driver for the whole workflow: fill in a platemap, submit, come back to
finished summaries. It replaces the hand-edited `batch_inference*.py`,
`batch_analysis.py`, `submit_inference.sh` and `submit_analysis.sh`.

The analysis module itself is documented in [README.md](README.md); this file
is only about running it in batch. For why it is built the way it is, see
[PIPELINE_DESIGN.md](PIPELINE_DESIGN.md).

---

## The five minutes that matter

```bash
cd ~/cellapp-analysis
conda activate img-env                # enough for init/check/submit/status

python pipeline.py init   --root /scratch/ajitj_root/ajitj99/ajitj/20251025
$EDITOR /scratch/ajitj_root/ajitj99/ajitj/20251025/platemap.csv
python pipeline.py check  --root /scratch/ajitj_root/ajitj99/ajitj/20251025
python pipeline.py submit --root /scratch/ajitj_root/ajitj99/ajitj/20251025 --sbatch
python pipeline.py status --root /scratch/ajitj_root/ajitj99/ajitj/20251025
```

`--root` is the folder holding the `*phs.tif` stacks. Everything the pipeline
writes goes under it.

---

## 1. `init` — write the platemap

```bash
python pipeline.py init --root <folder>
```

Scans the folder for `*phs.tif`, pulls the well and site out of each file name
(`20251009_HeLa-siRNA_G03_s8_phs.tif` → well `G03`, site 8), and writes
`<folder>/platemap.csv` with one row per well, already filled in as HeLa. The
instructions live in comment lines at the top of that file, so you never have
to come back here to remember the format.

`--force` overwrites an existing platemap. `--pattern` changes the glob if your
phase files are not `*phs.tif`.

## 2. Fill it in

```
well_ids,role,celltype,transfection,drug,channels,model,confluency_est,conf_threshold,skip,notes
B05,,HeLa,pEN2,DMSO,,,,,,
G03,,rpe1,siBUB1,DMSO,Texas Red,,,,,
G03_s9,,rpe1,siBUB1,STLC 5uM,GFP;Texas Red,,,,,this site got drug
H03-H06,,HeLa,none,DMSO,,,,,y,not imaged this run
A01,fluorobrite,,,,Texas Red,,,,,blank - background map
A02,dmem,,,,Texas Red,,,,,blank - intensity map
```

| column | what it does |
| --- | --- |
| `well_ids` | which wells the row describes (see below) |
| `role` | blank for wells with cells. `fluorobrite`/`background` or `dmem`/`intensity` for the blank-media wells (see below) |
| `celltype` | `hela` (default), `u2os`, `rpe1`, `ht1080`. Case-insensitive. Ignored on blank-media rows. |
| `transfection` | free text, e.g. `pEN2`, `siBUB1`, `none` |
| `drug` | free text, e.g. `DMSO`, `STLC 5uM` |
| `channels` | channels to measure, `;`-separated. **Blank = measure every channel found next to the phase file.** Valid: `GFP`, `Texas Red`, `Cy5`. On a blank-media row, the channels that well provides a map for |
| `model` | optional override of the inference model |
| `confluency_est` | optional, `(0, 2000]`, default 1800 |
| `conf_threshold` | optional, `(0, 1)`, default 0.25 |
| `skip` | `y`/`yes`/`true` leaves the well out entirely |
| `notes` | free text, ignored |

**`celltype` is the one that matters.** It sets *both* the inference model and
the `analysis_pars` cell type, so the two stages cannot disagree about what is
in a well:

| celltype | inference model | analysis parameters |
| --- | --- | --- |
| `hela` | `HeLa_focal` | hela (vanilla tracking, 20 px, max cell 4000) |
| `u2os` | `U2OS_focal` | u2os (predictive, 30 px, max cell 9500) |
| `rpe1` | `RPE1_focal` | rpe1 (predictive, 30 px, max cell 9500) |
| `ht1080` | `HT1080_focal` | ht1080 (predictive, 30 px, max cell 9000) |

Put a plain (non-focal) model in the `model` column if you want one.

**`well_ids` accepts** a single well (`G03`, `g3`, `G3` — all normalize to
`G03`), a comma/space separated list (`A01,A02 B01`), a range within one plate
row (`B01-B06`), or one exact position (`G03_s9`). A position-level row beats
the well-level row for that site alone, which is how you handle the one site
that got something different.

**Blank-media wells** are how the fluorescence corrections get built. Mark the
fluorobrite well `background` (or `fluorobrite`) and the DMEM well `intensity`
(or `dmem`), and list the channels they cover. These wells hold no cells, so
they are **never segmented** — they take no GPU time and produce no tracks.
Before the analysis runs, one map per channel per role is built from them and
written into the root folder, where `cellaap_analysis` finds them on its own
and applies them to every measured signal.

Only a `dmem` well can produce an intensity map and only a `fluorobrite` well a
background map — and the pipeline confirms that against the wells themselves
rather than trusting the labels. DMEM autofluoresces and FluoroBrite is
formulated not to, so the DMEM well is always the brighter of the two. If the
labels say otherwise, they are swapped, and `check` and the maps job both say
so in as many words:

```
! role check Texas Red: the well labeled dmem was DIMMER than the well labeled
  fluorobrite. DMEM autofluoresces and FluoroBrite does not, so these labels are
  swapped in the platemap. SWAPPING THEM: intensity map from ..._A02_s1_Texas
  Red.tif, background map from ..._A01_s1_Texas Red.tif (means now 949 vs 119).
  Fix the platemap to make this permanent.
```

The maps are then built the right way round regardless of what the platemap
said, so a mislabeled run still produces correct corrections — but fix the
platemap, because it is what everything else reads. Wells within 5% of each
other are reported as too close to call and left exactly as labeled. The check
is a few sampled frames per stack, not a full read, so it costs nothing.

Declare none and nothing breaks: the analysis runs with no background or
intensity correction, exactly as it does today when no maps are present.
`check` says which of the two corrections you will get.

**Wells you leave out** are run as HeLa, but `check` flags them with a `*` and
`submit` refuses to run until you either add them or pass `--allow-unmapped`.

The `celltype`, `transfection` and `drug` columns are the grouping keys the
downstream compilation uses, so this same file drives that too — you never
transcribe the plate layout twice:

```python
import cellaap_aggregate as agg

whole_df = agg.load_experiment(root, expt_length=150, delta_t=10)
```

`cellaap_aggregate` reads the platemap through this module's own parser, so the
groups it compiles apply the same precedence, `skip` and blank-media rules the
run did. See [README.md](README.md#compiling-a-plate-cellaap_aggregatepy).

If you want the platemap as a DataFrame yourself, use
`pipeline.read_platemap_df()` rather than a bare `pd.read_csv()` — the comment
header needs `comment='#'`.

## 3. `check` — validate before spending GPU hours

```bash
python pipeline.py check --root <folder>
```

```
Positions: 4 to run (1 skipped by the platemap)

position                             well  celltype  model        channels        infer    analysis
---------------------------------------------------------------------------------------------------
20251009_HeLa-siRNA_B05_s1_phs       B05   HeLa*     HeLa_focal   GFP;Texas Red   pending  pending
20251009_HeLa-siRNA_G03_s8_phs       G03   HeLa      HeLa_focal   Texas Red       done     stale
20251009_HeLa-siRNA_G03_s9_phs       G03   HeLa      HeLa         GFP;Texas Red   failed   pending
20251009_HeLa-siRNA_H03_s6_phs       H03   rpe1      RPE1_focal   Texas Red       pending  pending

inference: done 1, pending 2, failed 1
analysis:  pending 3, stale 1
* well not in the platemap; defaulted to HeLa

Correction maps (built once, before the analysis):
  background Texas Red  from 20251009_expt_A01_s1_Texas Red.tif    pending
  intensity  Texas Red  from 20251009_expt_A02_s1_Texas Red.tif    pending
```

It checks every cell type and model name, every channel name, that each phase
stack opens and has more than one page, that every well in the folder is
accounted for, and that the correction maps are unambiguous — two background
maps for one channel, or a stray file with `background`/`intensity` in its
name, both break every analysis in the folder, so they are caught here. It exits non-zero if anything would stop the run. Nothing is
submitted and nothing is written.

The four states:

| state | meaning |
| --- | --- |
| `pending` | not run yet — will run |
| `done` | finished with the parameters the platemap currently asks for |
| `stale` | finished, but under different parameters — will re-run |
| `failed` | ran and raised — will be retried, traceback in the log |

## 4. `submit` — the SLURM jobs

```bash
python pipeline.py submit --root <folder>            # write scripts, submit nothing
python pipeline.py submit --root <folder> --sbatch   # actually submit
```

A GPU array (`cellaap-env`) for inference and a CPU array (`img-env`) for
analysis, one array task per position, the second held on
`--dependency=aftercorr` against the first — so **analysis of a position starts
as soon as that position's inference lands**, not when the whole plate is done.

**Positions that already have an inference folder skip the GPU entirely.**
Their analysis goes into its own array with no inference dependency, so it
starts as soon as the CPU queue allows instead of waiting behind a GPU array
that would do nothing for them. If *every* position already has one — a folder
you segmented last month, or a re-analysis after changing tracking parameters —
no GPU job is submitted at all and `submit` says so. A plate that is part
segmented and part not gets both: the GPU array plus a paired analysis array
for the new positions, and a second analysis array running immediately over the
folders that already exist.

If the platemap declares blank-media wells whose maps are not built yet, a
small CPU job builds them first and every analysis array waits on that
(`afterok`). The maps are plate-wide and every analysis task reads them at
start-up, so they are built once, in one place, rather than by 40 array tasks
racing to write the same files.

Everything is recorded in `<root>/pipeline/jobs/<timestamp>/`: the two sbatch
scripts, the task list, a frozen copy of the platemap, the job ids, and the
SLURM logs. Editing `platemap.csv` afterwards does not change a run that is
already queued.

Without `--sbatch` it writes the scripts and prints the two `sbatch` commands,
which is what you want if you would rather submit them yourself or read them
first.

Positions that are already finished are left out of the arrays; `--all-positions`
includes them.

Common overrides (full list under `submit --help`):

| flag | default | |
| --- | --- | --- |
| `--infer-time` / `--analysis-time` | `0-03:00:00` / `0-06:00:00` | walltime per position |
| `--infer-mem` / `--analysis-mem` | `16g` per GPU / `24g` | |
| `--maps-mem` / `--maps-time` | `48g` / `0-02:00:00` | the background map median-filters the whole blank stack |
| `--infer-concurrent` / `--analysis-concurrent` | 4 / 12 | array tasks running at once |
| `--gpus` / `--analysis-cpus` | 1 / 1 | |
| `--account` / `--mail-user` | `ajitj99` / `$USER@umich.edu` | the address is resolved from your login at submit time |
| `--gpu-partition` / `--cpu-partition` | `gpu` / `standard` | |
| `--infer-env` / `--analysis-env` | `cellaap-env` / `img-env` | conda envs |
| `--job-name` | the root folder name | suffixed `_inf` and `_ana` |
| `--allow-unmapped` | off | run wells missing from the platemap as HeLa |
| `--analysis-only` | off | submit no GPU job at all; refuses if a position has no inference folder |
| `--force-analysis` | off | re-analyze every position, including ones already analyzed |

### Re-analyzing existing inference folders

The common cases need no flags — `submit` works out that the segmentation is
already there and runs the analysis only:

```bash
python pipeline.py submit --root <folder> --sbatch
```

Two flags cover the rest:

* `--analysis-only` refuses to submit a GPU job under any circumstances, and
  fails loudly naming the positions if any of them has no inference folder.
  Use it when you know the segmentation is complete and want that assumption
  checked rather than silently patched with a GPU run.
* `--force-analysis` re-analyzes positions whose analysis is already finished.
  This is what to use after editing `analysis_pars.py` or anything else the
  platemap cannot see — those changes do not mark positions stale, so an
  ordinary submit would report everything done and submit nothing.

```bash
python pipeline.py submit --root <folder> --analysis-only --force-analysis --sbatch
```

## 5. `status` — where things stand

```bash
python pipeline.py status --root <folder>
python pipeline.py status --root <folder> --csv        # also writes pipeline/status.csv
```

Each position's `*_summary.xlsx` also gains a **`corrections` sheet**: for each
role, the blank well the map came from, the map file, the measured DMEM and
FluoroBrite means, whether the role check passed or had to swap the labels, and
whether that correction was actually applied to this position. A position
measured without corrections says so. The summary is the file that outlives the
run, so the provenance of its numbers travels with it.

Counts per stage, and for anything that failed, the last line of its traceback
and where the full log is. The CSV has one row per position with well, cell
type, treatment, model, channels, both stage states, run times and the
inference folder name — it is the thing to hand to a plotting notebook or to
skim when a plate looks odd.

---

## Where everything ends up

```
<root>/
├── platemap.csv                                     <- you edit this
├── 20251009_..._G03_s8_phs.tif                      <- your data
├── 20251009_..._G03_s8_Texas Red.tif
├── 20251009_..._A01_s1_Texas Red_background_map.tif <- from the blank wells
├── 20251009_..._A02_s1_Texas Red_intensity_map.tif  <- named for the well it came from
├── 20251009_..._G03_s8_phs_HeLa_focal_1800_0.25_inference/
│   ├── ..._semantic.tif  ..._instance.tif  ..._scores.tif    <- from infer
│   └── ..._tracks.xlsx   ..._analysis.xlsx  ..._summary.xlsx <- from analyze
└── pipeline/
    ├── state/{inference,analysis}/<position>.json    <- what finished, with which parameters
    ├── state/maps/corrections.json                   <- which map came from which well
    ├── superseded/                                   <- maps built from a well later swapped away
    ├── logs/{inference,analysis,maps}/*.log          <- one log per position
    ├── status.csv                                    <- written by `status --csv`
    └── jobs/20251025_142233/
        ├── platemap.csv  tasks.txt  tasks_ready.txt  <- frozen at submit time
        ├── infer.sbatch  maps.sbatch                  <- only when needed
        ├── analyze.sbatch      (positions being segmented now)
        ├── analyze_ready.sbatch (positions already segmented)
        ├── jobids.json
        └── slurm/                                    <- SLURM stdout/stderr
```

The inference folder name matches what `batch_inference_APJ.py` produced, so
folders from earlier runs are picked up as-is.

---

## Running things by hand

Neither `infer` nor `analyze` needs SLURM. Run them in the right conda env and
they process every pending position in the folder sequentially, or just the
ones you name:

```bash
conda activate cellaap-env
python pipeline.py infer   --root <folder> --stem G03          # one well
python pipeline.py infer   --root <folder> --stem 20251009_..._G03_s8_phs

conda activate img-env
python pipeline.py analyze --root <folder> --stem G03_s8       # one position
python pipeline.py analyze --root <folder> --semantic-gap 5    # override gap closing
```

`--force` redoes positions that are already finished. This is the quickest way
to try tracking parameters on a single movie before committing the plate.

## Correction maps

Normally you do not run this: `submit` builds the maps the platemap's
blank-media wells describe, in their own job, before the analysis array starts.
Run it directly to build them ahead of time, or when working outside SLURM:

```bash
conda activate img-env
python pipeline.py maps --root <folder>
```

It reads the `role` column, builds one map per channel per role, skips maps
that already exist (`--force` rebuilds), and writes them into the root folder.
`cellaap_analysis` picks them up from there on its own — there is nothing to
point at them.

If the blanks were acquired in a **separate experiment**, they have no well id
on this plate; name the stacks directly instead, which overrides the platemap
for that run:

```bash
python pipeline.py maps --root <folder> \
    --background "/other/expt/20251009_blank_G03_s8_Texas Red.tif" \
    --intensity  "/other/expt/20251009_dmem_H03_s6_Texas Red.tif"
```

Both flags are repeatable, one per channel. The stack file name must end in
`_<channel>.tif` — that is where the channel name is read from — and the
pipeline says so rather than silently writing a map called `s8_intensity_map`.

A map is applied to every analysis started *after* it exists. Positions already
analyzed without it keep their uncorrected numbers until re-run with `--force`,
which is why the maps job is ordered ahead of the analysis array.

Nothing here is mandatory. With no blank wells declared and no flags given, the
analysis runs uncorrected — `measure_signal` reports raw means, with the
correction columns left at their neutral values (`_bkg_corr` 0, `_int_corr` 1).

---

## When something goes wrong

**`MKL_INTERFACE_LAYER: unbound variable`, or any other "unbound variable" from
a conda `activate.d` script.** The job died during `conda activate`, before
python ran, so nothing was recorded and every position still reads `pending`.
This was a `set -u` in the generated sbatch scripts, fixed on 2026-08-15.
Job directories written before that fix still contain the broken scripts:
cancel whatever is still queued from them (`scancel <jobid>`, including the
analysis array, which SLURM will otherwise hold forever on a dependency that
can never be satisfied), then re-run `submit` to regenerate and resubmit.

**The correction maps take far longer than a minute or two.** Building them is
about a minute of CPU per 150-frame 2048x2048 blank stack — median-filtering
every frame, plus a Gaussian for the intensity map. If the job runs for an
hour, the time is not going into the filtering. Two things to check:

* `seff <jobid>` on the maps job. Low CPU efficiency (single digits) means the
  job spent its time waiting, not computing — either memory pressure or the
  filesystem. Near 100% means it really was computing, and the stack is much
  larger than assumed. A measured example: 57 min wall, 2:53 CPU (5%), 8 GB of
  48 used — i.e. no memory problem and ~54 min of pure I/O wait, moving about
  6 GB, an effective 2 MB/s. That is a filesystem problem, not a pipeline one.
* Where you ran it. `gen_background_correction_map` holds its result as int64
  while writing int16, so the working buffer is **four times** the size of the
  blank stack: ~5 GB for 150 frames at 2048x2048, ~15 GB for 450. That fits the
  48 GB the maps job requests, but not a login node, where it will swap for
  hours. Run it through `submit`, or with `sbatch`, not interactively.

The maps log (`pipeline/logs/maps/corrections.log`) records the size of each
source stack, the working buffer it implies, which filesystem the root sits on,
and the elapsed time with an effective MB/s, so the next run answers this
directly. If the throughput is a few MB/s, check whether `--root` is on shared
NFS rather than scratch. The reader is not the culprit: `skimage.io.imread`
delegates to tifffile for multipage TIFFs and reads at the same speed.

**A position failed.** `status` prints the last line; the full traceback is in
`pipeline/logs/<stage>/<position>.log` and in
`pipeline/state/<stage>/<position>.json`. Fix the cause and re-run the same
command or resubmit — only the failed positions run again.

**trackpy `SubnetOversizeException`.** Too many candidate links. Lower the
search radius for that cell type in `analysis_pars.py`, or drop the cell type
to a `vanilla` track mode. Then `--force` the affected positions.

**The analysis array died at the walltime limit.** Dense positions can exceed
an hour. Resubmit with `--analysis-time 0-10:00:00`; finished positions are
skipped, so it picks up where it stopped.

**Out of memory in analysis.** The zoomed instance mask is the big allocation.
Resubmit with `--analysis-mem 48g`.

**Everything says `pending` after a run.** The state markers live under
`<root>/pipeline/`; if the root path changed (e.g. `/scratch/...` vs a
symlink), the pipeline is looking at a different folder. Use the same `--root`
you submitted with.

**A well ran as HeLa that should not have.** It was missing from the platemap.
Add it; the affected positions turn `stale` and re-run.

**The corrections were not applied.** The analysis log's first line says
`corrections: background map yes/no, intensity map yes/no`, and `check` lists
the maps it expects. The usual causes are a blank well left unmarked in the
platemap, a `channels` entry on the blank row that does not match the channel
being measured, or positions analyzed before the maps job ran — those keep
their uncorrected numbers until re-run with `--force`.

**`check` reports two maps for one channel.** `cellaap_analysis` loads whichever
one the filesystem walk happens to return last, so this is refused rather than
guessed at. Delete the stale map — usually one left from an earlier blank well
— and re-run.

**You changed `analysis_pars.py` itself.** That is *not* visible to the staleness
check, which only tracks what the platemap says, so an ordinary submit reports
everything done. Re-run the affected positions with `--force` locally, or the
whole plate with `submit --force-analysis` — no GPU job is submitted, since the
inference folders are already there.
