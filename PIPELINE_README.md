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
python pipeline.py report --root /scratch/ajitj_root/ajitj99/ajitj/20251025
```

`status` says what finished. `report` reads the SLURM logs and says what did
not — in particular which stacks were cut off at the walltime limit, which
`status` cannot see.

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

Declare none and nothing breaks. The background correction is measured from
each position's own frames and always applies; the blank pair is only needed
for the **flat field**, so a plate without blanks is corrected for background
alone — which is the larger of the two errors anyway. `check` says which of the
two you will get. See [SIGNAL_CORRECTION_README.md](SIGNAL_CORRECTION_README.md).

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

If the platemap declares a blank pair whose flat field is not built yet, a
small CPU job builds it first and every analysis array waits on that
(`afterok`). The flat field is plate-wide — it is a property of the optics, not
of a position — so it is built once, in one place, rather than by 40 array
tasks racing to write the same file. The per-position background needs no such
job: each analysis task measures its own.

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
| `--infer-time` / `--analysis-time` | `0-00:40:00` / `0-01:00:00` | walltime per position |
| `--infer-mem` / `--analysis-mem` | `12g` per GPU / `25g` | |
| `--maps-mem` / `--maps-time` | `20g` / `0-00:30:00` | the flat-field job; it reads a dozen frames of each blank |
| `--correction-dilation` | `121` | how far from a cell the background is measured, in fluorescence pixels |
| `--correction-frames` | `24` | frames sampled per position for the background |
| `--infer-concurrent` / `--analysis-concurrent` | 4 / 12 | array tasks running at once |
| `--gpus` / `--analysis-cpus` | 1 / 1 | |
| `--account` / `--mail-user` | `ajitj99` / `$USER@umich.edu` | the address is resolved from your login at submit time |
| `--gpu-partition` / `--cpu-partition` | `gpu` / `standard` | |
| `--infer-env` / `--analysis-env` | `cellaap-env` / `img-env` | conda envs |
| `--job-name` | the root folder name | suffixed `_inf` and `_ana` |
| `--stem` | — | submit only these positions: a full stem, a stub like `B12_s5`, or a whole well |
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

### Re-running only what was lost

`report` names the positions SLURM cut off; `--stem` submits exactly those and
leaves the rest of the plate alone:

```bash
python pipeline.py submit --root <folder> \
    --stem B12_s12 B12_s2 B12_s5 \
    --analysis-time 0-04:00:00 --sbatch
```

`report` prints this command already filled in. Without `--stem` the whole plate
is considered, which is correct but re-runs everything the platemap has since
made stale as well — on a 54-position plate that was 51 positions of needless
work to redo 3.

A killed position is picked up by an ordinary submit too, because a walltime
kill writes no state marker and the position reads `pending`. There is one case
where that was not true and now is, worth knowing about because it was silent:

> **Outputs left over from an earlier run no longer count as finished if SLURM
> killed the position after they were written.** A folder with outputs and no
> marker is normally *adopted* as done — that is what makes folders from the old
> `batch_*` scripts usable. But a position killed mid-run also has no marker,
> and the previous run's `*_summary.xlsx` is still sitting there, so the killed
> position looked finished and was skipped by every subsequent submit. Outputs
> older than the kill are now treated as `pending`; outputs newer than it are
> still adopted, since something must have written them afterwards.

## 5. `status` — where things stand

```bash
python pipeline.py status --root <folder>
python pipeline.py status --root <folder> --csv        # also writes pipeline/status.csv
```

Each position's `*_summary.xlsx` also gains a **`corrections` sheet**: the
background level this position measured and how far it drifted, how close to
cells it could work, which blank pair the flat field came from, and whether any
legacy `*_map.tif` files were also on disk. A position measured without a flat
field says so. The summary is the file that outlives the run, so the provenance
of its numbers travels with it.

Counts per stage, and for anything that failed, the last line of its traceback
and where the full log is. The CSV has one row per position with well, cell
type, treatment, model, channels, both stage states, run times and the
inference folder name — it is the thing to hand to a plotting notebook or to
skim when a plate looks odd.

## 6. `report` — what the run actually did

```bash
python pipeline.py report --root <folder>
python pipeline.py report --root <folder> --all-jobs   # every submit, not just the last
python pipeline.py report --root <folder> --csv        # also writes pipeline/report.csv
```

`status` reads state markers. `report` reads the **logs** — the `.out` files
SLURM wrote, joined back to positions through the frozen task lists — and it
exists because of one blind spot the markers cannot cover:

> **A task killed for exceeding its walltime is SIGKILLed, so it never writes a
> marker.** `status` therefore reports it as `pending`: the same word it uses
> for a position that was never submitted. Nothing distinguishes the two except
> the job's `.out` file.

That is the first thing `report` prints, and the reason to run it after every
plate:

```
==============================================================================
CUT OFF BY THE WALLTIME LIMIT - 2 task(s)
==============================================================================
These were killed mid-run, so they wrote no state marker and
`status` reports them as `pending`, not `failed`. Nothing is
corrupt; the work simply did not finish.

  stage      position                          asked for   ran for   state now
  ---------- --------------------------------- ----------- --------- ----------
  inference  20251009_HeLa_G03_s8_phs          0-00:40:00  39m52s    no marker
  inference  20251009_HeLa_G03_s9_phs          0-00:40:00  39m52s    no marker

  To finish them, raise the limit and resubmit. Positions that
  already finished are skipped, so this only re-runs what died:

    python pipeline.py submit --root <folder> --infer-time 0-01:20:00 --sbatch
```

`ran for` is **measured**, not the limit restated: the position log's opening
timestamp against SLURM's `CANCELLED AT`. That distinction matters, because a
task that used 39 of its 40 minutes needs a bigger limit, while one reaped after
four minutes was stuck on something else and a bigger limit changes nothing.

The remaining sections, each omitted when empty:

| section | what it means |
| --- | --- |
| **Other tasks SLURM or python ended** | OOM kills, preemption, node failures, `scancel`, a conda activation that failed before python started, and tracebacks — each with the log line that proves it |
| **Started but never finished** | the position log opens and stops, with no marker and no SLURM verdict: either still running, or the job vanished. Check `squeue` before resubmitting |
| **Failed with an exception** | from the state markers — these *are* `failed` in `status` and the next submit retries them |
| **Completed, and how much headroom is left** | per stage: how many finished, median and slowest runtime, and the slowest as a percentage of the walltime that was requested |
| **How the plate stands now** | marker counts, with a reminder that walltime-killed tasks are not in them |

The headroom section is the one to read when nothing has failed *yet*:

```
  inference  13 done, median 12m09s, slowest 16m43s (..._A08_s1_phs) - 42% of the 40m00s limit
  analysis   28 done, median 41m10s, slowest 57m30s (..._B08_s2_phs) - 96% of the 1h00m limit
             ! the slowest position used 96% of its walltime; the next plate will lose positions here
```

A stage that is already brushing its limit will start losing positions as soon
as a plate gets denser, and this says so before it happens rather than after.

### The position-wise table

```bash
python pipeline.py report --root <folder> --txt          # pipeline/report.txt
python pipeline.py report --root <folder> --txt out.txt  # somewhere else
python pipeline.py report --root <folder> --txt --sacct  # recover memory for older runs
```

One row per position, both stages side by side — state, run time and peak
memory — as a plain text file:

```
                                    inference                   analysis
position                            state    time     memory    state    time     memory
------------------------------------------------------------------------------------------
20250213_HT1080 pPS18_A08_s1_phs    ok       16m43s   4.1 GB    ok       7m51s    9.8 GB
20250213_HT1080 pPS18_A08_s2_phs    ok       11m58s   4.0 GB    TIMEOUT  59m48s   -
20250213_HT1080 pPS18_B08_s1_phs    ok       10m00s   3.9 GB    OOM      -        -
------------------------------------------------------------------------------------------
ok                                  3                          1
lost to slurm                       0                          2
```

The file carries its own legend, so it stands on its own when mailed to someone
or dropped next to the data.

**States** are `ok`, `failed`, `TIMEOUT`, `OOM`, `CANCELLED`, `running?` and
`pending`. A SLURM kill outranks the marker: a position carrying a `done` marker
from an earlier attempt still reads `TIMEOUT` if that is what happened in the
run being reported on.

**Time** comes from the state marker, or — for a `TIMEOUT`, which has no marker
— from the position log's first timestamp against SLURM's cancellation time.

**Memory** is peak resident set size, recorded into the state marker by the run
itself. Three things worth knowing:

* Positions analyzed before this was added show `-`. `--sacct` recovers them
  from SLURM accounting, on a login node, while the job ids are still in the
  accounting window. Without sacct available the flag degrades quietly.
* `-` means *not measured*, never zero.
* Under SLURM one array task is one position, so the figure belongs to that
  position. A plate run locally in one process reports the high-water mark of
  the process to that point instead, which over-reports every position after
  the heaviest; the marker records which case it was.

`report` exits non-zero when it found tasks that SLURM or python ended, the same
convention `check` uses, so `report && echo clean` is a usable gate.

Everything is read-only: `report` opens logs and state files and writes nothing
except the optional CSV.

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
    ├── state/flatfield/{GFP,Texas_Red}.npz           <- the plate's flat fields
    ├── state/surfaces/<stem>_<channel>_bkg.tif       <- background actually subtracted
    ├── superseded/                                   <- maps built from a well later swapped away
    ├── logs/{inference,analysis,flatfield}/*.log     <- one log per position
    ├── status.csv                                    <- written by `status --csv`
    ├── report.csv                                    <- written by `report --csv`
    ├── report.txt                                    <- written by `report --txt`
    └── jobs/20251025_142233/
        ├── platemap.csv  tasks.txt  tasks_ready.txt  <- frozen at submit time
        ├── infer.sbatch  flatfield.sbatch             <- only when needed
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
python pipeline.py analyze --root <folder> --frame-interval 4  # override the metadata
```

`--force` redoes positions that are already finished. This is the quickest way
to try tracking parameters on a single movie before committing the plate.

`--frame-interval` is in **minutes** and is only needed when the acquisition
metadata is wrong or missing. Otherwise the interval is read from the
`Time interval:` line of the `*_metadata.txt` files beside the image stacks,
and `analyze` refuses the position rather than guessing. It matters more than
it looks: the interval converts `min_mitotic_duration` (30 min) into frames, so
it decides what counts as a mitotic episode at all. Both flags pass through
`submit` to the generated cluster command.

## Signal correction

Two corrections are applied, and they are estimated in different places
because they are different kinds of quantity. The full account is in
[SIGNAL_CORRECTION_README.md](SIGNAL_CORRECTION_README.md); the short version:

| | scope | when | needs blanks? |
| --- | --- | --- | --- |
| **background** `D + A(t)F` | per position, per frame | inside each analysis task | no |
| **flat field** `F` | one per plate per channel | the `flatfield` job, before the analysis array | yes |

The background is measured from each position's own cell-free pixels, because
the medium's brightness depends on the well, the position and the time — on the
20250213 plate it runs 115–140 counts across positions and climbs 9–15% through
a movie, which no single blank well can represent. The flat field is a property
of the optics, so one is enough for the plate, and it comes from the
*difference* of the two blank wells: both are `D + A·F` with the same camera
offset, so subtracting one from the other cancels an offset that neither alone
can be separated from.

Normally you do not run this: `submit` builds the flat fields before the
analysis array starts. Run it directly to build them ahead of time, or when
working outside SLURM:

```bash
conda activate img-env
python pipeline.py flatfield --root <folder>
```

It reads the `role` column, builds one flat field per channel from the DMEM
(`intensity`) and FluoroBrite (`background`) wells, skips channels already
built (`--force` rebuilds), and caches them under
`<root>/pipeline/state/flatfield/`. A channel with only one blank declared gets
no flat field and is corrected for background only.

**Check that the channels agree.** The flat field is optics, so every channel
should measure the same vignette; the job prints the centre-to-edge ratio per
channel and warns if they differ by more than 0.05. On the 20250213 plate GFP
gives 1.220 and Texas Red 1.229.

### What lands in the tables

| column | meaning |
| --- | --- |
| `<ch>` | raw mean inside the cell mask |
| `<ch>_corrected` | **the corrected signal** — `(raw − background) / flat field` |
| `<ch>_bkg_corr`, `<ch>_int_corr` | legacy `*_map.tif` corrections, only if such files are still in the folder |

The background each position subtracted is also kept, at grid resolution, in
`pipeline/state/surfaces/` — ~490 kB per position-channel, ~12 MB for a plate.
Open it in Fiji, or read it with `correction_tools.read_background_stack()`.
The `corrections` sheet names the file for each channel. Everything you might
want to do with a correction after the run — inspect it, apply it to numbers
already in a workbook, lend its shape to a crowded position — is in
[CORRECTION_TOOLS_README.md](CORRECTION_TOOLS_README.md).

`<ch>_corrected` appears **per frame** in the `cell_data` sheet of
`*_analysis.xlsx`, so corrected signal *dynamics* can be recovered, and
**per track** in the summary, averaged over each track's mitotic window
alongside `<ch>_corrected_std`.

The correction runs after every channel is measured and **before** the summary
is written, so the summary always reflects it. Since `check` treats a position
as analyzed only when its `*_summary.xlsx` exists, a position whose correction
failed reports `pending` and is re-run, rather than sitting there finished but
uncorrected.

### The old `*_map.tif` maps

The pipeline no longer builds them. `pipeline.py maps` still exists for
reproducing an old run by hand, but nothing submits it, and the aggregate
notebook's finding stands: applying those maps made the field *less* flat, not
more, because the intensity map keeps the camera offset in it. See
[SIGNAL_CORRECTION_DESIGN.md](SIGNAL_CORRECTION_DESIGN.md) §1.

`cellaap_analysis` still loads any such file it finds anywhere under the root,
so a folder from an earlier run will keep populating `<ch>_bkg_corr` and
`<ch>_int_corr`. That is harmless — they are extra columns — but they are a
**different correction** and must not be mixed with `<ch>_corrected`. `check`
warns when it finds them.

A correction applies to every analysis started *after* the flat field exists.
Positions already analyzed keep their numbers until re-run with
`--force-analysis`, which is why the flat-field job is ordered ahead of the
analysis array. Changing `--correction-dilation`, `--correction-frames` or
`--correction-block` marks finished positions `stale`, so `check` and `submit`
will offer to re-run them.

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
  larger than assumed. A measured example: 57 min wall, 2:53 CPU (5%), 8 GB
  used — i.e. no memory problem and ~54 min of pure I/O wait, moving about
  6 GB, an effective 2 MB/s. That is a filesystem problem, not a pipeline one.
* Where you ran it. `gen_background_correction_map` holds its result as int64
  while writing int16, so the working buffer is **four times** the size of the
  blank stack: ~5 GB for 150 frames at 2048x2048, ~15 GB for 450. The default
  `--maps-mem 20g` covers the first comfortably and the second not at all, so
  raise it for long blank stacks. A login node will swap for hours whatever the
  figure says — run it through `submit`, or with `sbatch`, not interactively.

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

**A stage died at the walltime limit.** Run `report` — it names the positions,
says how long each actually ran before it was killed, and prints the resubmit
command with a raised limit:

```bash
python pipeline.py report --root <folder>
```

The defaults are 40 min per inference position and 1 h per analysis position,
which are typical costs; a dense position can take several times that. Finished
positions are skipped on resubmit, so it picks up where it stopped rather than
starting over.

**`status` says `pending` for positions you know were submitted.** That is what
a walltime kill looks like: the task was SIGKILLed before it could write a state
marker, and `status` cannot tell "killed" from "never ran". `report` reads the
SLURM `.out` files and can.

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
