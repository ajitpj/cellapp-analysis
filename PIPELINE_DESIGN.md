# `pipeline.py` — design notes

Why the batch driver is built the way it is. For how to use it, see
[PIPELINE_README.md](PIPELINE_README.md); for the analysis module it drives,
[README.md](README.md).

---

## 1. The problem

The workflow is two programs that cannot run in one process:

| | inference | analysis |
| --- | --- | --- |
| environment | `cellaap-env` (`cell_AAP`, `detectron2`) | `img-env` |
| hardware | 1 GPU | 1 CPU, ~20 GB |
| runtime | 30–40 min per position (~150 frames) | 30–60 min, varies with cell count |
| input | `*phs.tif` | inference folder + the raw stacks |
| output | `*_inference/` with semantic/instance/scores | `*_tracks/_analysis/_summary.xlsx` |

Before this, each was a script edited in place before every run
(`batch_inference_APJ.py`, `batch_analysis.py`) and a matching `submit_*.sh`.
That arrangement has four specific problems, and they set the requirements:

1. **The two scripts each carried their own copy of the truth.** Inference had
   `model_name = 'HT1080_focal'` at the top; analysis had `cell_type = "rpe1"`
   and a glob for `*_HeLa_focal*_inference`. Nothing tied them together, and
   the failure is silent: a plate segmented with one model and analyzed with
   another cell type's tracking parameters produces numbers, just wrong ones.
2. **A plate was one job.** One position that failed 8 hours in — a trackpy
   subnet blow-up, a walltime kill — took the rest with it, and the rerun
   started from the top.
3. **Interrupted work left traps.** `batch_inference_APJ.py` created the output
   directory and *then* wrote three TIFFs into it. A job killed in between left
   a directory that every later run treated as finished (`"{path} exists;
   abort!"`), so the position was silently never segmented.
4. **One cell type per run.** The scripts had a single `cell_type` variable, so
   a plate with HeLa in some wells and RPE1 in others meant two passes with a
   hand-edited glob in between.
5. **The fluorescence corrections were a hand-maintained dict.**
   `batch_analysis.py` carried a `map_dict` of literal file paths and a
   `create_correction_maps` boolean, commented in and out between runs. The
   blank wells are part of the plate layout, so they belong in the same
   document as the rest of it — and, being wells like any other, they were
   also being segmented for no reason.

### Requirements

* One place to say what is in each well; both stages read it.
* HeLa is the default and must need no configuration.
* Mixed cell types on one plate.
* A failure costs one position, not a plate.
* Re-running is always safe and always cheap — the natural way to top up a
  plate as more positions arrive, or to retry failures.
* The blank-media wells that produce the fluorescence corrections are declared
  the same way as everything else, are never segmented, and are optional: with
  none declared the analysis runs uncorrected.
* Nothing destructive, nothing submitted, without the user asking for it.

### Non-requirements

Deliberately out of scope: per-track curation (that is `particle_browser.py`),
cross-well compilation and plotting (`cellaap_aggregate.py`), model training,
and any change to the analysis algorithms themselves. The pipeline is
scheduling and bookkeeping around unchanged code.

`cellaap_aggregate` is out of scope in the sense that the pipeline does not call
it — but it calls the pipeline, importing `read_platemap`, `discover_positions`
and `build_tasks` to resolve the plate exactly as a run resolves it. That makes
this module's platemap parsing a published interface, not a private one; see §6.

---

## 2. Shape of the solution

```
                      platemap.csv        (the user writes this, once)
                           |
        well/position -> role, celltype, treatment, channels, model
                           |
      +--------------------+--------------------+
      |                    |                    |
 role = sample        role = sample     role = background/intensity
      |                    |                    |
 pipeline.py infer    pipeline.py analyze   pipeline.py maps
 (cellaap-env, GPU)   (img-env, CPU)        (img-env, CPU, once per plate)
 inference.configure  cellaap_analysis      cellaap_analysis
   run_inference        .files                .create_correction_maps
      |                  .track_centroids           |
      |                  .measure_signal            |
      v                  .summarize_data            v
 <stem>_..._inference/      |    ^          *_background_map.tif
 semantic/instance/scores --+    +--------- *_intensity_map.tif
      |                     |               (found by _load_maps, applied
      |                     v                to every signal measured after)
      |          tracks/analysis/summary .xlsx
      |                     |
      +---------> pipeline/state/*.json  (what finished, with which parameters)
                  pipeline/logs/*.log    (one per position per stage)
```

`submit` wraps these in SLURM jobs: two arrays plus, when the platemap declares
blank wells whose maps are not yet built, a single job for the maps. `init`,
`check` and `status` never touch the data.

### One file, subcommands, no heavy imports at module scope

The script has to be runnable under *both* conda environments — it is the entry
point of both array jobs — so at module scope it imports only the standard
library. `inference`, `cellaap_analysis`, `numpy`, `tifffile` and `pandas` are
imported inside the subcommands that need them. In particular the platemap is
parsed with `csv`, not `pandas`, because `img-env` guarantees pandas and
`cellaap-env` does not.

The alternative — a shared config module plus two thin scripts — was rejected
because the shared module is the whole thing: parsing, validation, state,
naming and submission are one body of logic, and splitting it across files
would mean three files to keep in step and copy to the cluster instead of one.

---

## 3. The platemap

### Why a document, not flags

The mapping from well to condition is per-experiment data with a natural
tabular shape, it is the thing most likely to be wrong, and it is the thing a
person needs to look at later to know what a number means. That argues for a
file that lives next to the data, is written once, and is read by everything —
rather than command-line flags that vanish into shell history.

### Why the columns are what they are

`celltype`, `transfection`, `drug`, `well_ids` are not chosen freely: they are
what the downstream grouping (then `cellaap_utils.create_wellmap_dict()`, now
`cellaap_aggregate`) already requires. Matching them means one document covers
inference, analysis *and* the compilation of summaries into groups — the user
never transcribes the plate layout twice.

That compilation now goes further and reuses `build_tasks` itself rather than
re-deriving groups from the same columns. Re-deriving was where the two drifted:
a well-level match on the string `_G03_` also matches `_G03_s9_`, so a site the
platemap deliberately moved to another condition was counted under both. The
pipeline resolved that precedence correctly and the compilation did not.
The rest (`channels`, `model`, `confluency_est`, `conf_threshold`, `skip`,
`notes`) are pipeline additions.

### `celltype` as the single knob

The mapping `hela -> (HeLa_focal, hela)` is the design's core move: **one
column feeds both the inference model registry and `analysis_pars`.** They are
resolved from the same row at the same moment, in the same function, so
problem 1 above cannot recur. Everything else is an override on top of it.

Defaults were chosen to be the ones already in use: the `*_focal` models, since
those are what the recent runs used, and HeLa as the fallback because it is the
common case.

### Precedence

Exact position (`G03_s9`) → well (`G03`) → HeLa default. The position level
exists because a single site in a well genuinely does get treated differently
sometimes, and forcing the user to split a well across rows to express that
would push them back to editing scripts.

An unmapped well is *not* an error at parse time — the default has to work —
but `check` marks it with `*` and `submit` refuses to run without
`--allow-unmapped`. The reasoning: defaulting quietly is right for someone
exploring, wrong for someone about to spend 40 GPU-hours.

### `role`, and why the blank wells live in the same document

Background (fluorobrite) and excitation-intensity (DMEM) corrections come from
wells on the same plate, imaged the same way, differing only in having no cells
in them. That makes them plate layout, not configuration: they belong in the
platemap next to everything else, and the previous arrangement — a dict of
literal paths inside `batch_analysis.py` — was the same duplicated-truth
problem as the cell type, one step further downstream.

The column names the *medium*, not the correction. `fluorobrite` and `dmem` are
what is actually pipetted into the well and what the person filling in the map
knows; `background` and `intensity` are accepted too, for whoever thinks of it
by what it produces. Both spellings resolve to the same role.

Three consequences fall out of marking a well this way:

* **It is never segmented.** A blank well has no cells, so inference on it is
  35 minutes of GPU time spent finding nothing. `build_tasks` returns these
  positions on a separate path that neither stage ever sees as work.
* **It is optional.** Declare none and the pipeline says so and carries on;
  `measure_signal` already treats missing maps as no correction, so nothing
  special is needed to support the uncorrected case.
* **One map per channel per role.** A second fluorobrite well for the same
  channel would produce a second map file, and `_load_maps` keeps whichever the
  directory walk yields last. Extra sources are refused at resolution time, and
  duplicates already on disk are a `check` failure rather than a coin toss.

### CSV with comment lines

The template is generated from the folder — one row per well found, HeLa
pre-filled, site counts in `notes` — so the first edit is small. The format
instructions are `#` comment lines at the top of the file itself, where the
user is already looking, rather than in a document they have to find.

The cost is that `pd.read_csv` needs `comment='#'`, which is exactly the kind
of thing that gets forgotten; `pipeline.read_platemap_df()` exists to absorb
it. This was judged a better trade than a sidecar instructions file that drifts
out of date, and much better than YAML/JSON, which do not open in Excel.

---

## 4. Execution model

### The unit of work is a position

Not a plate, not a well. A position is 30–60 minutes of work with one input
file and one output folder, which makes it exactly the right granularity for a
SLURM array task: failures are isolated (requirement 4), progress is
observable, and the resume rule is trivially "skip the ones that are done".

### Two arrays, coupled by `aftercorr`

`submit` writes a GPU array and a CPU array over the *same frozen task list*,
and submits the second with `--dependency=aftercorr:<jobid>`. `aftercorr` pairs
array elements by index: analysis element *N* starts when inference element *N*
succeeds, and never runs if it failed.

That pairing is why the task list is frozen to a file at submit time rather
than recomputed per task. Index *N* must mean the same position in both arrays,
and it must keep meaning that even if someone drops a file into the folder or
edits the platemap while the jobs are queued. For the same reason `submit`
copies the platemap into the job directory and the array tasks read *that* copy.

### Two analysis arrays, because inference is often already done

Re-analysis is not an edge case: tracking parameters get tuned, a plate gets
re-summarized months later, half a plate gets segmented before the rest is
acquired. In all of those the inference folders already exist, and the GPU has
nothing to do.

`submit` therefore splits the selected positions in two by asking a single
question of each — is its inference folder already there?

* **Needs inference** → the paired arrays described above, coupled by
  `aftercorr`.
* **Already segmented** → its own analysis array, `analyze_ready.sbatch`, with
  no inference dependency at all.

Either group can be empty, and the empty one produces no job: a plate that is
fully segmented gets a single CPU array and no GPU allocation is ever
requested; a fresh plate behaves exactly as before.

The first design here kept one array list and let the inference tasks no-op for
positions already done. It is simpler and it is wrong in a way that costs real
time: `aftercorr` holds each analysis element until its inference element has
*run*, so on a busy GPU partition an analysis that needs no GPU still waits in
the GPU queue — potentially hours — to be told there was nothing to do.

Two flags cover what detection cannot infer:

* `--analysis-only` turns the detection into an assertion. It submits no GPU
  job under any circumstances and fails, naming the positions, if any of them
  lacks an inference folder — for when the user believes segmentation is
  complete and wants that belief checked rather than quietly patched with a
  40-minute GPU run.
* `--force-analysis` re-analyzes positions whose analysis is already finished.
  It exists because of a known limit of the staleness model (§8): edits to
  `analysis_pars.py` are invisible to it, so without this flag an ordinary
  submit reports everything done and submits nothing.

### Trusting the wells over the labels

`role` says which medium is in a well, and a mislabeled pair would invert both
corrections — the background map built from autofluorescent DMEM, the intensity
map from a well with almost no signal — while everything downstream looked
normal. Nothing in the file names or the metadata can catch that, but the data
can: FluoroBrite is formulated for low background and DMEM is not, so the DMEM
well is severalfold brighter, every time.

So before any map is built, a few frames of each blank stack are sampled and
the two means compared. Brighter-is-DMEM decides which well feeds which map,
and the platemap labels are treated as a claim to be checked rather than an
instruction to be followed.

Three outcomes, and the reasoning behind each:

* **Consistent** — build as labeled.
* **Inverted** — swap the two sources and say so loudly. Swapping rather than
  failing, because a mislabeled pair still contains both media; only the names
  are wrong, and refusing to proceed would cost a queue cycle to fix something
  the data has already resolved. The alert is deliberately unmissable, since
  the platemap remains wrong for everything else that reads it.
* **Within 5%** — too close to call, so nothing is swapped. A bare `>` test
  would flip the labels on noise if someone ever imaged two similar blanks,
  which is a worse failure than the one being guarded against: it would be
  silent and it would be wrong half the time.

Sampling a few frames rather than reading the stacks is what makes this free,
and it is enough — the media differ by a factor of several, not a few percent.

A swap has a consequence on disk. Map files are named after their source stack,
so the maps from a previous, wrongly-labeled run sit alongside the new correct
ones, and `_load_maps` would choose between them by directory order. Those
files are moved to `pipeline/superseded/` — moved rather than deleted, since
they are data, and the move is logged.

The verdict is recorded three times over, because the point is that a person
later can tell what happened: in the maps log, in the ledger, and in a
`corrections` sheet appended to every summary workbook. The last one matters
most — the summary is what gets opened months later, and it now carries which
well each correction came from and whether the labels had to be corrected.

### The maps job

Cost, measured rather than assumed (2048x2048 frames): a 3x3 median filter is
~0.29 s per frame, so a 150-frame blank stack is ~45 s of filtering; the
intensity map's slab mean and sigma=45 Gaussian together are under a second.
Roughly a minute of CPU for a plate.

What is *not* small is the memory. `gen_background_correction_map` allocates
`np.zeros_like(stack, dtype=int)` — int64, four times the int16 it eventually
writes: 4.7 GB for 150 frames, 14 GB for 450. That buffer, not the arithmetic,
is why this job asks for 48 GB and why it must not be run on a login node. The
upstream one-word fix (`dtype=np.int16`, matching what the function already
writes) would cut it by 4x; it lives in `cellaap_utils`, so it is recommended
rather than made here.

The other cost is self-inflicted and was removed: `create_correction_maps()`
ends by calling `_load_maps()`, which re-reads every map file in the folder.
Calling it once per map — the first version here — made a four-channel plate
re-read gigabytes it had just written, once per channel. The batch is now a
single call, and the session is built with `plotting_only=True` so the initial
scan does not read them a further time.

The correction maps are plate-wide, every analysis task reads them at start-up,
and building one is a median filter over an entire blank stack — a few GB of
working set and minutes of CPU. Building them inside the analysis array would
have every task computing the same maps and racing to write the same files, so
they get their own single job, submitted ahead of the analysis array, which
then carries a second dependency: `aftercorr:<inference>,afterok:<maps>`.

It is submitted only when the platemap declares blank wells whose maps are not
already built, so the common case of a second run over the same plate is still
two jobs.

Ordering matters here in a way that is easy to miss: a correction map affects
analyses started *after* it exists. Positions analyzed before the maps landed
keep their uncorrected numbers, which is why the dependency is on the job
rather than a best-effort "build them if you can".

Alternatives considered:

* **`afterok` on the whole array** — analysis waits for the slowest position.
  On a 40-position plate that idles the CPU stage for hours.
* **One job per position running both stages** — needs both conda environments
  and a GPU allocation held through the 45-minute CPU analysis. Wasteful, and
  the two environments are separately installed for a reason.
* **Building the maps inside `submit`** — simplest to write, but it turns a
  bookkeeping command that runs on a login node into one that allocates several
  GB and runs a median filter there.
* **A workflow engine (Snakemake/Nextflow)** — correct in the abstract, but it
  is a dependency, a DSL and an install for the user to maintain, for a DAG
  with two nodes. The state directory here does the same job in code the user
  can read.

The trade-off accepted: analysis elements sit `PENDING` in the queue for hours
holding queue slots (not resources) while inference runs.

### Concurrency and resources

Defaults come from the scripts this replaces (`--account ajitj99`, gpu/standard
partitions, 12 GB per GPU rounded to 16, 20 GB analysis memory rounded to 24;
the notification address is `$USER@umich.edu`, resolved at submit time because
SBATCH directives are not shell-expanded)
with array throttles — 4 concurrent GPU tasks, 12 CPU — chosen to be polite on
a shared allocation rather than optimal. All are flags; none require editing
the script, which was the point.

---

## 5. State, idempotency and resume

Every position × stage has a JSON marker under `pipeline/state/<stage>/`
recording status, timing, host, SLURM job id, outputs, any traceback, and —
importantly — the **parameters it ran with**. Status is then computed by
comparing three things: the marker, the parameters the platemap asks for now,
and whether the outputs are actually on disk.

| marker | outputs | current params | → |
| --- | --- | --- | --- |
| none | absent | — | `pending` |
| none | present | — | `done` (adopted) |
| done | present | match | `done` |
| done | present | differ | `stale` → re-runs |
| done | absent | — | `pending` (outputs were removed) |
| failed | — | — | `failed` → retried |

Three deliberate choices here:

**Adoption.** A finished inference folder with no marker is treated as done,
not redone. Folders produced by the old `batch_inference*.py` scripts are
therefore usable immediately, which is what makes this adoptable mid-experiment
rather than something that starts by recomputing a month of GPU time.

**Parameters, not content hashes.** Staleness compares the small set of values
that actually change the output — model, confluency, threshold for inference;
cell type, channels and source folder for analysis. Hashing the input stacks
would be more thorough and would cost minutes per position for a case (someone
silently replacing a raw stack) that does not happen. The honest limit is
recorded in §8.

**The outputs-absent row.** The marker alone is not trusted, because the
cheapest way for a user to force a redo is to delete a folder, and that should
just work.

### The correction-map ledger, and a name that had to be avoided

Maps are recorded differently from positions, because the map file *is* the
record: its name is derived from the source stack, so asking for a map from a
different blank well asks for a different file, and "already built" is just a
file existence test. The JSON ledger alongside it is provenance — which well,
which channel, when — not control flow.

That ledger is a single file called `corrections.json`, and the name is not
cosmetic. `cellaap_analysis._load_maps()` walks the *entire* root folder at the
start of every analysis, treats any file whose name contains `background` or
`intensity` as a correction map, and reads it with `imread`. A per-role state
file named `Texas Red_background.json` would therefore be loaded as an image
and crash every analysis on the plate. The same hazard is what
`check_map_files()` looks for on the user's behalf: a stray file matching those
words with no channel name in it makes `_load_maps` call `.group()` on a failed
match and die before any analysis starts.

### Atomic inference writes

Inference writes into `<dir>.partial` and renames only after every stack is on
disk *and* the frame count of the output matches the page count of the phase
stack. A job killed at the walltime limit therefore leaves either nothing or a
`.partial` the next run clears — never something that looks finished. This is
the direct fix for problem 3, and the frame-count check catches the subtler
version where inference returns a short movie.

`shutil.rmtree` is called in exactly two places, both on `.partial`
directories the pipeline itself created. Real inference folders are never
deleted; a parameter change writes to a *new* folder name and leaves the old
one alone.

### Logging

Each position gets `pipeline/logs/<stage>/<position>.log`, appended to, in
addition to the SLURM `.out` file. Logs keyed by position rather than by array
task id are what make a failure two weeks later findable: you know the well,
not the job number.

---

## 6. Contracts with the rest of the module

The pipeline calls unchanged code and depends on these facts about it. Each is
a place where a change elsewhere breaks the pipeline quietly, so each is worth
knowing:

| depends on | where | if it changes |
| --- | --- | --- |
| model names | `inference.get_model()` registry, duplicated in `VALID_MODELS` | `check` rejects a valid model, or passes an invalid one through to the GPU node. Update both. |
| cell types | `analysis_pars.cell_types`, mirrored in `CELLTYPES` | same; `CELLTYPES` must stay a subset |
| channel names | `cellaap_analysis.files()` regex `GFP|Texas Red|Cy5`, mirrored in `VALID_CHANNELS` | channel auto-detection and validation drift |
| position stub | `[A-H]nn_s<site>`, from `cellaap_analysis.files()` | positions stop being discovered |
| inference folder name | `<stem>_<model>_<confluency>_<threshold>_inference`, from `batch_inference_APJ.py` | old folders stop being adopted; `analysis.files()` finds channel stacks via the stub, so the folder must stay a sibling of the raw stacks |
| output file names | `*_semantic.tif`, `*_instance.tif` (`files()` globs for the words) | completion detection breaks |
| analysis API | `analysis(root, plotting_only)`, `.files()`, `.track_centroids()`, `.measure_signal()`, `.summarize_data()` | the analyze stage breaks loudly, which is fine |
| map discovery | `_load_maps()` walks the whole root, matching `background\|intensity` and a channel name anywhere in the file name | the ledger name, the duplicate-map check and the stray-file check in `check` all exist because of this rule |
| map naming | `create_correction_maps()` derives both the channel and the output name from the source stem (`..._<channel>.tif`) | `MapPair.output` would predict the wrong file, and "already built" would always be false |
| tifffile API | `create_correction_maps()` used `tifffile.imsave`, an alias removed from recent releases | it raised `AttributeError` on any current install; fixed in `cellaap_analysis` (see below) rather than worked around here |

The `tifffile.imsave` case is worth recording, because it was a real break
rather than a theoretical one. `create_correction_maps()` called `imsave`, an
alias tifffile removed; on tifffile 2026.8 the call raises
`AttributeError: module 'tifffile' has no attribute 'imsave'`, which was
reproduced against the real module. Correction maps could not be built at all
on a current install.

This driver carried a shim for a while (`imsave = imwrite`, which is exactly
what the alias was) on the principle that the pipeline is a caller of
`cellaap_analysis` and not a patch on it. That was the wrong place for it once
the break was confirmed: the shim hid a bug from anyone using the module
directly. Both call sites in `create_correction_maps` now use `imwrite`, and
the shim is gone.

The reads in the same function moved from `skimage.io.imread` to
`tifffile.imread` at the same time. Measured on scikit-image 0.26, the two are
equivalent for multipage TIFFs — `skimage.io.imread` delegates to tifffile,
same data, same time to within noise — so this is a consistency and
version-robustness change, not a speed-up: the function already read intensity
maps through tifffile, and skimage's plugin dispatch has varied across
releases. The stack reads in `files()` and `measure_signal()` still go through
skimage and were left alone; they are the analysis hot path, not map
generation.

`VALID_MODELS` and `VALID_CHANNELS` are duplicated rather than imported on
purpose: `check` has to validate a platemap in `img-env`, where importing
`inference` (and therefore detectron2) is impossible.

### What depends on *this* module

One thing does, and it is the reverse of every row above.
`cellaap_aggregate.py` imports `read_platemap`, `discover_positions`,
`build_tasks`, `normalize_well`, `expand_well_token`, `platemap_path`,
`STUB_RE`, `ROLES` and `MAP_ROLES` from here, so that a compiled group is
resolved by the same code that resolved the run. The alternative — re-deriving
groups from the platemap's columns — is what produced the double-counting bug
described in §3.

Two consequences. First, these names are now interface: renaming `build_tasks`
or changing what `Task` carries breaks compilation, loudly. Second, the
no-heavy-imports-at-module-scope rule (§2) is what makes this import cheap
enough to sit at the top of a notebook module; keep it.

The analysis stage also reuses one `analysis` session across positions when run
sequentially, deleting `summaryDF` and `tracked` between them — the same guard
`batch_analysis.py` used, for the same reason.

---

## 7. Safety

* `submit` writes scripts and prints commands; it submits only with `--sbatch`.
  Spending a shared allocation is the user's decision, made explicitly.
* If the inference array submits but the analysis array does not, the script
  says so, prints the exact `sbatch --dependency=aftercorr:<id>` command to
  recover, and exits non-zero — it does not leave the user thinking both went.
* `init` refuses to overwrite an existing platemap without `--force`. That file
  is hand-written data.
* Paths are shell-quoted in generated sbatch scripts.
* Validation is front-loaded into `check` (readable stacks, real models, real
  channels, mapped wells) so that the expensive failure mode — finding out on
  the GPU node, 40 minutes in — mostly cannot happen.

---

## 8. Known limits

Things the design does not do, stated so they are not discovered the hard way:

* **Editing `analysis_pars.py` does not mark anything stale.** Staleness only
  tracks what the platemap says. Change tracking parameters and you must re-run
  the affected positions with `--force`, or the plate with
  `submit --force-analysis`.
* **`--semantic-gap` is likewise not part of the staleness key.** Changing it
  needs `--force`.
* **State is keyed by file stem.** Renaming a raw stack orphans its markers and
  the position looks unstarted.
* **Changing `conf_threshold` shows as `pending`, not `stale`,** because the
  threshold is part of the inference folder name: the new parameters point at a
  folder that does not exist yet. The old folder stays on disk, deliberately.
* **A correction map only affects analyses started after it exists.** Positions
  analyzed before the maps job ran keep their uncorrected numbers; they are not
  automatically marked stale, because the maps are not part of the analysis
  staleness key. Re-run them with `--force`. The job ordering makes this rare
  rather than impossible.
* **One map per channel per role, plate-wide.** There is no notion of a
  correction that applies to half a plate.
* **No per-position resource tuning.** One walltime and memory figure for the
  whole analysis array, sized for the worst position.
* **Inference runs one position per GPU task**, reloading the model each time
  (~seconds against ~35 minutes, so not worth batching, but it is a choice).
* **`check` does not verify the model downloads.** The first GPU task fetches
  from Zenodo via pooch; no network on the compute node fails there, not in
  `check`.
* **The dead-label augmentation (`augment_dead_label.py`) is not wired in.** It
  remains a separate manual step.

---

## 9. How this was verified

The two stages call code that needs a GPU and a 4 GB conda environment, so
end-to-end verification on this workstation was not possible. What was done:

* A synthetic folder (4 positions, 3 wells, two channels) exercised `init`,
  `check`, `status`, `status --csv` and `submit` for real, including well
  ranges, `skip`, position-level overrides, unmapped-well refusal, warning
  aggregation, and shell quoting of a root path containing a space.
* `infer`, `analyze` and `maps` were driven end-to-end against stub
  `inference`, `cellaap_analysis`, `numpy` and `tifffile` modules, covering: a
  fresh run, the skip-when-done path, stale detection and re-run after a
  platemap edit, the frame-count mismatch path (partial directory discarded,
  failure marker written, no leftover folder), and `--index` / `--stem`
  selection.
* For the re-analysis split: all four combinations — every position already
  segmented (no GPU job written or submitted), a mixed plate (paired arrays
  plus a second immediate analysis array over the existing folders, with
  disjoint task lists), `--analysis-only` refusing a plate with a missing
  inference folder, and `--force-analysis` carrying `--force` into the array
  task, which then re-ran a position whose analysis was already recorded done.
* For the correction maps specifically: blank wells excluded from both stages,
  maps built once and skipped on re-run, the ledger written, the analysis
  session reporting both maps present afterwards, the no-blank-wells path
  reporting "uncorrected" and submitting no maps job, and both `check`
  failures — two maps for one channel, and a stray `background`-named file.
* The real `inference` and `cellaap_analysis` calls are therefore exercised
  only through their signatures. The first cluster run is the real test of
  those two call sites.
