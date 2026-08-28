# cellapp-analysis module

A python module to track cells and measure fluorescence and mitotic duration using cellapp inference files and raw fluorescence images. It runs trackpy on the instance segmentation file first, and then measures cell state (mitotic/non-mitotic) from the semantic segmentation and raw mean fluorescence value from the specific fluorescence channel. The module also summarizes the analysis by calculating the mean fluorescence for each cell over the duration of mitosis, the total mitotic duration, and correction factors based on intensity and background correction maps (must be acquired using empty wells with DMEM and fluorobrite respectively).
This analysis and summary are saved as separate excel spreadsheets.

## Installation

```bash
conda env create -f environment.yml
conda activate img-env
```

This covers the analysis module and the napari curation browser. The upstream
cellaap inference itself needs `cell_AAP` and `detectron2` and is installed
separately.

Batch runs across a whole plate (platemap, SLURM array jobs, resume after
failures) are handled by `pipeline.py`, documented separately in
[PIPELINE_README.md](PIPELINE_README.md) and [PIPELINE_DESIGN.md](PIPELINE_DESIGN.md).

## Usage

1. **Specify the root folder**. This folder must contain cellapp-generated inference folders and the raw intensity stacks and is input as a Path object.

```python
experiment1 = cellapp_analysis.analysis(Path(root_folder), plotting_only: False)
```

If the boolean input is set False, the experiment1 object will look for and read in correction maps (if they are present; not used otherwise). Otherwise, the object waits for the path to a root_folder containing cellapp inference folders.

2.**Measurement mode**: Set the plotting_only mode to False. With this, analysis object instatiation will detect image stacks with either "channel_background" and "channel_intensity" in the filename.

```python
exp_analysis = cellaap_analysis.analysis(Path(root_folder), False)
```

If present, these stacks will be loaded. The analysis object can also create the maps by reading the images stacks from a user-provided dictionary. The dictionary must have "channel_intensity" and "channel_background" as the keys and Path objects pointing to the corresponding image stacks (DMEM and Fluorobrite respectively) that are present in the root folder. Multiple channels can be specified as long as the "intensity" and "background" keywords are in the dictionary keys. e.g.,

```python
map_dict = {'GFP_background': Path(GFP_fluorobrite_stk),
            'GFP_intensity' : Path(GFP_DMEM_stk), 
            'TRed_background': Path(TRed_fluorobrite_stk),
            'TRed_intensity':Path(TRed_DMEM_stk)}

exp_analysis.create_correction_maps(map_dict)
```

These correction maps are optional; analysis will progress without them.

**A note about default analysis parameters:**
These are set by the **analysis_pars** class. If necessary, you can change them to achieve satisfactory results. The object is stored as self.defaults and can be reset as such. Some of the tracking parameters can be changed on the fly; refer to the self.track_centroids function for details.

The main parameters relate to trackpy configuration:

1. max_cell_size = 9500 (pixels) - larger segments are filtered out. This needs to be changed for large cells or cells that tend to spread out.
2. max_pixel_movement (pixels) - the search radius for trackpy. Depends on whether or not the cells crawl. It's set at 20 pixels for Hela, 22 for U2OS/RPE1/HT1080.
3. track_mode - "vanilla" (for cells that don't move much at all, e.g. HeLa)
                "predictive" (for cells that move; all the rest)
                "adaptive" (could be used for cells that move)
4. min_track_length = 10: Only cells tracked for > 10 timepoints are analyzed.
5. memory = 1: tracking memory in timepoints
6. min_mitotic_duration = 30 **minutes** (unrelated to trackpy); mitotic
   events shorter than this are filtered out. It is held in minutes and
   converted to frames against the acquisition's own interval, so the number of
   frames it becomes depends on the experiment: 7 frames at 4 min/frame,
   3 at 10.
7. frame_interval — **no default**. Read from the microscope's metadata
   (`*_metadata.txt`, the `Time interval:` line) by `files()` and
   `from_analysis_file()`, validated, and recorded in the output's
   `parameters` sheet. Pass `frame_interval=` to either to override it when
   the metadata is wrong. It used to be hard-coded to 10, which silently made
   `min_mitotic_duration` 12 minutes rather than 30 on any 4 min/frame
   acquisition.
8. semantic_gap_closing = 3 **frames**, and deliberately *not* converted by
   the frame interval. It closes single-frame flicker in the segmentation's
   output, and a classifier error is one frame wide whether frames are 4 or 10
   minutes apart.

Be careful when using the "predictive" tracking mode. It's very powerful, but can be computationally costly if the tracking memory>1 and max. pixel movement is > 25. This will lead to trackpy exceeding the max. number of nodes in one or more subnetworks. Currenlty, trackpy just exits on this error, which can be problematic when analysis is being done in batch mode. If you come across this issue, gradually decrease the max. pixel movement parameter to get under this error.

**Step 1:** Point the analysis object to a specific inference folder by providing a path to it. This will read the instance segmentation file from this folder.

**Step 2:** Use the **track_centroids** method; it will erode the instance segmentation with the default footprint (needs to be customized for different cells), track the resultant masks using trackpy, and then determine the cell-state by reading the semantic stack. The intermediate padas dataframe can be saved to excel using the flag. See the notes above regarding optimal tracking parameters. 

*It seems that trackpy is using 32-bit integers for assigning labels to individual points. Overrun of this number leads to missing signal measurements, which shouldn't be a big issue. But this needs to be adddressed at some point.*

**Step 3:** Use the **measure_signal** function to measure the fluorescence from the specified channel. The channel string must match the channel name in the file names. The "id = -1" will make the function measure data for all cells that went through a complete mitosis during the time lapse. Optionally, one can provide a list with cell numbers (development only). Thus, cells that remained in interphase throughout the experiment are not measured. Their tracks are still reported.

**Step 4:** Use the **summarize_data** function to create the summary Excel file that lists the average signals measured for all channels, duraion of mitosis, and the correction factors to account for background and excitation intensity variation. Before computing the summary measurements, **gaps in the semantic label vector are filled by "closing" with a footprint (semantic_footprint) of width semantic_gap_closing = 3 frames, so only short gaps are filled. The median filter and the closing are applied per particle, over that track's own frames in frame order** — run over the whole table they would bleed across track boundaries, letting the end of one cell's trace close a gap at the start of the next.
### How a cell is summarized

Two rules, applied per track.

**Mitotic episodes** are the runs of mitotic semantic label at least
`min_mitotic_duration_in_frames` long, read straight off the trace. The first
episode is the mitosis that gets reported; `n_peaks` says how many there were.
A later episode is usually the same cell re-rounding, most often to die.

**Death** is called when `P(dead)` exceeds `death_proba_threshold` (0.8) for
`death_run_frames` (5) consecutive frames, and is dated to the first frame of
that run. Unscored frames (`NaN`) break a run — the model was never asked about
them, so they cannot manufacture a death.

Death is irreversible, so a run the cell **recovers** from is discarded:
`death_run_frames` consecutive scored frames back under `1 - threshold` after
the run mark it as a transient burst rather than a death. The classifier can
read a cell as dead for a few frames while it rounds up and then correct
itself, and without this guard such a cell is killed on the first frame of its
mitosis and reports a zero-length one.

Where the death frame falls relative to the first episode gives `fate_label`:

| where the death run starts | `fate_label` | reported |
|---|---|---|
| never | `mitotic_survived` | mitotic duration + signal |
| < `min_mitotic_duration_in_frames` after entry | — **excluded** | the cell never had a mitosis; it rounded up because it was dying |
| inside the first episode | `dead_in_mitosis` | death frame, `frames_to_death`, signal over the pre-death mitotic frames |
| after the first episode | `dead_post_mitosis` | full mitotic duration + signal, plus death frame |

Excluding the interphase deaths is what keeps every remaining row meaningful:
each one carries both a mitotic duration and a fluorescence measurement, rather
than a row of `NaN` for a cell that had no mitosis to measure.

**Summary columns**, in the order they are written. `track_start_frame`,
`mitotic_start_frame` and `death_frame` are absolute movie frames and compare
directly; the other frame counts are durations.

| column | meaning |
|---|---|
| `particle` | trackpy particle id |
| `track_length` | frames in the track |
| `n_peaks` | mitotic episodes in the track |
| `track_start_frame` | frame the track begins on |
| `mitotic_start_frame` | frame of mitotic entry (first episode) |
| `frames_to_death` | frames from mitotic entry to the death call; `NaN` if the cell survives |
| `sem_frames_in_mitosis` | every mitotic-labelled frame in the track |
| `corrected_frames_in_mitosis` | of those, the ones before the death call — time in mitosis while still alive |
| `n_scored` | frames the classifier actually scored |
| `dead_cell_score` | mitotic frames the classifier flagged dead |
| `fate_label` | see the table above |
| `death_frame` | frame of the death call; `NaN` if the cell survives |
| `obs_window_after_entry` | frames from mitotic entry to the end of the track — how long the cell *could* have been watched in mitosis |
| `duration_censored` | the track ends while the cell is still mitotic, mid-movie, so its duration is a **lower bound** |
| `<channel>`… | mean signal over `corrected_frames_in_mitosis` |

`sem_frames_in_mitosis` minus `corrected_frames_in_mitosis` is the time the
cell lay dead while still carrying the mitotic label; the two are equal when it
never dies. `frames_to_death` can exceed `corrected_frames_in_mitosis`, since a
post-mitotic death happens after the cell has already left mitosis.

Fluorescence is averaged over the corrected window only, so nothing is measured
from a cell already called dead.

### What never reaches the summary

Eight filters act before or during summarization. They are separate
mechanisms, worth checking in this order when a cell you expect is missing:

| filter | where | control |
|---|---|---|
| track shorter than 10 frames | `tp.filter_stubs` in `track_centroids` | `min_track_length` |
| mitotic detection near the frame edge | `summarize_data` | `border_margin`, `exclude_border_tracks` |
| already mitotic on the track's first frame | `summarize_data` | — |
| starts late in the movie and is mitotic at once | `summarize_data` | `late_track_start_fraction`, `early_mitosis_frames` |
| still mitotic on the movie's last frame | `summarize_data` | — |
| no mitotic run ≥ `min_mitotic_duration` | `summarize_data` | `min_mitotic_duration`, `frame_interval` |
| interphase death (see above) | `summarize_data` | `min_mitotic_duration`, `frame_interval` |
| too little of the track left after mitotic entry | `summarize_data` | `min_window_factor`, `reference_track_fraction`, `reference_duration_frames` |
| track opened late in the movie | `summarize_data` | `max_track_start_fraction` |

The last two are the spurious-track filter, described in [its own
section](#filtering-tracks-whose-duration-cannot-be-trusted) below. Both are
switched off together with `exclude_short_window_tracks = False`.

The first is the one that surprises: a cell that rounds up, dies and loses its
track inside 10 frames is discarded at tracking and cannot be recovered
downstream — it never appears in `*_analysis.xlsx` at all.

The two mitosis-at-the-edge filters use deliberately different tests.

**Starting mitotic disqualifies a track wherever in the movie it begins.** The
segmentation labels anaphase mitotic, so when a cell divides trackpy commonly
opens a fresh particle on a daughter that still carries the mitotic label; its
"mitosis" is the tail of the mother's division, with no entry of its own. These
tracks have characteristically high particle numbers and late start frames,
since they only exist after a division.

That test misses the ones whose first frame or two are not yet labelled, so a
track is **also** rejected when it begins after `late_track_start_fraction`
(default 1/3) of the movie *and* reaches mitosis within `early_mitosis_frames`
(default 3) of its own start. Both conditions are required: a late start on its
own is ordinary, and so is an early mitosis in a track followed from the
outset. This removes a further 10 tracks on a normally cycling position and 37
on an arrest-heavy one.

**Ending mitotic disqualifies a track only if it runs to the movie's last
frame,** where the acquisition cut the episode short. A track that merely stops
mid-movie is trackpy losing the cell, which says nothing about the mitosis.

Note that on a position where cells arrest in mitosis for a long time, the
starting-mitotic filter also removes genuinely arrested cells that were picked
up mid-arrest, not only mislabelled daughters. Their entry is unobserved either
way, so no duration is measurable, but the count of discarded tracks will be
higher than on a normally cycling population.

```python
exp_analysis.files(Path(to_inference_folder), cell_type = "HeLa")
exp_analysis.track_centroids(save_flag = False)
tracks = exp_analysis.measure_signal('GFP', save_flag = False, id = -1)
tracks = exp_analysis.measure_signal('Texas_Red', save_flag = False, id = -1)
tracks = exp_analysis.measure_signal('Cy5', save_flag = True, id = -1) #as needed
summary = exp_analysis.summarize_data(True)
```

### Filtering tracks whose duration cannot be trusted

Two things put a spuriously short mitosis in the summary, and neither is
visible in the duration itself.

**A track watched only briefly after mitotic entry cannot show a long
mitosis**, whatever the cell does. At an observation window of 50 frames, 96%
of tracks report a short mitosis — including ones that began at frame 0 and
were tracked cleanly for hundreds of frames. This is censoring, not a bad
track, but the number it produces is not a duration.

**Trackpy opens fresh particles late in the movie on fragments in crowded
areas.** Holding the window fixed, tracks starting in the last two thirds of
the movie report a short mitosis 24–63% of the time against 2–9% for the rest
— a step, not a gradient. These are the high-`particle`, short-track,
short-mitosis rows.

So `summarize_data` requires, per track:

```
obs_window_after_entry >= min_window_factor * D
track_start_frame      <= max_track_start_fraction * n_frames
```

where `D`, the reference mitotic duration, is **measured from the data**: the
median duration of tracks followed for `reference_track_fraction` of the movie.
It has to be measured rather than set, because no constant serves both a
mitotic arrest (median ~4 h) and an unperturbed control (~40 min) — and neither
does a constant scaled only by the frame interval. The same factors calibrate
themselves to each:

| dataset | reference `D` | → min window | max start | tracks kept |
|---|---|---|---|---|
| CycB OE, 542 frames @ 4 min | 230 fr (15.3 h) | 150 fr | 180 | 213 / 254 |
| unperturbed, 269 frames @ 4 min | 11 fr (44 min) | 8 fr | 89 | 144 / 156 |

Across the whole arrest plate this keeps 69% of tracks, retains 99% of the
well-observed ones, drops the fraction reporting an implausibly short mitosis
from 22.5% to 1.9% (against 2.4% among the well-observed tracks), and leaves
the ranking between wells intact. On the unperturbed plate the same rule keeps
85% and moves neither well's median — which is the point: an artifact filter
should reshape a contaminated distribution and stay out of the way otherwise.
It follows that this filter does little for short-mitosis experiments, where a
3-frame flicker and a real 10-frame mitosis are genuinely hard to tell apart.

**Censoring is reported, not excluded.** `duration_censored` marks a track that
ends while the cell is still mitotic. Dropping those costs a fifth of the data,
moves the median by 0.2 h and makes the short-mitosis rate slightly *worse* — a
censored duration is still a real mitosis watched for a while, unlike the ones
the window test removes.

`self.quality["track_filter"]` records what was calibrated and cut: the
reference set size, `D`, both thresholds, and the row counts before and after.
Worth a look on any position where the yield surprises you.

**Per-position calibration.** `summarize_data` runs one position at a time, so
`D` is measured per position and the threshold varies about ±20% across a
plate. That is noise — positions in a well share their biology. Set
`defaults.reference_duration_frames` to pin one value plate-wide if it matters.
A position with fewer than 20 reference tracks skips the filter and says so.

### Mitotic vs dead discrimination

cellaap labels a rounded-up dying cell mitotic, because in phase contrast it
looks like one. `track_centroids` therefore runs a second classifier over the
detections the semantic segmentation called mitotic and adds three columns to
the cell table:

| column | meaning |
|---|---|
| `mitotic_proba` | P(truly mitotic) |
| `dead_proba` | P(dead-like) = 1 - `mitotic_proba` |
| `dead_flag` | 1 when `dead_proba` > 0.5 |

**The probabilities are computed where the semantic label is mitotic** (100 or
101 - `analysis_pars.mitotic_semantic_values`) **plus
`analysis_pars.post_peak_frames` frames after each episode**. Every other row
keeps `NaN`, which is not the same as a low `dead_proba`: the model was never
asked about those cells. Gating uses the raw semantic label rather than
`semantic_smoothed`, since smoothing both adds frames the segmentation never
called mitotic and drops frames it did.

The tail is what makes post-mitotic death visible: a cell dying on mitotic exit
drops the mitotic label while still rounded, so the frames carrying the death
sit just past the episode.

**Nothing before mitotic entry is ever scored, and the tail must stay short.**
The model is binary *within the rounded-cell population* — class 1 mitotic
against class 0 dead — so on a flat cell class 0 means only "not rounded", not
"dying". Measured on a flat interphase cell that goes on to divide, and is
therefore alive, it returns a mean P(dead) of 0.785 with 82.5% of frames over
0.5, against 0.294 and 27.4% for genuinely mitotic frames. Widening the scored
set past the rounded state therefore manufactures deaths rather than finding
them.

The two probabilities are complementary by construction - this is one binary
model, not two independent scores. Both are reported so downstream code never
has to remember which way round the class encoding runs.

The current model (`models/dead_classifier_pooled.joblib`) is ResNet18-512 ->
PCA(32) -> logistic regression, trained on 826 hand labels pooled from the BUB1
and CycB-overexpression datasets across 15 microscope positions:
leave-one-position-out balanced accuracy 0.878, AUC 0.945, calibration error
0.018. An area-only baseline reaches 0.691, so the model is reading more than
cell size.

Transfer between cell backgrounds is the known weak spot: fitting on one
dataset and testing on the other gives 0.75-0.83 balanced accuracy. A new cell
line will probably need its own hand labels rather than inheriting this model.

The bundle declares which feature blocks it consumes under `feature_set`, and
`classify_dead` builds the design matrix from that declaration, so swapping in
a model with a different layout does not require editing the call sites. A
width mismatch raises rather than silently producing confident nonsense.

To re-score analysis files produced before this model, run
`python augment_dead_label.py <root_folder>`. It rewrites each
`*_analysis.xlsx` in place and then rebuilds the matching `*_summary.xlsx`
from the new labels, since the summary is derived entirely from them. Use
`--suffix` to write copies instead (both files take the suffix), and
`--cell-type` to pick the `analysis_pars` defaults for the rebuilt summary.

`--frame-interval` overrides what the acquisition metadata claims, for the
rebuilt summary.

#### Re-summarizing an analysis file

`from_analysis_file` rebuilds a summary from a saved `*_analysis.xlsx` without
redoing tracking or signal measurement — no image stacks are read, so it takes
seconds rather than the better part of an hour. This is the way to try new
filter settings.

```python
from pathlib import Path
from cellaap_analysis import analysis

f = Path('/path/to/..._inference/..._analysis.xlsx')
analysis.from_analysis_file(f).summarize_data(True)
```

`summarize_data(True)` writes `*_summary.xlsx` beside the analysis file; pass
`False` to just get the DataFrame back. Add `suffix='_test'` to write
`*_summary_test.xlsx` and leave the existing one alone — worth doing for a
first look.

The frame interval comes from the acquisition metadata automatically. To
override it, or to compare against the unfiltered summary:

```python
from pathlib import Path
from cellaap_analysis import analysis

f = Path('/path/to/..._inference/..._analysis.xlsx')
a = analysis.from_analysis_file(f, frame_interval=4)
a.defaults.exclude_short_window_tracks = False   # to compare unfiltered
a.summarize_data(True, suffix='_nofilter')
print(a.quality['track_filter'].T)
```

A whole folder at once:

```python
from pathlib import Path
from cellaap_analysis import analysis

root = Path('/path/to/20593')
for d in sorted(root.glob('*_inference')):
    for f in d.glob('*_analysis.xlsx'):
        analysis.from_analysis_file(f).summarize_data(True, suffix='_test')
```

Two things to check when you do: `a.quality['track_filter']` shows exactly what
was calibrated and cut, and the printed `minimum mitotic duration 30 min = N
frames` line confirms the interval was picked up correctly.

Note that **every summary written before the frame-interval fix used the wrong
minimum mitotic duration**, so re-summarizing is worth doing across the board,
not only where you want the new filter.

**Quality metrics** — `summarize_data` records two, in `self.quality` and in the
spreadsheet's "quality" sheet: a histogram of mitotic episodes per track, and
per-track cell-area standard deviation. Both flag tracking or segmentation
problems; a track with several episodes is usually a cell that divided and then
re-rounded, often to die.

3.**Plotting mode**: One can create multiple objects corresponding, e.g., to multiple repeats of an experiment.

```python
experiment1 = cellapp_analysis.analysis(Path(root_folder_1), plotting_only: True)

experiment2 = cellapp_analysis.analysis(Path(root_folder_2), plotting_only: True)
```

Compiling many positions into one dataframe per condition is no longer part of
this module. It lives in `cellaap_aggregate.py`, described in
[Compiling a plate](#compiling-a-plate-cellaap_aggregatepy) below, and it reads
the plate layout from the same `platemap.csv` that `pipeline.py` ran the plate
from rather than from a hand-written well list.

**Step 6:** Use the **fit_model** function to fit a 4-parameter Hill model to binned data. The function expects input data as a dataframe with the first column containing the fluorescence signal and the second column containing the time in mitosis. For this model to work, the 0 dosage response must be defined as a positive value. This value must be obtained from a -rapamycin well or otherwise supplied. If it is unavailable, perform a rough background subtraction as shown below on a temporary basis.

quant_fraction must a list that specifies the quantiles to be evaluated for the dosage. The bin range is based on the quantile values of the dosage values. Remember that the eSAC dosage distribution is asymmetric (it should be possible to fit it with a log-normal distribution). Therefore, the default quantile values (used below) are asymmetric.

bin_size is arbitrarily defined and can also be adjusted if necessary. Don't use lower values (the default shown below is empirically defined).

```python
# plotting and curve-fitting compiled data
# Note that the first column must be the "dosage" (fluorescence signal) and the second
# column must be the "response" (time in mitosis)
from cellaap_aggregate import fit_model

dose_response = compiled_data.loc[:, ('Texas_Red', 'mitosis')]
# approximate background subtraction if blank well data are unavailable
dose_response.Texas_Red = dose_response.Texas_Red - dose_response.Texas_Red.min()
xy_data, bin_means, bin_stderrs, bin_sizes, fit_values = fit_model(
    dose_response, plot=True, quant_fraction=[0.025, 0.85], bin_size=2.5)
```

## Compiling a plate: `cellaap_aggregate.py`

`cellaap_analysis` turns one position into a `*_summary.xlsx`. `cellaap_aggregate`
is the step after: it collects many of those into one dataframe keyed by what
was in the well, and fits and plots the result.

It does not ask you to describe the plate again. The layout comes from the same
`platemap.csv` that `pipeline.py` ran the plate from, parsed by `pipeline`'s own
code, so a compiled group cannot disagree with what was actually segmented and
analyzed. See [PIPELINE_README.md](PIPELINE_README.md) for the platemap format.

```python
import cellaap_aggregate as agg

df = agg.load_experiment(root_folder, expt_length=150, delta_t=10)
```

That is the whole thing for the common case. `df` holds every position's
summary rows, filtered, with `well`, `position`, `stem`, `celltype`,
`transfection`, `drug`, `code` and `storage_location` added. `code` is the group
label (`HeLa_siBUB1_DMSO`), which is what `export_to_excel_by_col` splits on.

The filters `load_experiment` applies, from `filter_summary`: a row survives
only if mitotic entry was actually observed (`mitotic_start_frame > 0`) and the
whole episode finished inside the movie
(`mitotic_start_frame + corrected_frames_in_mitosis < expt_length`).
`corrected_frames_in_mitosis` is then multiplied by `delta_t`. Pass
`filter_rows=False` for the unfiltered table.

Useful arguments:

| argument | |
| --- | --- |
| `by` | grouping keys, default `("celltype", "transfection", "drug")`. `by=("drug",)` pools every cell type under each drug |
| `suffix` | which summary variant to read: `""`, or `"_dead"` for the copies `augment_dead_label.py --suffix _dead` writes |
| `platemap` | a platemap somewhere other than `<root>/platemap.csv` |
| `verbose` | per-position progress and the platemap warnings, on by default |

### The pieces, if you want to intervene

```python
positions = agg.platemap_positions(root_folder)   # one PositionRef per position
groups    = agg.group_positions(positions)        # keyed by (celltype, transfection, drug)
raw       = agg.compile_positions(groups[('HeLa', 'pEN2', 'DMSO')])
filtered  = agg.filter_summary(raw, expt_length=150, delta_t=10)
```

`platemap_positions` applies exactly the precedence the pipeline applies —
an exact position row (`G03_s9`) beats a well row (`G03`) beats the HeLa
default — and drops `skip` rows and the blank-media wells. A well that is
missing from the platemap is reported and falls back to HeLa, with
`PositionRef.mapped` set to `False`.

If the raw `*phs.tif` stacks have been archived, positions are recovered from
the `*_inference` folder names instead, so a results-only folder still compiles.

### The older well-list API

`create_wellmap_dict`, `compile_summaries`, `import_filter_data_for_wells` and
`import_whole_expt_data` moved here from `cellaap_utils` unchanged in signature,
so existing notebooks need only their import line updated:

```python
from cellaap_aggregate import create_wellmap_dict, import_whole_expt_data

wellmap  = create_wellmap_dict(agg.read_platemap(root_folder))
whole_df = import_whole_expt_data(wellmap, root_folder, expt_length=150, delta_t=10)
```

They now understand the platemap's `well_ids` syntax — ranges (`B01-B06`)
expand, `g3` normalizes to `G03` — and skip the `skip` and blank-media rows.

One hazard remains in this path, and it is why `load_experiment` exists:
`compile_summaries` matches a well id as a substring of the folder name, so if
the platemap gives one site of a well its own row, that site is counted under
**both** conditions. `agg.wellmap_from_platemap(root_folder)` returns the same
mapping resolved to position stubs instead of wells, which fixes it while
keeping the rest of the old call chain:

```python
wellmap  = agg.wellmap_from_platemap(root_folder)
whole_df = import_whole_expt_data(wellmap, root_folder, expt_length=150, delta_t=10)
```

## Curating particles: `particle_browser.py`

A napari browser over the particles a `*_summary.xlsx` lists, for reviewing the
pipeline's calls and dropping the bad ones.

```bash
conda run -n img-env python particle_browser.py "/path/to/root_folder"
```

The folder argument is optional and can also be set from the two folder buttons
in the panel: the inference folder is where the `*_inference` directories live,
the image folder is where the `*.tif` stacks live. They are the same by default
and unticking the checkbox separates them; both are searched one level deep as
a fallback, so pointing at a parent works.

- **Well -> Site -> Particle.** Particles come from the `Summary` sheet; the
  pooled workbook pair is preferred when a position has both. Stacks and
  inference folders pair on the `well_site` key (`A12_s2`), never the file
  stem, since the acquisition and analysis dates can differ.
- **ROI.** A 100x100 crop follows the tracked centroid through every channel
  found for that position, phase at the bottom and fluorescence additive on
  top. `x`/`y` are doubled (half-res segmentation grid -> raw stacks) and
  `frame` indexes the stacks directly. Border ROIs are zero-padded so the cell
  stays centred.
- **Traces.** `semantic`, fluorescence and `dead_proba` overlaid on a shared
  0-1 axis, with each one's true range in the legend; fluorescence is scaled on
  its 1st-99th percentiles so one bright frame cannot flatten it.
  `mitotic_start_frame` and `death_frame` are marked. The cursor follows the napari
  frame slider, and clicking the plot jumps the viewer to that frame.
  - **Pick one fluorescence column or tick `plot all`** to overlay every one
    the analysis table carries, each in its own colour. A column whose range is
    under 1% of its own magnitude - a correction factor that never really moves
    - is drawn dotted and marked `~flat`, so rescaling cannot dress up noise as
    signal.
  - **Legacy workbooks are handled by omission.** A column that is absent, or
    present but all-NaN, is left off the plot rather than drawn as a flat line
    at zero, and the status line names what was dropped. Pre-classifier
    analysis files therefore show semantic plus fluorescence only, with no
    `dead_proba` trace and no `death_frame` marker.
- **Annotation.** A free-text box per particle. Edits are held as you type and
  saved when you leave the box or move on, so switching particles mid-word
  loses nothing; annotated particles are marked `[note]` in the list. Notes go
  out with the data as a `user_annotation` column, and are read back in when a
  previously exported summary is re-opened, so curation resumes where it
  stopped.
- **Exclude / Export.** Exclusions are per position and persisted immediately,
  so scoring survives a restart. Export prompts for a name and writes
  `<name>_summary.xlsx` + `<name>_analysis.xlsx`: the summary keeps every
  original sheet with `Summary` filtered and `user_annotation` added, plus an
  `excluded` sheet listing dropped particles with their notes; the analysis
  keeps the rows of the surviving particles, gaining `user_annotation` only if
  something was annotated. If other positions were also curated it offers to
  export those too, suffixing each pair with its well_site.

Reading a 30 MB analysis workbook takes ~55 s, so each is memoised as parquet
(keyed by size and mtime) and revisits cost ~0.1 s.

**ROI caching.** Selecting a group starts one background pass over the stacks
that fills every particle in it. A pass is dominated by reading planes, and
that cost is independent of how many particles are cut from them, so this
collapses N per-particle reads into one. A 100×100 crop already drags in
~400 KB of pages (100 rows against 16 KB pages) and costs a round trip per
frame — hence ~14 s per channel for a 450-frame track over SMB. A whole plane
costs ~17 crops locally but only ~2–3 over SMB, so a pass pays for itself past
a handful of particles; `BULK_MIN_PARTICLES` skips it for smaller groups, where
per-particle reads are still cheaper.

All channels advance together frame by frame, so a particle is complete the
moment its last frame is read and is handed over while the pass continues.
Anything the pass has not reached yet still loads through the ordinary
on-demand path, so browsing is never blocked.

ROIs are held in memory by total size (1.5 GB, about 80 full-length
two-channel particles) rather than by count, since track lengths vary hugely;
`ROI_CACHE_BYTES` is the dial. Finished particles are also written to
`roi/` under the cache directory, keyed by the stacks' names, sizes and mtimes,
so revisiting a position costs nothing in later sessions.

The worker only touches the disk cache; the GUI thread does every in-memory
insert, so the LRU needs no lock. Speculative entries enter at the evict-first
end, so a prefetch never displaces a particle you actually opened. Each pass
carries a generation, bumped on any position or group change, so superseded
results are discarded on arrival.

Both caches and the curation state live in `~/.cache/particle_browser/`,
outside the repo, since they are machine-local and regenerable.
