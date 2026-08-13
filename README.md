# cellapp-analysis module

A python module to track cells and measure fluorescence and mitotic duration using cellapp inference files and raw fluorescence images. It runs trackpy on the instance segmentation file first, and then measures cell state (mitotic/non-mitotic) from the semantic segmentation and raw mean fluorescence value from the specific fluorescence channel. The module also summarizes the analysis by calculating the mean fluorescence for each cell over the duration of mitosis, the total mitotic duration, and correction factors based on intensity and background correction maps (must be acquired using empty wells with DMEM and fluorobrite respectively).
This analysis and summary are saved as separate excel spreadsheets.

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
6. min_mitotic_duration = 3: (unrelated to trackpy); mitotic events smaller than 3 timepoints are filtered out.

Be careful when using the "predictive" tracking mode. It's very powerful, but can be computationally costly if the tracking memory>1 and max. pixel movement is > 25. This will lead to trackpy exceeding the max. number of nodes in one or more subnetworks. Currenlty, trackpy just exits on this error, which can be problematic when analysis is being done in batch mode. If you come across this issue, gradually decrease the max. pixel movement parameter to get under this error.

**Step 1:** Point the analysis object to a specific inference folder by providing a path to it. This will read the instance segmentation file from this folder.

**Step 2:** Use the **track_centroids** method; it will erode the instance segmentation with the default footprint (needs to be customized for different cells), track the resultant masks using trackpy, and then determine the cell-state by reading the semantic stack. The intermediate padas dataframe can be saved to excel using the flag. See the notes above regarding optimal tracking parameters. 

*It seems that trackpy is using 32-bit integers for assigning labels to individual points. Overrun of this number leads to missing signal measurements, which shouldn't be a big issue. But this needs to be adddressed at some point.*

**Step 3:** Use the **measure_signal** function to measure the fluorescence from the specified channel. The channel string must match the channel name in the file names. The "id = -1" will make the function measure data for all cells that went through a complete mitosis during the time lapse. Optionally, one can provide a list with cell numbers (development only). Thus, cells that remained in interphase throughout the experiment are not measured. Their tracks are still reported.

**Step 4:** Use the **summarize_data** function to create the summary Excel file that lists the average signals measured for all channels, duraion of mitosis, and the correction factors to account for background and excitation intensity variation. Before computing the summary measurements, **gaps in the semantic label vector are filled by "closing" with a footprint (semantic_footprint) of width min_mitotic_duration = 3, so only gaps < 3 frames are filled. The median filter and the closing are applied per particle, over that track's own frames in frame order** — run over the whole table they would bleed across track boundaries, letting the end of one cell's trace close a gap at the start of the next.
### How a cell is summarized

Two rules, per track. No smoothing-vs-classifier precedence, no state decode —
those were replaced by this scheme.

**Mitotic episodes** are the runs of mitotic semantic label at least
`min_mitotic_duration_in_frames` long, read straight off the trace. The first
episode is the mitosis that gets reported; `n_peaks` says how many there were.
A later episode is usually the same cell re-rounding, most often to die.

**Death** is called when `P(dead)` exceeds `death_proba_threshold` (0.8) for
`death_run_frames` (5) consecutive frames, and is dated to the first frame of
that run. Unscored frames (`NaN`) break a run — the model was never asked about
them, so they cannot manufacture a death.

Death is irreversible, so a run the cell **recovers** from is discarded —
`death_run_frames` consecutive scored frames back under `1 - threshold` after
the run means it was a transient burst, not a death. The classifier can read a
cell as dead for a few frames while it rounds up and then correct itself:
E10_s7 particle 87 died on the first frame of its mitosis at P(dead) > 0.85,
was back at 0.01 five frames later, and ran a second mitosis 300 frames on at
P(dead) = 0.00. Without the guard it reported `mitosis = 0`.

Where the death frame falls relative to the first episode gives `fate_label`:

| where the death run starts | `fate_label` | reported |
|---|---|---|
| never | `mitotic_survived` | mitotic duration + signal |
| < `min_mitotic_duration_in_frames` after entry | — **excluded** | the cell never had a mitosis; it rounded up because it was dying |
| inside the first episode | `dead_in_mitosis` | death frame, `time_to_death`, signal over the pre-death mitotic frames |
| after the first episode | `dead_post_mitosis` | full mitotic duration + signal, plus death frame |

Excluding the interphase deaths is what keeps every remaining row meaningful:
each one carries both a mitotic duration and a fluorescence measurement, rather
than a row of `NaN` for a cell that had no mitosis to measure.

**Summary columns.** `mito_start` and `death_frame` are absolute movie frames
and compare directly; everything else counted in frames is a duration.

| column | meaning |
|---|---|
| `mito_start` | frame of mitotic entry (first episode) |
| `mitosis` | length of the first episode, as observed, never truncated at death |
| `time_to_death` | frames from mitotic entry to the death call; `NaN` if the cell survives |
| `death_frame` | frame of the death call; `NaN` if the cell survives |
| `fate_label` | see the table above |
| `n_peaks` | mitotic episodes in the track |
| `n_sem_mitotic` | frames the segmentation called mitotic |
| `n_scored` | frames the classifier actually scored |
| `<channel>` | mean signal over the mitotic window, up to the death call |

### What never reaches the summary

Three filters act before or during summarization, and they are separate
mechanisms — worth checking in this order when a cell you expect is missing:

| filter | where | control |
|---|---|---|
| track shorter than 10 frames | `tp.filter_stubs` in `track_centroids` | `min_track_length` |
| mitotic detection near the frame edge | `summarize_data` | `border_margin`, `exclude_border_tracks` |
| mitotic on the movie's first or last frame | `summarize_data` | — |
| no mitotic run ≥ 3 frames | `summarize_data` | `min_mitotic_duration_in_frames` |
| interphase death (see above) | `summarize_data` | `min_mitotic_duration_in_frames` |

The first is the one that surprises: a cell that rounds up, dies and loses its
track inside 10 frames is discarded at tracking and cannot be recovered
downstream — it never appears in `*_analysis.xlsx` at all. On E10_s7 the border
and short-episode filters removed 51 and 53 of 497 mitotic tracks; the 53 had a
median of 2 mitotic-labeled frames, i.e. sub-threshold roundings.

A cell mitotic on the **movie's** first or last frame had its entry or exit
clipped by the acquisition, so no duration is measurable and it is dropped (61
tracks on E10_s7). The test is deliberately against the movie bounds and not
the track's own ends: tracks routinely start and stop mid-movie when trackpy
loses a rounding cell and re-acquires it as a new particle, and excluding those
would discard 325 of 497 mitotic tracks rather than 73.

```python
exp_analysis.files(Path(to_inference_folder), cell_type = "HeLa")
exp_analysis.track_centroids(save_flag = False)
tracks = exp_analysis.measure_signal('GFP', save_flag = False, id = -1)
tracks = exp_analysis.measure_signal('Texas_Red', save_flag = False, id = -1)
tracks = exp_analysis.measure_signal('Cy5', save_flag = True, id = -1) #as needed
summary = exp_analysis.summarize_data(True)
```

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
"dying". Measured on A12_s2: frames before mitotic entry, from cells that go on
to divide and are therefore alive, score mean P(dead) = 0.785 with 82.5% over
0.5, against 0.294/27.4% for mitotic frames. Scoring every frame to the end of
the track was tried on that basis and reverted — it moved `dead_post_mitosis`
1 → 37 and `mitotic_survived` 114 → 27, almost entirely artifact.

The two probabilities are complementary by construction - this is one binary
model, not two independent scores. Both are reported so downstream code never
has to remember which way round the class encoding runs.

The current model (`models/dead_classifier_pooled.joblib`) is ResNet18-512 ->
PCA(32) -> logistic regression, trained on 826 hand labels pooled from the BUB1
and CycB-overexpression datasets across 15 microscope positions:
leave-one-position-out balanced accuracy 0.878, AUC 0.945, calibration error
0.018. An area-only baseline reaches 0.691, so the model is not simply
measuring cell size - the failure mode of the first version of it.

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

The same path is available directly for re-summarizing an analysis file
without redoing tracking or signal measurement — no image stacks are read, so
it takes seconds rather than the better part of an hour:

```python
from cellaap_analysis import analysis
analysis.from_analysis_file(path_to_analysis_xlsx, cell_type="hela").summarize_data(True)
```

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

In this mode, the module is used to compile all data corresponding to positions and wells belonging to the same treatment/cell line into one dataframe. e.g.,

**Step 5:** Use the **compile_summaries** function to collect multiple wells and/or positions that represent the same experiment. The function requires a list as the input. Each entry in the list must be a string encoding the well and position identifier. Notice the capitalization and well number convention used in the example below. The output is a dataframe with an additional column for the well+position designation. In the future, one more column indicating the experiment will be added.

```python
HeLa_wells = ["A02", "H12"] # Note the exact formant
HeLa_data = experiment1.compile_summaries(HeLa_wells)  
```

**Step 6:** Use the **fit_model** function to fit a 4-parameter Hill model to binned data. The function expects input data as a dataframe with the first column containing the fluorescence signal and the second column containing the time in mitosis. For this model to work, the 0 dosage response must be defined as a positive value. This value must be obtained from a -rapamycin well or otherwise supplied. If it is unavailable, perform a rough background subtraction as shown below on a temporary basis.

quant_fraction must a list that specifies the quantiles to be evaluated for the dosage. The bin range is based on the quantile values of the dosage values. Remember that the eSAC dosage distribution is asymmetric (it should be possible to fit it with a log-normal distribution). Therefore, the default quantile values (used below) are asymmetric.

bin_size is arbitrarily defined and can also be adjusted if necessary. Don't use lower values (the default shown below is empirically defined).

```python
# plotting and curve-fitting compiled data
# Note that the first column must be the "dosage" (fluorescence signal) and the second
# column must be the "response" (time in mitosis)
dose_response = compiled_data.loc[:, ('Texas_Red', 'mitosis')]
# approximate background subtraction if blank well data are unavailable
dose_response.Texas_Red = dose_response.Texas_Red - dose_response.Texas_Red.min()
xy_data, bin_means, fit_values = fit_model(dose_response, plot: True, 
                                           quant_fraction = [0.025, 0.85], 
                                           bin_size = 2.5)
```

Wellmap and experiment import utilities
(Generated by GPT5-mini)
Two helper functions in `cellaap_utils.py` simplify importing multiple wells and
grouping them by experimental metadata:

- `create_wellmap_dict(imported_wellmap: pd.DataFrame, wellid_col_name='well_ids')`
- `import_whole_expt_data(wellmap_dict: dict, analysis_object, expt_length: int, delta_t: int)`

Below are short usage examples and notes.

### create_wellmap_dict

Purpose: turn a user-provided well-map table (pandas DataFrame) into a
dictionary that maps grouping keys (typically `(celltype, transfection, drug)`)
to an ordered list of well identifiers.

Expected input DataFrame columns: `celltype`, `transfection`, `drug`, and a
column (default name `well_ids`) containing well identifiers. The `well_ids`
column may contain lists/tuples or strings with wells separated by commas,
semicolons, or whitespace.

Example:

```python
import pandas as pd
from cellaap_analysis.cellaap_utils import create_wellmap_dict

# a small example table with two rows describing the same group
df = pd.DataFrame([
    {'celltype': 'HeLa', 'transfection': 'plasmidA', 'drug': 'DMSO', 'well_ids': 'A01,A02'},
    {'celltype': 'HeLa', 'transfection': 'plasmidA', 'drug': 'DMSO', 'well_ids': ['B01']},
])

wellmap = create_wellmap_dict(df, wellid_col_name='well_ids')
# example key: ('HeLa', 'plasmidA', 'DMSO') -> ['A01', 'A02', 'B01']
```

Notes:

- The function preserves the first-seen order of wells and removes duplicates.
- It raises `KeyError` if required columns are missing.

### import_whole_expt_data

Purpose: iterate over a `wellmap_dict` (from `create_wellmap_dict`) and import
and filter summary spreadsheets for every well group. This produces one
concatenated `pandas.DataFrame` suitable for downstream analysis and plotting.

Inputs:

- `wellmap_dict`: mapping produced by `create_wellmap_dict`
- `analysis_object`: your analysis/session object (the same object used
  elsewhere in the package; it must provide the `root_folder` used by
  `compile_summaries`)
- `expt_length`: total number of frames in the movie (used to filter events)
- `delta_t`: time per frame (used to convert `mitosis` frames into time)

Example:

```python
from cellaap_analysis.cellaap_utils import import_whole_expt_data

# wellmap was produced above via create_wellmap_dict
whole_df = import_whole_expt_data(wellmap, analysis_object=exp_analysis,
                                  expt_length=600, delta_t=10)

# `whole_df` is a DataFrame containing filtered summary rows and a `code`
# column (e.g. 'HeLa_plasmidA_DMSO') that identifies the group.
```

Notes:

- Empty well lists are skipped and exceptions during a group's import are
  printed but do not stop processing of other groups.

- If no data is found the function returns an empty DataFrame.

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
  `mito_start` and `death_frame` are marked. The cursor follows the napari
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
(keyed by size and mtime) and revisits cost ~0.1 s. ROI reads memory-map the
stacks and pull only the crop -- ~14 s per channel for a 450-frame track over
SMB -- so they are cached in memory by total size (1.5 GB, about 80 full-length
two-channel particles) rather than by count, since track lengths vary hugely;
a revisit is then instant. `ROI_CACHE_BYTES` at the top of the file is the dial.
Both the parquet cache and the curation state live in
`~/.cache/particle_browser/`, outside the repo, since they are machine-local
and regenerable.
