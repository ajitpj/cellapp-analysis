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

**Step 4:** Use the **summarize_data** function to create the summary Excel file that lists the average signals measured for all channels, duraion of mitosis, and the correction factors to account for background and excitation intensity variation. Before computing the summary measurements, **any gaps in the semantic label vector are filling by "closing" with a footprint (semantic_footprint) with width equal to the minimum mitotic duration (min_mitotic_duration = 3). Only gaps < 3 frames are filled.**
Any cell that shows multiple peaks in the semantic label vector (after median filtering) is also not summarized.

**Important**: When cells are moving around, it is quite common to lose track of a mitotic cell right after it divides. This leads to a significant number of tracks that end in mitosis. These need not be discarded, especially given that cellapp labels anaphase cells as mitotic making it unlikely that a track ending in mitosis is somehow erroneous. The find_peaks function from scipy ignores peaks that persist till the end. Thus, tracks that end in mitosis are ignored by the find_peaks function. To avoid this, I am adding fantom semantic values at the end of each track (set to 0 - i.e. non-mitotic) before the values are input into the find_peaks function. This recovers the mitotic duration from these tracks.

```python
exp_analysis.files(Path(to_inference_folder), cell_type = "HeLa")
exp_analysis.track_centroids(save_flag = False)
tracks = exp_analysis.measure_signal('GFP', save_flag = False, id = -1)
tracks = exp_analysis.measure_signal('Texas_Red', save_flag = False, id = -1)
tracks = exp_analysis.measure_signal('Cy5', save_flag = True, id = -1) #as needed
summary = exp_analysis.summarize_data(True)
```

**Quality metrics** - The *summarize_data* function calculate wo simple quality metrics: the number of peaks per cell track and fluctuations in cell area (standard deviation). These are stored in the *self.quality* dictionary. The number of peaks per cell track is summarized as a histogram in the excel spreadsheet in the sheet labeled "Quality". Cell area standard deviation is reported as a column vector.

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
