import numpy.typing as npt
from skimage.filters import gaussian
from skimage.morphology import closing
import numpy as np
import scipy.ndimage as ndi
from scipy.signal import medfilt
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns
import re

def projection(im_array: np.ndarray, projection_type: str):
    """
    Compute a projection of a 3D image stack along the first axis.

    The function selects a central slab of slices (half of the central index)
    and computes either the maximum, minimum, or average projection across
    that slab.

    Parameters
    ----------
    im_array : np.ndarray
        A 3D image stack with shape (z, y, x).
    projection_type : str
        One of "max", "min", or "average" to control the type of
        projection applied to the central slab.

    Returns
    -------
    np.ndarray
        2D projected image.
    """

    if im_array.shape[0] % 2 == 0:
        center_index = im_array.shape[0] // 2 - 1
    else:
        center_index = im_array.shape[0] // 2

    range = center_index // 2

    try:
        assert projection_type in ["max", "min", "average"]
    except AssertionError:
        print("Projection type was not valid, valid types include: max, min, mean")

    if projection_type == "max":
        projected_image = np.max(
            im_array[center_index - range : center_index + range], axis=0
        )
    elif projection_type == "average":
        projected_image = np.mean(
            im_array[center_index - range : center_index + range], axis=0
        )
    elif projection_type == "min":
        projected_image = np.min(
            im_array[center_index - range : center_index + range], axis=0
        )

    return np.array(projected_image)
    
def gen_intensity_correction_map(image: npt.NDArray) -> npt.NDArray:
    """
    From Anish
    Computes the intensity map for flouresence microscopy intensity normalization if the input is a blank with flourescent media
    ----------------------------------------------------------------------------------------------------------------------------
    INPUTS:
        image: npt.NDArray
    OUTPUTPS:
        intensity_map: npt.NDArray
    """
    mean_plane = projection(image, "average")
    # med_filtered_mean_plane = ndi.median_filter(mean_plane, 9)
    smoothed_mean_plane = gaussian(mean_plane, 45)
    intensity_correction_map = smoothed_mean_plane / (np.max(smoothed_mean_plane))

    return intensity_correction_map

def gen_background_correction_map(background_stack: npt.NDArray) -> npt.NDArray:
    '''
    Newly written to avoid too much smoothing. The cMOS camera has a persistent noise pattern
    therefore, it is better to keep the corrections local. 
    '''

    background_correction_map = np.zeros_like(background_stack, dtype=int)
    footprint = footprint=np.ones((3,3))
    for i in np.arange( background_stack.shape[0]):
        background_correction_map[i,:,:] = ndi.median_filter(background_stack[i,:,:], footprint=footprint)

    return background_correction_map


def mean_signal_from_mask(img: npt.NDArray, mask: npt.NDArray):
    """
    Compute the mean intensity of pixels within a boolean mask.

    Parameters
    ----------
    img : npt.NDArray
        2D image from which to sample pixel values.
    mask : npt.NDArray
        Boolean mask of the same shape as `img`. True values indicate
        pixels to include in the mean.

    Returns
    -------
    float
        Mean intensity of the selected pixels. Returns NaN if the mask
        selects no pixels.
    """
    pixels = img[np.nonzero(mask)]
    if pixels.any():
        mean_signal = np.mean(pixels)
    else:
        mean_signal = np.nan

    return mean_signal


def calculate_signal(semantic, signal, bkg_corr, int_corr, area, footprint):
    '''
    utility function for calculating signal from the given semantic, signal, and bkg traces
    '''
    
    if signal.any():
        signal_mean = np.nanmean(signal[np.where(semantic)])
        signal_std = np.nanstd(signal[np.where(semantic)])
    else:
        signal_mean = 0
        signal_std = 0
    
    if bkg_corr.any():
        bkg_corr_mean = np.nanmean(bkg_corr[np.where(semantic)])
        bkg_corr_std = np.nanstd(bkg_corr[np.where(semantic)])
    else:
        bkg_corr_mean = 0
        bkg_corr_std = 0

    if int_corr.any():
        int_corr_mean = np.nanmean(int_corr[np.where(semantic)])
        int_corr_std = np.nanstd(int_corr[np.where(semantic)])
    else:
        int_corr_mean = 1
        int_corr_std = 0

    if area.any():
        area_mean = np.nanmean(area[np.where(semantic)])
        area_std = np.nanstd(area[np.where(semantic)])
    else:
        area_mean = 1
        area_std = 0

    return signal_mean, bkg_corr_mean, int_corr_mean, area_mean, signal_std, bkg_corr_std, int_corr_std, area_std


def calculate_displacement(coords: pd.DataFrame) -> pd.Series:
    '''
    Function calculates the absolute displacement from frame to frame
    '''
    pixel_shift = coords.diff()
    displacement = np.sqrt(pixel_shift.iloc[:, 0] ** 2 + pixel_shift.iloc[:, 1] ** 2)
    #

    return pd.Series(displacement, index=coords.index)


def fit_model(xy_data: pd.DataFrame, plot: True, quant_fraction = None, bin_size = None) -> (pd.DataFrame, dict): # type: ignore
    '''
    Function to fit the dose-response data with a 4-parameter sigmoid.
    Bin range is determined by quantiles. Default is 0.025 and 0.85. The data
    typically contain outliers on the high side, but not the low side. Hence the 
    default values are aysmmetric. For the model to be applicable, the fluorescence 
    signal must be background subtracted. A simple method is to subtract the smallest
    signal value from all values.
    
    Inputs:
    xy_data        - dataframe w/ dose as the first column and response as 
                     the second column
    plot           - Boolean to enable plotting
    quant_fraction - quantiles to determine bin range; 
    bin_size       - size of each bin, default is 2.5 (empirical)
    Outputs:
    xy_data        - the input dataframe with bin labels added as a new column
    fit_pars       - dictionary containing fit parameters
    '''

    xy_data.dropna(inplace=True)
    if quant_fraction is None:
        quant_fraction = [0.025, 0.85]
    quants = np.round(xy_data.iloc[:,0].quantile(quant_fraction)).tolist()

    # 
    if bin_size is None:
        bin_size = 2.5
    bins   = np.arange(0.5*quants[0], 1.5*quants[-1], bin_size).tolist()

    labels, _ = pd.cut(xy_data.iloc[:, 0], bins, retbins=True)
    xy_data["bins"] = labels
    
    bin_means = xy_data.groupby("bins").mean()
    bin_sizes = xy_data.groupby("bins").size()
    bin_stderrs = xy_data.groupby("bins").std()
    bin_stderrs.iloc[:,0] /= bin_sizes**0.5
    bin_stderrs.iloc[:,1] /= bin_sizes**0.5
    bin_means.dropna(inplace=True) # Some of the bins may not have any data
    bin_stderrs.dropna(inplace=True)
    

    fits, _ = curve_fit(sigmoid_4par, bin_means.iloc[:,0], bin_means.iloc[:,1], 
                        p0 = [bin_means.iloc[:,1].min(), bin_means.iloc[:,1].max(), 
                              5, (quants[0] + quants[-1])/ 4
                             ],
                        # sigma = bin_stderrs.iloc[:,1].to_numpy(),
                        maxfev = 10000
                       )

    fit_values = { 'min_duration' : fits[0],
                   'max_duration' : fits[1],
                   'Hill_exponent': fits[2],
                   'EC50'         : fits[3]
                 }
    if plot:
        fig, ax = plt.subplots(1,1, figsize=(8,6))
        sns.scatterplot(x=xy_data.iloc[:,0], y=xy_data.iloc[:,1], 
                        ax=ax, alpha=0.1, 
                        color="gray", edgecolor="None", size=1, 
                        )

        sns.scatterplot(x = bin_means.iloc[:,0], y = bin_means.iloc[:,1], 
                        color='w', edgecolor="blue", marker='s', linewidth=1,
                        label = "binned mean values")
        
        x_range = np.arange(0,1.5*quants[-1])
        sns.lineplot(x=x_range, y=sigmoid_4par(x_range,
                                               fit_values['min_duration'],
                                               fit_values['max_duration'],
                                               fit_values["Hill_exponent"],
                                               fit_values["EC50"]),
                                               ax = ax,
                                               markers='',
                                               color='b',
                                               label="Hill sigmoid fit")
        ax.set_xlabel("eSAC dosage (a.u.)")
        ax.set_ylabel("Time in mitosis (x 10 min)")
        ax.set_xlim(xmax=x_range[-1], xmin=x_range[0])
        y_quant = float(np.round(xy_data.iloc[:,1].quantile(0.99)))
        ax.set_ylim(0.0, y_quant)
    
    return xy_data, bin_means, bin_stderrs, bin_sizes, fit_values

def sigmoid_4par(x, base, top, exponent, ec50):

    """
    4-parameter Hill (sigmoid) function.

    This implements a common dose-response parameterization where
    `base` and `top` are the lower and upper asymptotes, `exponent`
    is the Hill coefficient, and `ec50` is the x value at half-max.

    Parameters
    ----------
    x : array-like or float
        Independent variable(s).
    base : float
        Minimum (baseline) value of the function.
    top : float
        Maximum (top) value of the function.
    exponent : float
        Hill coefficient (controls slope/steepness).
    ec50 : float
        Half-maximal effective concentration (EC50).

    Returns
    -------
    array-like or float
        The evaluated sigmoid at `x`.
    """

    return base + (top - base)*(x**exponent)/(x**exponent+ec50**exponent)


def create_wellmap_dict(imported_wellmap: pd.DataFrame, wellid_col_name = 'well_ids'):
    """Create a mapping from (celltype, transfection, drug) -> list of well ids.
    generated by GPT5-mini
    The function accepts `wellmap_df` (a pandas DataFrame) and a `wellid_col` column name which can contain
    either a Python list/tuple of wells or a string with comma/semicolon/space-separated wells.

    Returns: dict keyed by (celltype, transfection, drug) -> list of unique well ids (strings).
    """
    import re as _re

    # validate required columns
    for c in ('celltype','transfection','drug'):
        if c not in imported_wellmap.columns:
            raise KeyError(f'Required column not found in imported_wellmap: {c}')
    if wellid_col_name not in imported_wellmap.columns:
        raise KeyError(f'{wellid_col_name} column not found in imported_wellmap')

    def _parse_wells(v):
        # return list of wells from various possible representations
        if pd.isna(v):
            return []
        if isinstance(v, (list, tuple)):
            return [str(x).strip() for x in v if str(x).strip()]
        s = str(v)
        # split on comma, semicolon, or whitespace
        parts = [p.strip() for p in _re.split(r'[,;\s]+', s) if p.strip()]
        return parts

    wellmap_dict = {}
    grouped = imported_wellmap.groupby(['celltype','transfection','drug'])
    for key, group in grouped:
        wells = []
        for v in group[wellid_col_name]:
            wells.extend(_parse_wells(v))
        # keep order but unique
        seen = {}
        unique_wells = []
        for w in wells:
            if w not in seen:
                seen[w] = True
                unique_wells.append(w)
        wellmap_dict[key] = unique_wells

    return wellmap_dict

def compile_summaries(cellapp_expt, wells: list) -> pd.DataFrame:
        '''
        Collects and concatenates the summary xslx spreadsheets from the designated well_position list.
        
        Inputs:
        well : list with entries of the form r"[A-Z][dd]+_+[a-z][d+]"
        Output: 
        data_summary  : dataframe with the data concatenated; well - column designating well_pos
        '''
        if wells:
            if not hasattr(cellapp_expt, 'inf_folder_list'):
                # Assemble the folder list when function called for the first time
                cellapp_expt.inf_folder_list = [f for f in cellapp_expt.root_folder.glob('*_inference')]

            df_list = []
            storage_location = cellapp_expt.root_folder.parent

            pattern = re.compile(r'_(\w\d+)_(\w\d+)_')

            for well in wells:
                wp_string = '_'+well+'_'
                for f in cellapp_expt.inf_folder_list:
                    if wp_string in f.name:
                        xls_file_name = [file for file in f.glob('*_summary.xlsx')]
                        if xls_file_name:
                            df = pd.DataFrame()
                            df = pd.read_excel(xls_file_name[0]) #assumes only one
                            df["well"] = well #assign well-position identifier
                            position = pattern.findall(xls_file_name[0].name)[0][1]
                            df["position"] = position
                            df["storage_location"] = storage_location
                            df_list.append(df)
                            print(f"{well}_{position} loaded")
                            del df
                        else:
                            print(f"No summary file found for {well}")

            data_summary = pd.concat(df_list)

        return data_summary

def import_filter_data_for_wells(analysis_object, expt_label: str, expt_length: int, delta_t: int, well_list: list) -> pd.DataFrame:
    """
    Import summary spreadsheets for a list of wells, apply basic quality filters,
    and normalize time units.

    This function wraps `compile_summaries` to load per-well summary spreadsheets
    (assembled from the analysis object's inference folders), then applies the
    following processing steps:

    - Removes any records where `mito_start` is less than or equal to 0.
    - Removes records where the mitosis event would extend beyond the
      experiment length (i.e. `mito_start + mitosis >= expt_length`).
    - Converts the `mitosis` duration from frames to real time by multiplying
      by `delta_t`.
    - Adds a `code` column set to `expt_label` so rows can be identified by
      experiment.

    Parameters
    ----------
    analysis_object : object
        An analysis/session object that provides experiment folder information
        consumed by `compile_summaries`. Expected to have at least a
        `root_folder` attribute; see `compile_summaries` for details.
    expt_label : str
        Short label used to tag the imported rows in the `code` column.
    expt_length : int
        Number of frames in the full experiment/movie. Used to remove events
        that end after the movie finishes.
    delta_t : int | float
        Time per frame (units depend on the dataset). `mitosis` values in the
        source summaries are assumed to be in frames; they will be multiplied
        by `delta_t` to convert to time units.
    well_list : list
        Sequence of well identifiers (strings) to import. Passed to
        `compile_summaries`.

    Returns
    -------
    pandas.DataFrame
        Filtered summary table for the requested wells with columns at least
        including `mito_start`, `mitosis`, `well`, `position`, `storage_location`,
        and `code`. If no data is found for the provided wells, an empty
        DataFrame is returned.

    Notes
    -----
    The actual file I/O and summary assembly is performed by
    `compile_summaries`; this function only post-processes that result.
    """
    well_data = compile_summaries(analysis_object, well_list)
    well_data = well_data[well_data["mito_start"]>0]
    well_data = well_data[well_data["mito_start"]+well_data["mitosis"] < expt_length].copy()
    well_data["mitosis"] = well_data["mitosis"] * delta_t
    well_data["code"] = expt_label
    return well_data


def import_whole_expt_data(wellmap_dict: dict, analysis_object, expt_length: int, delta_t: int) -> pd.DataFrame:
    """
    Load and aggregate filtered summary data for every well group described by a
    well-map dictionary.

    This function iterates over the `wellmap_dict` which should map keys like
    (celltype, transfection, drug) to a sequence of well identifiers. For each
    key it calls `import_filter_data_for_wells` to import and filter the
    per-well summary tables, tags rows using a generated experiment name, and
    concatenates all results into a single DataFrame.

    Parameters
    ----------
    wellmap_dict : dict
        Mapping from grouping keys (typically tuples such as
        `(celltype, transfection, drug)`) to an iterable of well identifier
        strings. If `None` or empty, an empty DataFrame is returned.
    analysis_object : object
        Analysis/session object passed down to `import_filter_data_for_wells`.
        Expected to provide experiment folder information (see
        `compile_summaries`).
    expt_length : int
        Number of frames in the experiment; used to filter events that extend
        past the end of the movie.
    delta_t : int | float
        Time per frame used to convert `mitosis` from frames to time units.

    Returns
    -------
    pandas.DataFrame
        Concatenated and filtered summary table for all wells in the map. The
        returned DataFrame contains an added `code` column constructed from the
        grouping key (e.g. "celltype_transfection_drug"). If no wells are
        available, an empty DataFrame is returned.

    Notes
    -----
    - The function prints progress and skips keys with empty well lists.
    - Exceptions during import of a particular group are caught and logged
      (printed) so processing of other groups continues.
    """
    
    whole_expt_data = pd.DataFrame()

    for key, wells in (wellmap_dict or {}).items():
        try:
            if not wells:
                print('Skipping', key, '- no wells')
                continue
            # construct an experiment name for logging / code column
            expt_name = f'{key[0]}_{key[1]}_{key[2]}'
            print(expt_name)
            temp_df = import_filter_data_for_wells(analysis_object, expt_name, expt_length, delta_t, wells)
            # print(temp_df.shape)
            whole_expt_data = pd.concat([whole_expt_data, temp_df], ignore_index=True)
            print(f'Loaded {wellmap_dict}', key, '->', len(temp_df), 'rows')
            del temp_df

        except Exception as e:
            print(f'Failed to load {wellmap_dict}', key, type(e).__name__, e)

    return whole_expt_data