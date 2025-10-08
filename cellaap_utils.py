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


def calculate_signal(semantic, signal, bkg_corr, int_corr, footprint):
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

    return signal_mean, bkg_corr_mean, int_corr_mean, signal_std, bkg_corr_std, int_corr_std


def calculate_displacement(coords: pd.DataFrame) -> float:
    '''
    Function calculates the absolute displacement from frame to frame
    '''
    pixel_shift = coords.diff()
    displacement = np.sqrt(pixel_shift.x**2 + pixel_shift.y**2)

    return displacement


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
    bins   = np.arange(0.5*quants[0], 1.5*quants[-1], bin_size)

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
        y_quant = np.round(xy_data.iloc[:,1].quantile(0.99))
        ax.set_ylim([0, y_quant])
    
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


def filter_summary_df(summary_df: pd.DataFrame, end_frame: int) -> pd.DataFrame:
    '''
    Function removes summary entries wherein mitosis starts at 0 or
    ends in the last frame of the movie.
    '''
    filtered_df = summary_df[summary_df["mito_start"] > 0].copy()
    filtered_df = filtered_df[(filtered_df["mitosis"] + filtered_df["mito_start"] < end_frame)].copy()

    return filtered_df