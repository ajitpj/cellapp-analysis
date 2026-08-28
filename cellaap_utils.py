"""Measurement primitives for `cellaap_analysis`.

Everything here is per-image or per-track arithmetic that the analysis stage
calls while it works through a single position. `cellaap_analysis` pulls it in
with `from cellaap_utils import *`, so `__all__` below is the contract between
the two modules; adding to it widens that contract.

Deliberately light on imports. This module is loaded by every analysis array
task on the cluster, so it holds nothing that needs matplotlib, seaborn or
scipy.optimize. Compiling a plate into per-condition tables, fitting
dose-response curves and plotting them live in `cellaap_aggregate.py`, which is
a notebook module and imports what it likes.
"""

import re
import numpy.typing as npt
from skimage.filters import gaussian
# from skimage.morphology import closing
import numpy as np
import scipy.ndimage as ndi
# from scipy.signal import medfilt
import pandas as pd
from pathlib import Path

__all__ = [
    "projection",
    "gen_intensity_correction_map",
    "gen_background_correction_map",
    "mean_signal_from_mask",
    "window_stats",
    "calculate_signal",
    "calculate_displacement",
    "read_frame_interval",
]


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


def window_stats(window, values, empty_mean=0.0):
    '''
    Mean and standard deviation of one per-frame trace over a track's window.

    The generic form of `calculate_signal` below, which is fixed at four
    named traces. `summarize_data` needs to average whatever per-channel
    columns happen to be present - the raw signal, the two map corrections,
    and now the position-specific corrected signal - so it loops over columns
    and calls this instead of taking a fixed tuple apart.

    `empty_mean` is what an all-zero trace reports: 0 for something that is
    added or subtracted, 1 for something that divides, which is what keeps a
    missing correction from zeroing a signal downstream.
    '''
    values = np.asarray(values)
    if values.any():
        selected = values[np.where(window)]
        return np.nanmean(selected), np.nanstd(selected)
    return empty_mean, 0.0


def calculate_signal(semantic, signal, bkg_corr, int_corr, area, footprint):
    '''
    utility function for calculating signal from the given semantic, signal, and bkg traces
    '''

    signal_mean, signal_std     = window_stats(semantic, signal, 0)
    bkg_corr_mean, bkg_corr_std = window_stats(semantic, bkg_corr, 0)
    int_corr_mean, int_corr_std = window_stats(semantic, int_corr, 1)
    area_mean, area_std         = window_stats(semantic, area, 1)

    return signal_mean, bkg_corr_mean, int_corr_mean, area_mean, signal_std, bkg_corr_std, int_corr_std, area_std


def calculate_displacement(coords: pd.DataFrame) -> pd.Series:
    '''
    Function calculates the absolute displacement from frame to frame
    '''
    pixel_shift = coords.diff()
    displacement = np.sqrt(pixel_shift.iloc[:, 0] ** 2 + pixel_shift.iloc[:, 1] ** 2)
    #

    return pd.Series(displacement, index=coords.index)


def read_frame_interval(folder, pattern: str = "*metadata*.txt") -> float:
    """Minutes between frames, read from the acquisition metadata.

    The interval is a property of the acquisition, not of the analysis, so it
    is read from the microscope's own metadata rather than carried as a default
    in `analysis_pars`. A wrong value is silent and expensive: it feeds
    `min_mitotic_duration_in_frames`, so it decides what counts as a mitotic
    episode at all, and every duration reported downstream is scaled by it.

    The files are written one per wavelength (`*_w1_metadata.txt`, `_w2`, ...)
    and carry a `Time interval:4 min` line. All of them are read and required to
    agree - they describe one acquisition, so a disagreement means the folder
    holds more than one and no single interval is correct.

    The value is validated rather than trusted. Real acquisitions have written
    `Time interval:-26 min` here, and a negative interval propagates quietly:
    `30 // -26` is -2, which makes every run of mitotic frames long enough to
    count as an episode. Anything not finite and positive raises, naming the
    file, so the caller passes the interval in explicitly instead.

    Inputs:
    folder  : directory holding the metadata files, i.e. the one with the image
              stacks in it (`analysis.data_dir`)
    pattern : glob for the metadata files

    Returns the interval in minutes. Raises ValueError/FileNotFoundError with
    the offending file named.
    """
    folder = Path(folder)
    files = sorted(folder.glob(pattern))
    if not files:
        raise FileNotFoundError(
            f"no acquisition metadata matching {pattern!r} in {folder}; pass "
            f"frame_interval explicitly")

    found = {}
    for f in files:
        match = re.search(r"^Time interval:\s*(-?[\d.]+)\s*(\w*)",
                          f.read_text(errors="replace"), re.MULTILINE)
        if match is None:
            continue
        value, unit = float(match.group(1)), match.group(2).lower()
        if unit.startswith("s"):            # seconds, occasionally written
            value /= 60.0
        elif unit and not unit.startswith("min"):
            raise ValueError(f"{f.name} gives the time interval in {unit!r}, "
                             f"which is not a unit this understands; pass "
                             f"frame_interval explicitly")
        if not np.isfinite(value) or value <= 0:
            raise ValueError(
                f"{f.name} gives 'Time interval: {match.group(1)} "
                f"{match.group(2)}', which is not a usable frame interval. "
                f"The metadata is wrong; pass frame_interval explicitly.")
        found[f.name] = value

    if not found:
        raise ValueError(
            f"none of {[f.name for f in files]} in {folder} carries a "
            f"'Time interval:' line; pass frame_interval explicitly")
    if len(set(found.values())) > 1:
        raise ValueError(
            f"the metadata in {folder} disagrees about the frame interval "
            f"({found}); pass frame_interval explicitly")
    return next(iter(found.values()))
