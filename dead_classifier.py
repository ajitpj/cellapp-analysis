'''
Mitotic vs dead-like discrimination for cells carrying the mitotic semantic label.

A HistGradientBoosting classifier (models/dead_classifier_hgb.joblib) operates on
16 handcrafted texture/shape features computed from an instance-masked phase
crop of each detection. Class 1 = live mitotic (large, smooth, round), class
0 = dead-like (small, irregular, granular). Trained on 8,495 track-deduplicated
crops from the 2025-07-10 BUB1 experiment; leave-one-stack-out accuracy 0.994.

Feature order must not change - it matches the trained model.
'''
import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.ndimage as ndi
from skimage.feature import graycomatrix, graycoprops
from skimage.filters import sobel, laplace
from skimage.measure import shannon_entropy, regionprops

FEATURE_NAMES = ['std', 'entropy', 'glcm_contrast', 'glcm_homog', 'glcm_corr',
                 'glcm_ASM', 'grad_mean', 'grad_std', 'lap_mean', 'iqr80',
                 'rim_minus_inner', 'solidity', 'ecc', 'roughness',
                 'granule_density', 'area']

CROP_SIZE = 96  # full-resolution crop edge; must match classifier training


def _crop_features(im: npt.NDArray, mask: npt.NDArray) -> npt.NDArray:
    '''
    15 texture/shape features from a normalized [0,1] phase crop and its
    boolean particle mask (16th feature, area, is appended by the caller).
    '''
    inner = im[mask]
    q = (im * 31).astype(np.uint8)
    glcm = graycomatrix(q, [2, 4], [0, np.pi / 2], levels=32,
                        symmetric=True, normed=True)
    grad = sobel(im)
    lap = laplace(im)
    rp = regionprops(mask.astype(np.uint8))[0]
    rim = ndi.binary_dilation(mask) ^ ndi.binary_erosion(mask)
    return np.array([
        inner.std(),
        shannon_entropy((im * 255).astype(np.uint8)[mask]),
        graycoprops(glcm, 'contrast').mean(),
        graycoprops(glcm, 'homogeneity').mean(),
        graycoprops(glcm, 'correlation').mean(),
        graycoprops(glcm, 'ASM').mean(),
        grad[mask].mean(), grad[mask].std(),
        np.abs(lap)[mask].mean(),
        np.percentile(inner, 90) - np.percentile(inner, 10),
        im[rim].mean() - inner.mean(),
        rp.solidity, rp.eccentricity,
        rp.perimeter / (2 * np.sqrt(np.pi * rp.area)),
        (ndi.minimum_filter(im, 5) == im)[mask].sum() / mask.sum(),
    ], dtype=np.float32)


def classify_dead(phase_stack: npt.NDArray, instance_stack: npt.NDArray,
                  tracking_df: pd.DataFrame, model) -> pd.DataFrame:
    '''
    Classify every detection labeled mitotic by the semantic segmentation as
    live mitotic or dead-like.

    Inputs:
    phase_stack    : full-resolution phase stack (T, 2048, 2048)
    instance_stack : instance segmentation at analysis scale (T, 1024, 1024)
    tracking_df    : tracking dataframe; needs semantic_smoothed, frame, x, y,
                     label, area columns (x = row, y = col, analysis scale)
    model          : fitted classifier with predict_proba

    Returns a dataframe indexed like tracking_df (mitotic rows that could be
    cropped only) with columns dead_flag (1 = dead-like) and dead_proba
    (probability of the dead-like class).
    '''
    half = CROP_SIZE // 2
    T, H, W = phase_stack.shape
    rows = tracking_df[tracking_df["semantic_smoothed"] == 1]

    feats, kept_index = [], []
    for index, row in rows.iterrows():
        ar, ac = int(row['x']), int(row['y'])       # analysis scale (row, col)
        r0, c0 = ar * 2 - half, ac * 2 - half
        frame = int(row['frame'])
        if r0 < 0 or c0 < 0 or r0 + CROP_SIZE > H or c0 + CROP_SIZE > W or frame >= T:
            continue

        # particle mask: 48x48 analysis-scale window, upsampled x2, dilated
        window = instance_stack[frame,
                                max(0, ar - 24):ar + 24,
                                max(0, ac - 24):ac + 24] == row['label']
        canvas = np.zeros((48, 48), dtype=bool)
        rr = 0 if ar - 24 >= 0 else 24 - ar
        cc = 0 if ac - 24 >= 0 else 24 - ac
        canvas[rr:rr + window.shape[0], cc:cc + window.shape[1]] = window
        mask = ndi.binary_dilation(np.kron(canvas, np.ones((2, 2), dtype=bool)))
        if mask.sum() < 50:
            continue

        crop = phase_stack[frame, r0:r0 + CROP_SIZE, c0:c0 + CROP_SIZE].astype(np.float32)
        lo, hi = np.percentile(crop, [1, 99.5])
        im = np.clip((crop - lo) / (hi - lo + 1e-6), 0, 1)
        im = np.where(mask, im, im[mask].mean())

        feats.append(np.append(_crop_features(im, mask), row['area']))
        kept_index.append(index)

    if not feats:
        return pd.DataFrame(columns=["dead_flag", "dead_proba"])

    p_live = model.predict_proba(np.stack(feats))[:, 1]
    return pd.DataFrame({"dead_flag": (p_live < 0.5).astype(int),
                         "dead_proba": 1.0 - p_live},
                        index=kept_index)
