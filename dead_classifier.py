'''
Mitotic vs dead discrimination for cells carrying the mitotic semantic label.

Each detection is reduced to an instance-masked phase crop, from which two
feature blocks are computed and concatenated in this order:

    [ 512 ImageNet ResNet18 embedding ] + [ 16 handcrafted texture/shape ]

A Platt-calibrated random forest (models/dead_classifier_handlabeled.joblib) maps that
to P(mitotic). Class 1 = mitotic, class 0 = dead.

Trained on 345 hand-labeled crops spanning 8 microscope positions;
leave-one-stack-out balanced accuracy 0.873, AUC 0.923, calibration error 0.050.

An earlier version of this model was trained on unsupervised pseudo-labels and
collapsed onto an area threshold (~639 px), which systematically called small
but healthy mitotic cells dead. Hand labels replaced those pseudo-labels.
Feature order and preprocessing must not change - they match the trained model.
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

CROP_SIZE = 96          # full-resolution crop edge; must match training
_DILATE = np.ones((3, 3), dtype=bool)   # matches the training-time dilation
_net = None


def _resnet():
    '''Lazily build the frozen ImageNet ResNet18 used as a feature extractor.'''
    global _net
    if _net is None:
        import torch
        from torchvision.models import resnet18, ResNet18_Weights
        net = resnet18(weights=ResNet18_Weights.IMAGENET1K_V1)
        net.fc = torch.nn.Identity()
        dev = 'mps' if torch.backends.mps.is_available() else 'cpu'
        _net = (net.eval().to(dev), dev)
    return _net


def _embed(masked_imgs: npt.NDArray) -> npt.NDArray:
    '''ResNet18 penultimate features for a stack of masked [0,1] crops.'''
    import torch
    net, dev = _resnet()
    out = []
    masked_imgs = np.asarray(masked_imgs, dtype=np.float32)
    with torch.no_grad():
        for i in range(0, len(masked_imgs), 256):
            b = torch.from_numpy(masked_imgs[i:i + 256]).unsqueeze(1)
            b = ((b - 0.449) / 0.226).repeat(1, 3, 1, 1).to(dev)
            out.append(net(b).cpu().numpy())
    return np.concatenate(out)


def _crop_features(im: npt.NDArray, mask: npt.NDArray) -> npt.NDArray:
    '''15 texture/shape features; the 16th (area) is appended by the caller.'''
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


def prepare_crop(phase_frame: npt.NDArray, instance_frame: npt.NDArray,
                 ar: int, ac: int, label: int):
    '''
    Build the masked, normalized crop for one detection.

    ar, ac are the analysis-scale centroid (row, col); the phase frame is at
    twice that resolution. Returns (masked_image, mask) or None if the crop
    would fall outside the frame or the mask is too small.
    '''
    half = CROP_SIZE // 2
    H, W = phase_frame.shape
    r0, c0 = ar * 2 - half, ac * 2 - half
    if r0 < 0 or c0 < 0 or r0 + CROP_SIZE > H or c0 + CROP_SIZE > W:
        return None

    window = instance_frame[max(0, ar - 24):ar + 24,
                            max(0, ac - 24):ac + 24] == label
    canvas = np.zeros((48, 48), dtype=bool)
    rr = 0 if ar - 24 >= 0 else 24 - ar
    cc = 0 if ac - 24 >= 0 else 24 - ac
    canvas[rr:rr + window.shape[0], cc:cc + window.shape[1]] = window
    mask = ndi.binary_dilation(np.kron(canvas, np.ones((2, 2), dtype=bool)),
                               structure=_DILATE)
    if mask.sum() < 50:
        return None

    crop = phase_frame[r0:r0 + CROP_SIZE, c0:c0 + CROP_SIZE].astype(np.float32)
    # np.percentile returns float64; keep everything float32 to match training
    # and because torch on MPS rejects float64.
    lo, hi = np.percentile(crop, [1, 99.5]).astype(np.float32)
    im = np.clip((crop - lo) / (hi - lo + np.float32(1e-6)), 0, 1)
    return np.where(mask, im, im[mask].mean()).astype(np.float32), mask


def classify_dead(phase_stack: npt.NDArray, instance_stack: npt.NDArray,
                  tracking_df: pd.DataFrame, model) -> pd.DataFrame:
    '''
    Classify every detection labeled mitotic by the semantic segmentation.

    Inputs:
    phase_stack    : full-resolution phase stack (T, 2048, 2048)
    instance_stack : instance segmentation at analysis scale (T, 1024, 1024)
    tracking_df    : needs semantic_smoothed, frame, x, y, label, area
                     (x = row, y = col, analysis scale)
    model          : fitted classifier with predict_proba; class 1 = mitotic

    Returns a dataframe indexed like tracking_df (only the mitotic-labeled rows
    that could be cropped) with dead_flag (1 = dead) and dead_proba.
    '''
    T = len(phase_stack)
    rows = tracking_df[tracking_df["semantic_smoothed"] == 1]

    imgs, feats, kept = [], [], []
    for index, row in rows.iterrows():
        frame = int(row['frame'])
        if frame >= T:
            continue
        got = prepare_crop(phase_stack[frame], instance_stack[frame],
                           int(row['x']), int(row['y']), row['label'])
        if got is None:
            continue
        im, mask = got
        imgs.append(im)
        feats.append(np.append(_crop_features(im, mask), row['area']))
        kept.append(index)

    if not feats:
        return pd.DataFrame(columns=["dead_flag", "dead_proba"])

    X = np.column_stack([_embed(np.stack(imgs)), np.stack(feats)])
    p_mitotic = model.predict_proba(X)[:, 1]
    return pd.DataFrame({"dead_flag": (p_mitotic < 0.5).astype(int),
                         "dead_proba": 1.0 - p_mitotic},
                        index=kept)
