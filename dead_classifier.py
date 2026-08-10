'''
Mitotic vs dead discrimination for cells carrying the mitotic semantic label.

ONLY mitotic-labeled detections may be scored. The model is binary within the
rounded-cell population - class 1 "mitotic" against class 0 "dead" - so on a
flat interphase cell class 0 means nothing more than "not rounded". Measured
on 20260624_CycB oe_A12_s2: frames BEFORE mitotic entry, from cells that go on
to divide and are therefore alive by construction, score mean P(dead) = 0.785
with 82.5% over threshold, against 0.869/90.2% for post-mitotic frames and
0.294/27.4% for mitotic ones. Extending scoring past the mitotic label to
catch post-mitotic death was tried and reverted for exactly this reason: it
moved dead_post_mitosis 1 -> 37 and mitotic_survived 114 -> 27, almost all
artifact. Post-mitotic death is instead recovered in summarize_data, where a
dying cell that re-rounds carries the mitotic label again and so stays inside
this model's competence - see summarize_multi_peak_tracks.

Each detection is reduced to an instance-masked phase crop, from which two
feature blocks are computed and concatenated in this order:

    [ 512 ImageNet ResNet18 embedding ] + [ 16 handcrafted texture/shape ]

A model bundle maps that to P(mitotic). Class 1 = mitotic, class 0 = dead.
Which of the two blocks the model actually consumes is declared by the
bundle's 'feature_set' key and assembled by _feature_matrix() - the current
pooled model uses the embedding alone (512 columns), while the previous
hand-labeled model used both (528). Passing the wrong width silently
misclassifies everything, so the width is checked before predicting.

Current model, models/dead_classifier_pooled.joblib: ResNet18-512 -> PCA(32)
-> logistic regression, trained on 826 hand labels pooled from the BUB1
dataset (345) and the CycB-oe 20576 dataset (481) across 15 microscope
positions. Leave-one-position-out balanced accuracy 0.878, AUC 0.945,
calibration error 0.018.

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

# Semantic values that mean "mitotic". Different cellaap versions write
# different ones (101 in the current pipeline, 100 in others); a detection is
# classified mitotic-vs-dead only if its semantic label is one of these.
MITOTIC_SEMANTIC_VALUES = (100, 101)

# Feature layout of the pre-pooled model, used when a bundle does not say.
DEFAULT_FEATURE_SET = 'resnet18 + handcrafted'

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


def unpack_model(model):
    '''Accept either a loaded bundle dict or a bare estimator.

    Returns (estimator, feature_set). Bare estimators predate the
    'feature_set' key and always used the 528-column layout.
    '''
    if isinstance(model, dict):
        return model['model'], model.get('feature_set', DEFAULT_FEATURE_SET)
    return model, DEFAULT_FEATURE_SET


'''Which feature blocks each known layout consumes: (embedding, handcrafted).'''
_BLOCKS = {
    'resnet18 + handcrafted': (True, True),
    'resnet18 (512-d)':       (True, False),
    'handcrafted, all 16':    (False, True),
}


def feature_blocks(feature_set: str) -> tuple:
    '''(needs_embedding, needs_handcrafted) for a declared feature set.

    Callers use this to skip work the model will not look at. The handcrafted
    block costs a GLCM per crop, which dominates runtime over a whole stack,
    and the current model does not use it.
    '''
    try:
        return _BLOCKS[feature_set]
    except KeyError:
        raise ValueError(
            f"unknown feature_set {feature_set!r}; the model bundle must "
            f"declare one of {sorted(_BLOCKS)}") from None


def _feature_matrix(masked_imgs: npt.NDArray, handcrafted,
                    feature_set: str) -> npt.NDArray:
    '''Assemble the design matrix the trained model expects.

    handcrafted may be None when feature_blocks() says it is not needed.
    '''
    use_emb, use_hand = feature_blocks(feature_set)
    if use_hand and handcrafted is None:
        raise ValueError(f"feature_set {feature_set!r} needs the handcrafted "
                         f"block but none was computed")
    if use_emb and use_hand:
        return np.column_stack([_embed(masked_imgs), handcrafted])
    if use_emb:
        return _embed(masked_imgs)
    return np.asarray(handcrafted)


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


def mitotic_rows(tracking_df: pd.DataFrame,
                 mitotic_values=MITOTIC_SEMANTIC_VALUES) -> pd.DataFrame:
    '''Rows whose semantic label marks them mitotic, and so worth classifying.

    Gating on the raw semantic label rather than on semantic_smoothed is
    deliberate: smoothing closes gaps and median-filters the trace, so it both
    adds frames the segmentation never called mitotic and drops frames it did.
    The probabilities should describe what the segmentation actually saw.
    Falls back to semantic_smoothed only for tables that predate the column.
    '''
    if 'semantic' in tracking_df.columns:
        return tracking_df[tracking_df['semantic'].isin(mitotic_values)]
    if 'semantic_smoothed' in tracking_df.columns:
        return tracking_df[tracking_df['semantic_smoothed'] == 1]
    raise KeyError("tracking_df needs a 'semantic' (or 'semantic_smoothed') column")


def rows_to_classify(tracking_df: pd.DataFrame,
                     mitotic_values=MITOTIC_SEMANTIC_VALUES,
                     post_peak_frames: int = 0) -> pd.DataFrame:
    '''Mitotic-labeled rows, plus post_peak_frames frames after each episode.

    The tail is what makes post-mitotic death visible: a cell dying as it
    leaves mitosis drops the mitotic label while still rounded, so the frames
    carrying the death sit just past the episode - close enough that the cell
    still looks like something this model was trained on.

    Keep the tail short. Scored far enough past the peak the cell is flat
    again, and a flat cell reads as "dead" to this model whether or not it is:
    on A12_s2, alive pre-mitotic frames scored mean P(dead) = 0.785. The tail
    length is analysis_pars.post_peak_frames.

    Episodes are runs of the raw mitotic label, for the same reason
    mitotic_rows gates on it rather than on semantic_smoothed. Tables with no
    particle/frame columns fall back to the mitotic rows alone.
    '''
    mitotic = mitotic_rows(tracking_df, mitotic_values)
    if post_peak_frames <= 0 or not {'particle', 'frame'}.issubset(tracking_df.columns):
        return mitotic

    keep = pd.Series(False, index=tracking_df.index)
    keep.loc[mitotic.index] = True
    is_mitotic = keep.to_numpy().copy()
    for _, g in tracking_df.groupby('particle', sort=False):
        rows = tracking_df.index.get_indexer(g.index)
        gm = is_mitotic[rows]
        if not gm.any():
            continue
        frames = g['frame'].to_numpy()
        # ends of each run of mitotic frames within this track
        padded = np.concatenate(([False], gm, [False]))
        edges = np.flatnonzero(padded[1:] != padded[:-1])
        for start, stop in zip(edges[::2], edges[1::2]):
            last = frames[stop - 1]
            tail = (frames > last) & (frames <= last + post_peak_frames)
            if tail.any():
                keep.iloc[rows[tail]] = True
    return tracking_df[keep]


def classify_dead(phase_stack: npt.NDArray, instance_stack: npt.NDArray,
                  tracking_df: pd.DataFrame, model,
                  phase_offset: int = 0,
                  mitotic_values=MITOTIC_SEMANTIC_VALUES,
                  post_peak_frames: int = 0,
                  chunk_size: int = 2048) -> pd.DataFrame:
    '''
    Classify the mitotic-labeled detections, plus a short post-episode tail.

    Do not widen this beyond the tail - see the module docstring. Far from the
    rounded-cell population the model's class 0 means "not mitotic", not "dead".

    Inputs:
    phase_stack    : full-resolution phase stack (T, 2048, 2048)
    instance_stack : instance segmentation at analysis scale (T, 1024, 1024)
    tracking_df    : needs semantic, frame, x, y, label, area
                     (x = row, y = col, analysis scale)
    model          : model bundle dict, or a bare fitted classifier with
                     predict_proba; class 1 = mitotic
    phase_offset   : analysis frame f corresponds to phase page f + phase_offset.
                     Some acquisitions keep leading phase frames that were not
                     segmented (e.g. a 361-frame phase stack against a 341-frame
                     segmentation needs phase_offset=20).
    mitotic_values : semantic values treated as mitotic (default 100 and 101).
    post_peak_frames : frames scored after each mitotic episode, so a death on
                     mitotic exit is still seen. See rows_to_classify.
    chunk_size     : crops buffered before a predict pass, to cap peak memory.
                     A long movie yields tens of thousands of mitotic crops,
                     and 96x96 float32 adds up (~37 MB per 1000).

    Returns a dataframe indexed like tracking_df, holding only the scored
    rows that could be cropped, with columns
    mitotic_proba, dead_proba (= 1 - mitotic_proba) and dead_flag (1 = dead).
    Rows outside that set are left unclassified by the caller rather than
    being assigned a probability the model was never asked for.
    '''
    estimator, feature_set = unpack_model(model)

    # A frame-count mismatch silently pulls crops from the wrong time point,
    # which is invisible in the output, so refuse rather than guess.
    n_phase, n_inst = len(phase_stack), len(instance_stack)
    if n_phase - phase_offset != n_inst:
        raise ValueError(
            f"phase stack has {n_phase} frames, instance stack has {n_inst}, "
            f"and phase_offset={phase_offset} does not reconcile them. "
            f"Analysis frame f must map to phase page f + phase_offset; pass "
            f"phase_offset={n_phase - n_inst} if the extra phase frames are "
            f"leading frames that were not segmented.")

    rows = rows_to_classify(tracking_df, mitotic_values, post_peak_frames)
    _, needs_handcrafted = feature_blocks(feature_set)

    # Scored in chunks rather than all at once. Holding every crop as 96x96
    # float32 costs ~37 MB per 1000, and np.stack then doubles it; a 450-frame
    # movie has tens of thousands of mitotic detections. A chunk is ~75 MB,
    # and only the probabilities are kept.
    kept, probas = [], []
    imgs, feats, idx = [], [], []

    def _score_chunk():
        '''Embed and predict the buffered crops, then release them.'''
        if not imgs:
            return
        X = _feature_matrix(np.stack(imgs),
                            np.stack(feats) if needs_handcrafted else None,
                            feature_set)
        # A width mismatch produces confident nonsense rather than an obvious
        # failure, so check it against what the estimator was fitted on.
        expected = getattr(estimator, 'n_features_in_', None)
        if expected is not None and X.shape[1] != expected:
            raise ValueError(
                f"model expects {expected} features but feature_set "
                f"{feature_set!r} produced {X.shape[1]}. The bundle's feature_set "
                f"does not match the estimator it was saved with.")
        probas.append(estimator.predict_proba(X)[:, 1])
        kept.extend(idx)
        imgs.clear(), feats.clear(), idx.clear()

    for index, row in rows.iterrows():
        frame = int(row['frame'])
        if frame < 0 or frame >= n_inst:
            continue
        got = prepare_crop(phase_stack[frame + phase_offset], instance_stack[frame],
                           int(row['x']), int(row['y']), row['label'])
        if got is None:
            continue
        im, mask = got
        imgs.append(im)
        # the GLCM in _crop_features dominates runtime over a whole stack;
        # skip it entirely when the model does not consume that block
        if needs_handcrafted:
            feats.append(np.append(_crop_features(im, mask), row['area']))
        idx.append(index)
        if len(imgs) >= chunk_size:
            _score_chunk()
    _score_chunk()

    if not kept:
        return pd.DataFrame(columns=["mitotic_proba", "dead_proba", "dead_flag"])

    p_mitotic = np.concatenate(probas)
    # The two probabilities are complementary by construction (binary model);
    # both are reported so downstream code never has to remember which way
    # round the class encoding runs.
    return pd.DataFrame({"mitotic_proba": p_mitotic,
                         "dead_proba": 1.0 - p_mitotic,
                         "dead_flag": (p_mitotic < 0.5).astype(int)},
                        index=kept)
