# `signal_correction.py` — position-specific background correction

Measures the fluorescence background from each position's **own** frames instead
of from a blank well, and takes the illumination correction from the
**difference** of two blank wells. It replaces the `<ch>_bkg_corr` and
`<ch>_int_corr` columns that `cellaap_analysis` writes, without re-running
inference.

For why it is built this way, and what was measured to justify it, see
[SIGNAL_CORRECTION_DESIGN.md](SIGNAL_CORRECTION_DESIGN.md). For reading a saved
correction back, applying one to numbers already in a workbook, or lending a
shape to a position too crowded to measure its own, see
[CORRECTION_TOOLS_README.md](CORRECTION_TOOLS_README.md). For the analysis
module whose output it corrects, see [README.md](README.md).

---

## The five minutes that matter

```bash
conda activate img-env
```

```python
from pathlib import Path
import tifffile
from signal_correction import (flatfield_from_blank_pair,
                               estimate_position_correction)

ROOT = Path('/path/to/20353')
instance_stack = next(ROOT.glob('*_A08_s1_*_inference/*instance.tif'))

# 1. one flat field for the whole plate, from the two blank wells
flat = flatfield_from_blank_pair(
    tifffile.imread(next(ROOT.glob('*_D08_s8_GFP.tif')), key=range(0, 137, 16)),
    tifffile.imread(next(ROOT.glob('*_E08_s6_GFP.tif')), key=range(0, 137, 16)))

# 2. one background per position, from that position's own frames
with tifffile.TiffFile(next(ROOT.glob('*_A08_s1_GFP.tif'))) as tf:
    corr = estimate_position_correction(
        tf.series[0].asarray(out='memmap'),
        labels=tifffile.imread(instance_stack),
        stem='A08_s1', channel='GFP', flatfield=flat)

# 3. corrected per-cell signal, straight from the analysis table
df['GFP_corrected'] = corr.correct_measurements(
    df['GFP'], df['x'] * 2, df['y'] * 2, df['frame'])
```

Three things to know before running it:

* **`x` and `y` must be in fluorescence-frame pixels.** cellaap writes centroids
  at segmentation scale, which is half the fluorescence resolution — hence the
  `* 2`. Getting this wrong silently samples the wrong part of the field.
* **Pass a memory map**, not `tifffile.imread(...)`. Only ~16 frames are read,
  and the whole point is not to pull 1.1 GB off the share per position.
* **The flat field is per plate; the background is per position.** Do not
  compute a flat field per position — see design doc §3.

---

## 1. `flatfield_from_blank_pair` — the illumination, once per plate

```python
flat = flatfield_from_blank_pair(bright, dim, block=64, min_contrast=20.0)
```

Two blank wells imaged through the same optics are `D + A_bright·F` and
`D + A_dim·F`. Their **difference** is `(A_bright − A_dim)·F`: the camera offset
cancels exactly, so no one has to know it. Returns a mean-1 grid.

| argument | meaning |
| --- | --- |
| `bright`, `dim` | the two blank stacks, `(t, y, x)` or a single frame. On the 20250213 plate: DMEM is bright, FluoroBrite dim |
| `min_contrast` | refuse if the two differ by fewer than this many counts on average — below that the difference is mostly noise |

Order does not matter (the sign is fixed internally), but **a bigger gap between
the two wells is better**. Subsample frames on the way in (`key=range(0, n, 16)`);
a dozen is plenty.

**Sanity check:** the centre-to-edge ratio should agree between channels, since
it is a property of the optics. On the 20250213 plate GFP gives 1.220 and Texas
Red 1.228. If two channels disagree by much, something is wrong with the blanks
— see design doc §4.

Without a usable blank pair, run the background correction alone
(`flatfield=None`). It is the larger of the two corrections by an order of
magnitude on this data.

---

## 2. `estimate_position_correction` — the background, once per position

```python
corr = estimate_position_correction(
    stack, stem='', channel='', labels=None, flatfield=None, darkfield=None,
    background_model='auto', grid_min_blocks=0.05, n_frames=24, block=64,
    dilation=121, min_usable_blocks=0.05, how=None, smooth=1.5, keep_grids=True)
```

| argument | default | meaning |
| --- | --- | --- |
| `stack` | — | the movie; a tifffile memmap is fine and preferred |
| `labels` | `None` | instance or semantic segmentation, for masking cells out. **Strongly recommended** |
| `flatfield` | `None` | from step 1. Without it, background is corrected and illumination is not |
| `dilation` | `121` | how far to stay from cells, in fluorescence pixels. **The most important parameter** |
| `min_usable_blocks` | `0.05` | back off `dilation` on any frame where fewer than this fraction of blocks survive |
| `n_frames` | `24` | frames sampled for the surfaces; the level is interpolated between them |
| `block` | `64` | block side; 64 on a 2048 frame gives a 32×32 grid |
| `background_model` | `'auto'` | `'grid'` fills and smooths the measured grid; `'flatfield'` fits `offset + level·F` per frame; `'auto'` uses the grid and falls back when a position is too crowded |
| `grid_min_blocks` | `0.05` | block fraction below which `'auto'` switches to the flat-field fit |

### `dilation` is the one to think about

Out-of-focus light around a cell is not background, and subtracting it deletes
signal. On the 20250213 plate the medium 10 px from a cell reads 132 counts and
30 px out 127, against only ~13 counts of actual cell signal. Widening the
exclusion from 21 to 121 px raised recovered GFP from 12.0 to 15.9 counts and
cut the residual radial trend from 34% to 9%.

The cost is measurable blocks — 75% at 21 px, 10% at 121 px. The function backs
the exclusion off **per frame** (121 → 60 → 30 → 15 → none) until
`min_usable_blocks` survive, so a crowded frame still gets a background. Check
`diagnostics['dilation_used']` to see how far it had to give.

If your cells are larger, more spread, or imaged at higher magnification than
HT1080 at 20×, raise `dilation` in proportion.

### When a position is too crowded to measure a grid

`background_model='auto'` (the default) drops to a two-parameter fit,
`offset + level·F`, which needs ten measurable blocks rather than a whole grid.
It engages when fewer than `grid_min_blocks` of blocks survive, or when any
sampled frame has none at all — the case the grid model cannot survive. The
switch and its reason are logged, stored in `diagnostics`, and written to the
corrections sheet.

The threshold is deliberately low, because the fit is a **fallback, not an
improvement**: cross-validated on the 20250213 plate, the grid predicts
held-out background blocks better from 5% of blocks upwards (5.7 vs 5.8 counts
at 5%, 4.5 vs 5.7 at 20%) and only loses below ~4% (6.7 vs 6.0 at 2%). The
fit's error plateaus near 5.6 counts however many blocks it gets — that is what
two parameters buys. Force it with `background_model='flatfield'` if you want
it everywhere.

With no flat field *and* no measurable blocks there is nothing to fall back on,
and the estimator raises rather than inventing a background.

---

## 3. Using the result

`PositionCorrection` stores a 32×32 background shape, a 32×32 flat field, and
two numbers per frame — about 50 kB, against a full-size map stack.

```python
corr.correct_measurements(raw, x, y, frame)   # (raw − B) / F, per cell
corr.subtract_background(raw, x, y, frame)    # raw − B only
corr.apply(image, frame)                      # a whole frame
corr.background(frame)                        # B at frame resolution
corr.flat()                                   # F at frame resolution
corr.set_flatfield(flat)                      # attach one after the fact
corr.save(path)                               # .npz + a readable .json
PositionCorrection.load(path)
```

`correct_measurements` samples the fields at each cell's centroid rather than
averaging them over its mask, the way `<ch>_bkg_corr` does. Both fields are
smooth on a 64-px scale, far larger than a cell, so the two agree to well under
a percent.

### Applying it to existing `*_analysis.xlsx` files

The analysis tables already carry everything needed — a raw masked mean per cell
per frame, and the centroid it came from — so no images have to be re-read:

```python
cells = pd.read_excel(analysis_xlsx, sheet_name='cell_data')
cells['GFP_corrected'] = corr.correct_measurements(
    cells['GFP'], cells['x'] * 2, cells['y'] * 2, cells['frame'])
```

To reproduce a **summary-level** per-track number, average the corrected
per-frame values over that track's mitotic window
(`mitotic_start_frame` … `+ corrected_frames_in_mitosis / DELTA_T`). That
reconstruction reproduces the summary's own raw column to r = 0.995 and −0.5%
median, so the corrected numbers are directly comparable to the stock ones.

---

## 3b. Inspecting the saved background surfaces

The pipeline writes the background it actually subtracted, per position and
channel, to `<root>/pipeline/state/surfaces/<stem>_<channel>_bkg.tif`:
`(n_frames, 32, 32)` float32 **in counts**, ~490 kB per position-channel
compressed — about 12 MB for a 13-position two-channel plate, against ~42 GB of
raw stacks. The `corrections` sheet of each summary names the file.

```python
from signal_correction import read_background_stack, upsample_stack

stack, meta = read_background_stack(path)              # (137, 32, 32), counts
stack, meta = read_background_stack(path, shape=(2048, 2048))   # frame-sized
one_frame   = upsample_stack(stack[68], (2048, 2048))  # just one, cheaply
```

`meta` carries the stem, channel, block size, frame shape and the full
diagnostics dict, so a stack found on disk months later explains itself.

It is stored at grid resolution because the surface has no structure finer than
one block — upsampling adds pixels, not information. Expanding a whole 137-frame
movie to 2048² is 2.3 GB, so take a slice unless you mean it.

**The file name must not contain `background` or `intensity`.**
`cellaap_analysis._load_maps` walks the entire root — the pipeline directory
included — and loads any file matching those words as a *plate-wide* correction
map. A per-position surface caught that way would be applied to every position,
silently reinstating the blank-well bug. `save_background_stack` refuses such a
name, and the pipeline names them `_bkg.tif`.

---

## 4. The building blocks

`estimate_position_correction` is a driver over these; reach for them directly
only when you want to do something it does not.

| function | |
| --- | --- |
| `cell_free_mask_from_labels(labels, shape, dilation)` | boolean mask, True where no cell is. Upsamples the segmentation to the fluorescence frame and dilates by a square |
| `block_reduce_robust(image, block, how, mask=...)` | image → grid of per-block statistics. `how` is `median`, `mean`, `std`, `quantile` or `clipped_mean`; NaN where too few pixels survive |
| `background_surface(image, cell_free, ...)` | one frame's background: block-reduce, fill, smooth. `full_resolution=True` returns it frame-sized |
| `background_surfaces(frames, cell_free, ...)` | the same over a sequence → `(n, gy, gx)` |
| `fit_background_to_flatfield(grid, flatfield)` | fits `offset + level·F` to one grid, returns `(offset, level)`. The `background_model='flatfield'` path |
| `upsample(grid, shape)` | bilinear resize of a coarse grid to a full frame |
| `upsample_stack(stack, shape)` | the same, frame by frame over an `(n, gy, gx)` stack |
| `save_background_stack(corr, path)` | write the per-frame background as a TIFF, with metadata |
| `read_background_stack(path, shape=None)` | read one back, optionally frame-sized |
| `centre_edge_ratio(grid)` | middle-to-rim ratio; a one-number vignette depth |

Note `block_reduce_robust`'s `clipped_mean`: it clips the **high side only**,
because cells only ever add signal, and it takes its width from a MAD measured
on the low half of the block so cells cannot inflate the threshold they are
about to be tested against. That is what makes it usable without a mask. The
`quantile` reduction is available but biased low by 12–25 counts on this data —
use it only where a constant offset does not matter.

---

## 5. What to check afterwards

Every result carries a `diagnostics` dict, also written to the `.json` sidecar:

| key | what a healthy value looks like |
| --- | --- |
| `background_mean_range` | the position's background in counts, first to last frame |
| `background_drift_percent` | 8–15% on this plate. Near 0% with a real drift means the level was not tracked |
| `dilation_used` | `[min, max]` actually applied. Far below `dilation_requested` means a crowded field |
| `unusable_block_fraction` | 0.7–0.9 is normal at 121 px. 1.0 means nothing was measurable |
| `background_shape_centre_edge` | ~1.06–1.09. Much steeper suggests cells leaking into the background |
| `fitted_offset_range` | only under `background_model='flatfield'`; should be stable across frames |

Then check the correction actually flattened the field, using the same metric as
the aggregate notebook:

```python
from signal_correction import radial_profile, flatness
_, prof = radial_profile(corrected, x, y, shape=(2048, 2048))
print(flatness(prof))       # rms_deviation_percent is the verdict
```

**Judge flatness in counts, not percent.** A stage whose median is 16 counts and
one whose median is 140 cannot be compared on a percentage — multiply by the
median first. See design doc §8.

---

## 6. Cost

| | |
| --- | --- |
| per position per channel | ~14 s, nearly all of it reading frames off the share |
| CPU proper | ~1 s |
| memory | one frame at a time (8 MB) plus 32×32 grids |
| stored result | ~50 kB |

Reads only `n_frames` of the movie. The expensive-looking mask dilation uses
`ndi.maximum_filter`, which is the separable form of dilating by a square and
gives a bit-identical result ~100× faster (0.03 s against 4.9 s for a 121-px
square on a 2048 frame) — that is what makes a per-frame dilation ladder
affordable.

---

## 7. Limits

* **It needs a segmentation to work well.** Without `labels` it falls back to a
  sigma-clipped block mean, which is decent but measurably worse, and it cannot
  keep its distance from cells.
* **It does not help with faint cells.** On the 20250213 plate cells sit ~10%
  above their own local background (140 counts inside the mask against 127 in
  the gaps). At that contrast a 1-count error in the background is a ~10% error
  in the cell, and ~5% of Texas Red tracks come out slightly negative — real
  cells at their local background, previously hidden behind the stock
  correction's under-subtraction. No estimator gets under that floor; a brighter
  reporter or a longer exposure would.
* **It corrects a plate, not an experiment.** The flat field is valid for one
  acquisition session. Re-derive it if the objective, filter or illumination
  path changed.
* **It does not touch `cellaap_analysis` or `pipeline.py`.** Corrections are
  computed and applied downstream; the analysis stage still writes its own
  `<ch>_bkg_corr` / `<ch>_int_corr` columns from whatever maps are on disk.
