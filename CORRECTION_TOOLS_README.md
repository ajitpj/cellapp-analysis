# `correction_tools.py` — using a correction after the pipeline has run

`signal_correction` estimates a correction *while* a position is analyzed.
This module is everything you might want to do with one **afterwards**: look at
the background that was subtracted, apply a correction to numbers already in a
workbook, or repair a position that was too crowded to measure its own.

The two are kept apart so neither clutters the other. Nothing you need during a
pipeline run is in here, and nothing in here is needed to run the pipeline.

* the estimator: [SIGNAL_CORRECTION_README.md](SIGNAL_CORRECTION_README.md)
* why it works that way: [SIGNAL_CORRECTION_DESIGN.md](SIGNAL_CORRECTION_DESIGN.md)
* running a plate: [PIPELINE_README.md](PIPELINE_README.md)

---

## What the pipeline leaves you

```
<root>/pipeline/state/flatfield/<channel>.npz            one per plate
<root>/pipeline/state/surfaces/<stem>_<channel>_bkg.tif  one per position
```

The surface is `(n_frames, 32, 32)` float32 **in counts** — the background that
was actually subtracted, at the correction's own grid resolution. About 490 kB
per position-channel, ~12 MB for a 13-position two-channel plate. Together with
the flat field it reconstructs the whole correction, so nothing here re-reads an
image or re-estimates anything.

The `corrections` sheet of each `*_summary.xlsx` names the surface file it used.

---

## 1. Look at a background surface

```python
from pathlib import Path
from correction_tools import read_background_stack, upsample_stack

root = Path('/path/to/20353')
p = root / 'pipeline/state/surfaces/20250213_HT1080 pPS18_A08_s1_phs_GFP_bkg.tif'

stack, meta = read_background_stack(p)
print(stack.shape, stack.min(), stack.max())      # (137, 32, 32) 114.0 149.7
print(meta['diagnostics']['background_model'])    # 'grid'
print(meta['diagnostics']['dilation_used'])       # [60, 121]
```

`meta` carries the stem, channel, block size, frame shape and the full
diagnostics dict, so a stack found on disk months later explains itself.

To overlay it on the raw images, expand it:

```python
frame_68 = upsample_stack(stack[68], tuple(meta['frame_shape']))   # (2048, 2048)
whole    = read_background_stack(p, shape=(2048, 2048))[0]         # 2.3 GB — careful
```

It is stored at grid resolution because the surface has no structure finer than
one block; upsampling adds pixels, not information. Take a slice unless you
really want the whole movie at full size.

The file is opened by Fiji directly — it is an ordinary float32 TIFF stack.

### The whole plate at once

`surface_diagnostics` reads only the TIFF headers, so it costs a fraction of a
second for a plate and no image data is touched:

```python
from correction_tools import surface_diagnostics

d = surface_diagnostics(root, channel='GFP')
d[['well', 'position', 'usable_blocks', 'dilation_min', 'drift_percent']]
```

`usable_blocks` is the fraction of the grid that had enough cell-free pixels to
measure; `dilation_min`/`dilation_max` are how far from the cells the estimator
managed to stay over the sampled frames; `drift_percent` is how much the mean
background moved over the movie. They are the three symptoms of a crowded
field, and they move together: as the field fills, the exclusion ring backs
off, more out-of-focus halo is counted as medium, the background rises through
the movie — and the corrected signal comes out too low.

This is the companion to `cellaap_aggregate.correct_wells`, which shifts each
well by the over-subtraction read off its negative cells; this says *why* a
field was over-subtracted. On the 20260826 CycB plate, over the ten positions
of the two wells with no GFP induced, the bottom of each position's
distribution tracked all three: r = −0.83 against peak background, −0.80
against drift, +0.79 against the narrowest dilation held. That is also why
`correct_wells` works per well and not per position: the error is real at the
position level, but a few hundred cells per position do not measure it well
enough to subtract.

Read it only across positions expected to hold the same fluorophore. A well
that is genuinely brighter has a genuinely higher floor, and this table cannot
tell you which you are looking at — it tells you whether the estimator was in
trouble.

---

## 2. Apply a correction to numbers already measured

```python
import pandas as pd
from correction_tools import Surface, correct_cell_table

surface = Surface.from_files(
    root / 'pipeline/state/surfaces/..._A08_s1_phs_GFP_bkg.tif',
    root / 'pipeline/state/flatfield/GFP.npz')          # omit for background only

cells = pd.read_excel(analysis_xlsx, sheet_name='cell_data')
cells = correct_cell_table(cells, {'GFP': surface})     # adds GFP_corrected
```

That is the whole operation: `(raw − B) / F`, sampled at each cell's centroid.
No images are read.

**Centroid scale.** cellaap writes centroids at segmentation scale, half the
fluorescence resolution, so they have to be doubled before the correction is
sampled. `correct_cell_table` infers the factor from the frame shape and checks
it, raising if the centroids fall outside the frame and warning if they only
reach one corner of it — getting this wrong samples the wrong part of the field
and produces numbers that look plausible. Pass `xy_scale=2.0` explicitly if you
would rather be sure.

### Carrying it to the summary

The summary's per-track signal is the mean over each track's mitotic window, so
a new per-frame column has to be averaged the same way:

```python
from correction_tools import track_means

summary = pd.read_excel(summary_xlsx, sheet_name='Summary')
summary = track_means(cells, summary, ['GFP_corrected'], delta_t=10)
```

`delta_t` converts `corrected_frames_in_mitosis` back to frames. It holds
**minutes** once `cellaap_aggregate.load_experiment` has been through the
summary and **frames** before that, so a summary read straight off disk needs
`delta_t=1`. Getting it backwards shortens every window and still returns
plausible numbers — on one position, r falls from 0.9996 to 0.970 with a −4%
bias — so `track_means` warns when the durations look like frames but you asked
for minutes.

**Check it first.** Ask for the raw column and compare against the summary's
own; a column the summary already has comes back suffixed rather than replacing
it:

```python
check = track_means(cells, summary, ['GFP'], delta_t=10)
# check['GFP'] is the summary's; check['GFP_recomputed'] is yours
```

On the 20250213 plate that reproduces the summary to **r = 0.999, −0.6%**. If
your plate does not, the window reconstruction is wrong and the corrected
numbers are not comparable to the stock ones.

---

## 3. Repair a position that could not measure its own background

For a field so crowded that almost no block has cell-free pixels. Borrow the
**shape** from other positions and keep the position's **own level**:

```python
from correction_tools import borrow_shape

fixed = borrow_shape(crowded, donors=[other, another], note='98% confluent')
```

**Only the shape is borrowed, and that asymmetry is the point.**

| | varies between positions | safe to borrow? |
| --- | --- | --- |
| shape | 1.6% (GFP) / 2.1% (TR) from the plate mean ≈ 2–3 counts | **yes** |
| level | 4.1 (GFP) / 7.3 (TR) counts | **no** |

The level is not merely noisier — it **correlates with confluence**
(r = +0.51, ~26 counts per unit confluence). A position that cannot measure its
own background is by definition the most crowded one, so its true level is the
highest on the plate while every donor is sparser and therefore darker.
Borrowing the level would under-subtract by ~5 counts against a ~13-count
signal, systematically, and worst in exactly the wells where cells grew densest
— which would read as a biological effect.

So the level is re-derived locally: the position's own measured block grids are
projected onto the borrowed shape, the same least-squares step the estimator
uses. On a real position this moves the mean background by **−0.01%** while the
shape changes from centre/edge 1.071 to 1.053 — the level stays at that
position's own 131.7 counts rather than being pulled toward the donors' 123.0.

Related notes:

* Donors must share the same `block`, so their grids are the same size.
* `average_shape(donors)` gives you the mean shape on its own.
* The result records `shape_borrowed_from`, `level_kept_local` and the method
  in its `diagnostics`, so a summary built from it says where it came from.
* Where `background_grids` were not kept (`keep_grids=False`), the level is
  instead rescaled to preserve the mean background — weaker, and flagged.
* `borrow_shape` returns a **new** correction; the input is untouched.

**Which donors?** Any positions that measured cleanly. Do not bother preferring
the same well: on the 20250213 plate the within-well spread of the level is
4.0 counts against 0.9 counts between well means, so essentially all the
variation is *within* wells and a same-well hierarchy buys nothing. Pool
whatever is available — more donors, less shape noise.

---

## 4. Writing a surface yourself

```python
from correction_tools import save_background_stack
save_background_stack(correction, out_dir / 'A08_s1_GFP_bkg.tif')
```

**Never name one `*_background.tif` or `*_intensity.tif`.**
`cellaap_analysis._load_maps` walks the *entire* root — the pipeline directory
included — and adopts any file matching those words as a **plate-wide**
correction map, applying one position's surface to every position on the plate.
`save_background_stack` refuses such a name, but it is worth knowing why the
convention is `_bkg.tif`.

---

## API

| | |
| --- | --- |
| `read_background_stack(path, shape=None)` | read a saved surface, optionally frame-sized → `(stack, meta)` |
| `surface_diagnostics(root, channel=None)` | every position's estimator diagnostics, from the headers alone |
| `upsample_stack(stack, shape)` | bilinear expand, frame by frame |
| `save_background_stack(corr, path)` | write one, with metadata; refuses unsafe names |
| `Surface.from_files(tif, npz)` | rebuild a corrector from what the pipeline saved |
| `Surface.from_correction(corr)` | the same from a live `PositionCorrection` |
| `Surface.correct(raw, x, y, frame)` | `(raw − B) / F`, centroids in frame pixels |
| `Surface.background_at(x, y, frame)` | just `B`, in counts |
| `correct_cell_table(cells, {ch: surface})` | add `<ch>_corrected` to a `cell_data` table |
| `track_means(cells, summary, cols)` | average per-frame columns over the mitotic window |
| `average_shape(donors)` | mean background shape, normalised to mean 1 |
| `borrow_shape(target, donors)` | donor shape, local level |
