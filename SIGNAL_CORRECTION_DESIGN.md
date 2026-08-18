# `signal_correction.py` — design notes

Why the module is built the way it is, what was measured to decide it, and which
alternatives were tried and rejected. For how to use it, see
[SIGNAL_CORRECTION_README.md](SIGNAL_CORRECTION_README.md).

Every number here was measured on the 20250213 HT1080 pPS18 plate
(`.../20250213/HT1080 pPS18/2025-02-13/20353`): 13 sample positions, two blank
wells, 137 frames at 10 min, GFP and Texas Red, and the 431k per-cell
measurements already sitting in the `*_analysis.xlsx` files.

---

## 1. The problem

`cellaap_analysis` corrects fluorescence with two maps built from blank wells: a
FluoroBrite well gives the background map (subtracted), a DMEM well gives the
intensity map (divided by). Both are built once and applied to every position on
the plate.

The aggregate notebook found the result unconvincing: after the full, correctly
ordered correction the field came out **less** flat than the raw signal — RMS
deviation from a flat radial profile rising from 6.2% to 9.0% in GFP and 7.2% to
10.5% in Texas Red, with Texas Red no longer even a monotonic gradient. That is
the observation this module started from.

There are three distinct errors, and they are worth separating because they have
different fixes.

### 1.1 The background level is per-position

Measured cell-free background, by position, GFP:

| | counts |
| --- | --- |
| FluoroBrite blank (what gets subtracted) | 116 |
| the 13 sample positions | 118 – 128 |

The blank sits **below every sample position**. Subtracting it leaves 2–12
counts of unremoved background in GFP and up to 21 in Texas Red — against a cell
signal of only 10–16 counts. The residual is position-dependent, so it does not
cancel in a between-position comparison.

### 1.2 The level drifts, and the blank does not drift with it

Over the 137-frame movie the sample background climbs **9–15%** (first frame
118.6 → last 129.7, GFP; 146.8 → 167.9, Texas Red). The FluoroBrite well stays
flat over the same period. The background map is stored per frame, so it *looks*
like it tracks time, but it tracks the wrong well's time course.

### 1.3 The intensity map keeps the camera offset in it

This is the subtle one, and it is what makes the correction actively harmful.

`gen_intensity_correction_map` computes `smoothed_mean / smoothed_mean.max()`.
What that actually estimates is

```
(D + A·F) / (D + A·F).max()
```

which is `F` flattened towards uniform by the fraction `D / (D + A)`. How much
flattening depends on how bright the medium in that blank well was. Measured
centre-to-edge ratios of the same optical field:

| source | centre/edge | mean level |
| --- | --- | --- |
| FluoroBrite blank alone (dim) | 1.090 | 115 counts |
| DMEM blank alone (bright) | 1.190 | 534 counts |
| difference of the two (offset cancelled) | **1.220** | — |
| implied by the raw cell signal | 1.121 (see §2) | — |

Three different answers to a question with one true value, ordered exactly as
the offset-dilution argument predicts: the dimmer the well, the flatter the map
looks. The pipeline divides by the DMEM version, which is closer to right than
FluoroBrite would be but still ~3% too flat — while the background map it
subtracts first is built from the *other* well at a different level. The
combination is what produced the notebook's reversed profiles.

---

## 2. The model, and the evidence for it

```
I(x, y, t) = D + A(t)·F(x, y) + S(x, y, t)·F(x, y)
```

| symbol | | varies with |
| --- | --- | --- |
| `D` | camera offset (dark level) | nothing |
| `F` | flat field: excitation profile × collection efficiency, mean 1 | position in the field only |
| `A(t)` | medium brightness: autofluorescence, stray light | well, position, time |
| `S` | the fluorophore in the cell | what we want |

so `S = (I − B)/F` with `B(x, y, t) = D + A(t)·F(x, y)`.

Two independent checks say the model holds:

**A blank well is `D + A(t)·F` and nothing else.** Regressing every pixel of a
blank stack on the frame mean across the movie's own drift, the linear model
reproduces the measured frames to **0.12–0.14% of their mean** (0.07% for the
Texas Red DMEM well). There is no meaningful residual structure to explain.

**The flat field predicts the raw cell profile.** Taking `F` from the blank-pair
difference (centre/edge 1.220) and the implied offset (~53 counts, §4), the
model predicts a raw cell-signal falloff of `(53 + 87×1.10)/(53 + 87×0.90)` =
**1.132**. Measured over 430k cells: **1.121**. That is a 1% agreement between a
number derived entirely from blank wells and a number measured entirely from
cells.

The model is therefore not the weak point. Every failure in §1 is an estimation
failure, not a modelling one.

---

## 3. What is per-position and what is not

This is the central design decision, and getting it backwards is what the stock
pipeline does.

| quantity | scope | why |
| --- | --- | --- |
| `F` flat field | **one per plate** | a property of the optics — objective, filter, illumination path. Every position measures the same thing, so pooling is what makes it precise |
| `A(t)` medium level | **per position, per frame** | depends on well, medium depth, meniscus, evaporation, time |
| `D` camera offset | one per camera | constant, but see §4 |

So the module estimates the background from each position's own frames, and the
flat field once, from the blanks. Estimating a flat field *per position* is a
mistake — it has neither the lever arm nor the dynamic range to support it
(§7.3), and there is nothing to gain because the answer is the same everywhere.

---

## 4. The camera offset, and why nothing depends on it

`D` looks unavoidable: anything that turns a measured *background* into a flat
field has to remove it first, and the result goes as `1/(background − D)`. With
a background of ~130 counts, `D = 20` amplifies the measured falloff by
130/110 = 1.2×, while `D = 100` amplifies it by 130/30 = 4.3×.

**`D` cannot be pinned down from this data.** Three estimators were implemented
and all were removed from the shipping module; here is what each gave.

| method | GFP | Texas Red | verdict |
| --- | --- | --- | --- |
| photon-transfer (variance vs mean, x-intercept) | −35 … +18 | +1 … +24 | extrapolates below the measured range; the GFP DMEM well bleaches 70% and gives a *negative* fitted gain in some windows |
| pixelwise level regression across positions | 106 | 143 | impossible — 143 exceeds the FluoroBrite blank's own level of ~118. Biased because rising background is partly cell haze, which is flatter than `F`, compressing the slope and inflating the extrapolated intercept |
| background-vs-signal shape consistency | ~80 | ~99 | inherits every assumption about cells being equally bright everywhere |

Consistent in order of magnitude (tens of counts), useless as a number to divide
by. So the architecture was changed to not need it:

* **The background is subtracted as the measured surface `B`.** It is never
  decomposed into `D` and `A(t)`, so nothing needs `D`.
* **The flat field comes from the difference of two blanks.** Both are
  `D + A·F` with the *same* `D`, so
  `bright − dim = (A_bright − A_dim)·F` — the offset cancels algebraically
  rather than being estimated.

That is `flatfield_from_blank_pair`, and it is the reason the module works. The
supporting evidence that it is right: **the two channels agree to 1%** (1.220
GFP, 1.228 Texas Red), as an optical property must, where the single-blank maps
disagree by 9%.

Regressing each blank surface back onto that `F` implies `D` = 53 (DMEM) and 61
(FluoroBrite) counts in GFP, 76 and 79 in Texas Red — two very different wells
agreeing within a channel. That is a good consistency check, and it is the only
use `D` has here.

**If you ever acquire a real dark frame** (shutter closed, same exposure and
gain), `D` becomes known and a single blank well is enough. That is the one case
where the removed `flatfield_from_blank` / `flatfield_from_surfaces` routes
would be preferable; they are at commit `85165af`.

---

## 5. Finding cell-free regions

The brief asked to identify cell-free regions from minimum-variance regions of
the phase image. That was implemented and benchmarked against the segmentation,
and **it does not work at this confluence.**

Cell coverage on this plate is 64–78% of the field (dilated); only 2–8% of
64-pixel blocks are strictly cell-free. Scoring each detector by what fraction of
the blocks it keeps are genuinely empty:

| detector | precision among its lowest-variance 10% |
| --- | --- |
| block standard deviation (the minimum-variance idea) | 8 – 25% |
| mean Sobel gradient | 6 – 25% |
| mean \|frame − local mean\| (top-hat-like) | 8 – 26% |
| temporal standard deviation across frames | 2 – 21% |
| **chance (random blocks)** | **2 – 8%** |

A factor of 2–4 over guessing. Phase contrast is the reason: the halo is bright
*next to* cells rather than on them, and the interior of a large spread cell is
low-variance, so both error modes are systematic rather than noise. Temporal
variance is worse still, because slowly-moving cells look quiet and the halo
sweeps through nominally empty regions.

Two conclusions followed, and both shaped the module:

1. **Use the segmentation.** It is already on disk beside every position, it is
   free, and it beats every image-based detector outright. This is
   `cell_free_mask_from_labels`, and it is the default.
2. **Do not require cell-free *regions* at all.** `background_surface` takes a
   robust low-side location estimate inside each block, so it needs *some*
   cell-free pixels per block rather than a contiguous empty area, and it
   degrades smoothly as coverage rises. Blocks that fall below the floor return
   NaN and are filled from neighbours.

The phase detectors and their benchmark were removed from the shipping module
(commit `85165af`). Re-run the benchmark before assuming any of this transfers
to a sparser cell line, where the minimum-variance idea may well work fine.

---

## 6. How `estimate_position_correction` works

Five stages. The whole thing runs on 32×32 grids, which is why it is cheap.

### 6.1 Subsample frames

16–24 of 137, evenly spaced, read from a memory map. Legitimate because the
background drifts *smoothly* and monotonically — there is no frame-scale
structure to miss. This is where the runtime goes: reading frames off the share,
not arithmetic.

### 6.2 Mask cells, as far out as the frame can afford

The mask is the complement of the instance segmentation, upsampled from 1024 to
2048 and dilated.

**The dilation is the single most important parameter in the module.**
Out-of-focus light around a cell is not background; subtracting it deletes
signal. Measured on one field:

| | GFP counts |
| --- | --- |
| inside the cell mask | 139.9 |
| medium within 10 px of a cell | 132.5 |
| medium more than 30 px from any cell | 127.2 |

That 5-count halo sits against ~13 counts of actual cell signal. Sweeping the
exclusion width over the whole plate, measured on 44k cells:

| exclusion | usable blocks | recovered GFP | residual radial trend |
| --- | --- | --- | --- |
| 21 px | 75% | 12.0 counts | 33.7% |
| 61 px | 35% | 13.5 counts | 27.7% |
| 121 px | 10% | **15.9 counts** | **8.7%** |

Texas Red behaves the same way (6.6 → 9.9 counts, 40.5% → 19.5%). So the default
is 121 px — but a crowded frame cannot spare that, and a background measured
from three surviving blocks is worse than one measured closer in. The function
therefore walks a ladder (121 → 60 → 30 → 15 → no mask) **per frame**, stopping
at the first width where `min_usable_blocks` survive. On this plate positions
land anywhere from `[121, 121]` to `[60, 60]`, the latter being the densest
field.

The mask uses `ndi.maximum_filter(size=d)`, the separable form of dilating by a
square: bit-identical to `binary_dilation` with a square footprint and ~100×
faster (0.03 s against 4.9 s at 121 px on a 2048 frame). Without that, a
per-frame ladder would not be affordable.

### 6.3 Reduce each frame to a 32×32 grid

Per-block **median of the unmasked pixels**, 64×64 blocks, NaN for blocks with
fewer than `min_pixels` usable. Median rather than mean so that anything the
segmentation missed — debris, a fragment — is ignored rather than averaged in.

Without a segmentation the reduction switches to `clipped_mean`: iterative
sigma-clipping of the bright tail only, with the width taken from a MAD measured
on the **low half** of the block so cells cannot inflate the threshold they are
about to be compared against. A plain low quantile was tried and rejected — the
10th percentile sits 12–25 counts below the true cell-free median, which is
10–17% of the background and fatal at this contrast.

This is a 4096× data reduction, and it is what makes everything downstream
instant.

### 6.4 Turn the grid into a surface

Default `'grid'`: nearest-neighbour fill of NaN blocks, then a Gaussian over the
32×32 grid. See §7.1 for the constrained alternative and why it lost.

### 6.5 Factor into a shape and a per-frame level

```
B(x, y, t) = background_offset[t] + background_level[t] · background_shape(x, y)
```

The shape is the mean surface normalised to mean 1. The level is a **least-
squares projection** of each frame's surface onto that shape, not a frame mean —
so filled-in blocks contribute according to how well they match the shape rather
than dragging the average. The 16 sampled levels are then linearly interpolated
across all 137 frames, which is safe precisely because the drift is smooth.

The `offset` term exists so the `'flatfield'` model can be stored exactly; under
`'grid'` it stays zero.

Result: a 32×32 shape, a 32×32 flat field, and two numbers per frame — ~50 kB
per position per channel, against a full-size map stack.

---

## 7. Alternatives tried and rejected

### 7.1 A constrained background: fit `offset + level·F` per frame

Two free parameters instead of ~1000 grid points. This *should* win when the
cells are faint — fewer parameters, less noise, no smoothing bias, no
nearest-neighbour fill at the corners where the flat field is furthest from 1,
and a fitted camera offset for free.

It lost, clearly:

| model | GFP residual radial trend | Texas Red |
| --- | --- | --- |
| grid + smooth | **8.7%** | **19.5%** |
| `offset + level·F` | 24.7% | 29.8% |

The reason is that the background is genuinely **not** proportional to `F`:
out-of-focus haze follows local cell density, which has its own spatial
structure, and a flexible grid can follow that where two parameters cannot.

So it is kept as a **fallback, not an improvement**: `background_model='auto'`
(the default) uses the grid and drops to the fit only when a position is too
crowded to measure a grid at all.

**Where the crossover sits was measured.** Fitting each model to a random
subset of the measured blocks and predicting the held-out ones:

| blocks kept | grid + smooth | `offset + level·F` |
| --- | --- | --- |
| 2% | 6.74 | **6.01** |
| 5% | **5.73** | 5.78 |
| 12% | **4.88** | 5.64 |
| 40% | **4.06** | 5.57 |

(counts RMS, median over 13 positions × 14 frames × 2 channels). The grid keeps
improving as blocks accumulate; the fit plateaus near 5.6 counts, which is what
two free parameters buys. Hence `grid_min_blocks = 0.05`.

**Measured on real positions.** Six positions spanning the confluence range,
both channels, 41k cell measurements:

* *The fallback is dormant on this plate.* All twelve default runs chose the
  grid, at 9.1–30.8% usable blocks — every one above the trip point. The
  dilation ladder gets there first: its own floor (`min_usable_blocks`) is also
  0.05, so the ladder opens the exclusion until 5% of blocks survive, and only
  a field too crowded for its widest rung can drop the model below the
  threshold. The switch is a genuine last resort, not a routine path.
* *Forced on at normal density the fit is worse*, reproducing the
  cross-validation on the endpoint that matters rather than on held-out blocks:

  | pooled, 41k cells | grid | flatfield |
  | --- | --- | --- |
  | GFP radial RMS | **11.4%** | 22.7% |
  | GFP spatial swing | **5.3 counts** | 7.9 counts |
  | Texas Red radial RMS | **81.5%** | 120.8% |
  | Texas Red spatial swing | **18.2 counts** | 22.0 counts |

* *The two models agree on the population but not on the cell.* Median
  difference −0.07 counts (GFP) and −0.52 (Texas Red) — no systematic offset,
  which is what makes a mixed plate safe in aggregate — but the IQR is ±2–3
  counts and 48–60% of cells move by more than 2. Against a 14-count (GFP) or
  9-count (Texas Red) signal that is real per-cell scatter.
* *The switch and the refusal both fire when starved.* Disabling the ladder at
  a 161 px exclusion drove 10 of 12 runs below 5% blocks and all ten switched,
  logging the reason. At 201 px six runs hit the terminal case and raised the
  explicit "no sampled frame has the 10 measurable blocks" error rather than
  inventing a surface.

One result points at a possible improvement rather than a problem. Under the
starved configuration a few positions came out markedly *better* — Texas Red
C08_s5 went from 9.3 to 26.4 counts of signal with negatives falling from 29%
to 0.2%, and A08_s1's residual radial trend fell from 368% to 76% — because a
161 px exclusion keeps far more halo out of the background than the 60–121 px
the ladder settles on. Others got worse, as ten fitted blocks should. It
suggests the ladder may be trading halo cleanliness for block count, and that
*wide exclusion + two-parameter fit* could beat *narrow exclusion + grid* on
crowded positions. That is a hint, not a result: it would need the same
held-out comparison run per rung before any default changed.

### 7.2 Fitting the flat field to the cells

If cells are on average equally bright everywhere in the field, the trend in
background-subtracted cell signal across the field *is* the illumination. This
needs no blank well, no cell-free region, and no camera offset — it is the most
attractive route on paper, and the right answer for a plate with no blanks.

It fails here for a mundane reason: **there is not enough contrast**. Cells sit
~10% above their local background (140 counts inside the mask against 127 in the
gaps), leaving ~12 counts of signal. The 8×8 bin medians then scatter over
0.30–1.75 and the "flat field" is noise; on Texas Red it produced negative bins,
which would flip the sign of the corrected signal. Cross-validation across
positions confirmed it does not generalise (centre/edge 1.065 fitted on one half
of the positions, 1.144 on the other — for a quantity that should be identical).

Removed from the shipping module at commit `85165af`. Worth restoring for a
plate with a bright reporter and no blank wells; it should carry an SNR guard
that refuses rather than returning noise.

### 7.3 Recovering the flat field from a position's own drift

Each position's background level moves 8–21% through the movie, and
`I = D + A(t)·F` means a pixelwise regression on the frame level recovers `F` as
the slope. On a blank well this works beautifully (model residual 0.1%, clean
slopes spanning 0.72–1.12).

On a sample position it produces garbage: fitted slopes from −0.6 to +3.2, where
the true range is roughly 0.7–1.2. The lever arm is too short relative to the
noise on per-block estimates, and cells growing and dividing change the block
statistics in step with the drift, confounding the regression. This is the
concrete evidence behind §3's claim that the flat field must not be estimated
per position.

### 7.4 Temporal low-quantile background

A low quantile through time at each pixel, after normalising each frame's level.
Needs no cell-free region in any single frame — only that each pixel is empty at
*some* point in the movie. Attractive for very high confluence.

Not used because the segmentation-masked block estimate is strictly better when
a segmentation exists, which it always does here, and because cells that do not
move leave their footprint in the result. Removed at commit `85165af`; the right
fallback if this is ever applied to data without segmentations.

---

## 8. What it buys, measured

Against the 431k per-cell measurements, 13 positions, both channels. The metric
is residual spatial bias **in counts** — how far a cell's reported signal moves
across the field purely because of where it sat.

| stage | median | radial swing | RMS |
| --- | --- | --- | --- |
| GFP raw | 139.6 | 15.4 counts | 8.2 |
| GFP blank-well maps | 24.9 | 4.2 counts | 2.5 |
| GFP this module | 15.9 | **4.0 counts** | **1.4** |
| Texas Red raw | 164.7 | 20.5 counts | 11.7 |
| Texas Red blank-well maps | 14.7 | 5.7 counts | 1.6 |
| Texas Red this module | 9.9 | **2.5 counts** | 1.9 |

Percentages are *not* comparable between stages whose medians differ tenfold —
the corrected stages are ~10× smaller in magnitude, so an equal absolute error is
a 10× larger percentage. This is why the aggregate notebook's RMS-percent metric
made the corrections look worse than they were.

On the spatial axis the module clearly beats the stock maps on Texas Red
(5.7 → 2.5 counts) and roughly ties on GFP (4.2 → 4.0 swing, though RMS improves
2.5 → 1.4). **The larger effect is the level**, per §1.1–1.2: the stock
correction leaves 2–21 counts of position-dependent background in place.

### 8.1 Effect on the biology

Applied to the 626 tracks the aggregate notebook analyses, per-track signal
averaged over each track's mitotic window:

| channel | correction | Spearman ρ vs mitotic duration | within-position ρ | median signal |
| --- | --- | --- | --- | --- |
| GFP | blank-well maps | +0.024 (p = 0.54) | +0.008 | 72.8 |
| GFP | this module | +0.006 (p = 0.89) | +0.039 | 55.8 |
| Texas Red | blank-well maps | +0.433 (p = 5e−30) | +0.447 | 56.6 |
| Texas Red | this module | **+0.440** (p = 5e−31) | **+0.463** | 42.0 |

The conclusions are unchanged: Texas Red predicts mitotic duration, GFP does
not. The correction slightly strengthens Texas Red and leaves GFP
non-significant. Spearman is rank-based, so modest per-position offsets barely
move it — but **median signal falls 23% (GFP) and 26% (TR)**, and anything using
absolute signal (the Hill fit's EC50, cross-plate comparisons) shifts by that
much.

The reconstruction of per-track values from the per-frame table was validated
against the summary's own raw column first: r = 0.995 / 0.996, median difference
−0.5% / −0.4%.

---

## 9. Known limits

* **~10% contrast is the binding constraint.** Cells sit 13 counts above their
  local background out of 140. A 1-count error in the background is a ~10% error
  in the cell, and residual systematics of ±1.5 counts remain. About 5% of Texas
  Red tracks now come out slightly negative — real cells at their local
  background, previously hidden behind the stock correction's under-subtraction.
  No estimator gets under this floor; a brighter reporter or longer exposure
  would.
* **Haze is handled bluntly.** Excluding a wide annulus around every cell is a
  crude way to deal with out-of-focus light. Modelling it — a density-dependent
  term convolved with the out-of-focus PSF — would recover the blocks currently
  discarded, which matters most in exactly the crowded fields where the ladder
  has to back off.
* **The dilation ladder is a heuristic.** It optimises usable-block count, not
  bias. A frame that backs off to 60 px has a measurably different background
  estimator from one that stayed at 121, and that difference is not currently
  propagated as an uncertainty.
* **No uncertainty is reported.** Every corrected value is a point estimate. The
  per-block scatter and the number of surviving blocks would support an error
  bar, and at this contrast an error bar would be genuinely useful.
* **BaSiC has not been compared against.** It estimates darkfield jointly — the
  one quantity this module has to design around — but its central assumption
  (sparse foreground, low-rank background) is violated at 64–78% confluence.
  Worth settling empirically; if its darkfield lands near the 53–79 counts
  implied here, that is strong mutual corroboration.

---

## 10. Where the removed code went

The module originally carried 21 public functions covering the investigation as
well as the result. It was trimmed to the 12 on the correction path. Everything
removed is at commit `85165af` and its findings are recorded above:

| removed | covered in |
| --- | --- |
| `cell_free_mask_from_phase`, `benchmark_cell_free_masks` | §5 |
| `flatfield_from_signal` | §7.2 |
| `flatfield_from_blank`, `flatfield_from_surfaces` | §4 (need a real dark frame) |
| `darkfield_from_ptc`, `darkfield_from_shape_consistency`, `darkfield_from_flatfield`, `level_regression` | §4, §7.3 |
| `temporal_background` | §7.4 |
