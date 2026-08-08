"""Compare old smoothing (medfilt + closing) with the constrained Viterbi
decode on exemplary tracks from A03_s1."""
import sys
import numpy as np
import pandas as pd
import tifffile
import joblib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pathlib import Path
from scipy.signal import find_peaks

CAP = Path('/Users/ajitjoglekar/Library/CloudStorage/GoogleDrive-ajitj@umich.edu/My Drive/ImageAnalysis/cellapp-analysis')
sys.path.insert(0, str(CAP))
from dead_classifier import classify_dead
from track_decode import decode_track, track_events, STATE_NAMES

DATA = Path('/Volumes/SharedHITSX/cdb-Joglekar-Lab-GL/Soubhagyalaxmi_Jema/BUB1_OE Project/BUB1 Mutant/HeLa_PN2_3_4_5_20250709/pEN2-pEN3-pEn4-pEN5/2025-07-10/20460')
STEM = '20250709_pEN2-pEN3-pEn4-pEN5_A03_s1'
HERE = Path(__file__).parent
CACHE = HERE / 'a03s1_decoded.csv'

inf = sorted(DATA.glob(f'{STEM}_phs_*_inference'))[0]

if CACHE.exists():
    df = pd.read_csv(CACHE)
else:
    df = pd.read_excel(inf / f'{STEM}_analysis.xlsx', sheet_name='cell_data', index_col=0)
    phase = tifffile.imread(DATA / f'{STEM}_phs.tif')
    instance = tifffile.imread(inf / f'{STEM}_phs_instance.tif')
    model = joblib.load(CAP / 'models' / 'dead_classifier_hgb.joblib')['model']
    label_df = classify_dead(phase, instance, df, model)
    df['dead_flag_new'] = 0
    df['dead_proba'] = np.nan
    df.loc[label_df.index, 'dead_flag_new'] = label_df.dead_flag
    df.loc[label_df.index, 'dead_proba'] = label_df.dead_proba
    del phase, instance
    df.to_csv(CACHE, index=False)

MIN_DUR = 3  # min_mitotic_duration_in_frames

rows, decoded = [], {}
for pid, g in df.groupby('particle'):
    g = g.sort_values('frame')
    mit_obs = (g.semantic == 101).to_numpy()
    states = decode_track(mit_obs, g.dead_proba.to_numpy())
    ev = track_events(states, g.frame.to_numpy())
    _, props = find_peaks(np.append(g.semantic_smoothed.to_numpy(), np.zeros(3)),
                          width=MIN_DUR)
    decoded[pid] = (g, states, ev)
    path = ''.join('IMPD'[s] for s in states)
    rows.append(dict(pid=pid, n=len(g),
                     raw_flips=int(np.abs(np.diff(mit_obs.astype(int))).sum()),
                     old_peaks=props['widths'].size,
                     old_dur=props['widths'][0] if props['widths'].size == 1 else np.nan,
                     new_dur=int((states == 1).sum()),
                     exits=int((states == 2).any()),
                     path=''.join(c for i, c in enumerate(path) if i == 0 or c != path[i-1]),
                     label=ev['fate_label'], death_frame=ev['death_frame']))
S = pd.DataFrame(rows)

print('=== decoded fate labels (all tracks with a 101 detection) ===')
print(S.label.value_counts().to_string())
print(f"\ntracks total: {len(S)}")
print(f"old method: {(S.old_peaks > 1).sum()} tracks give >1 mitotic episode "
      f"(discarded by summarize_data); {(S.old_peaks == 0).sum()} give none")
resc = S[(S.old_peaks > 1) & (S.new_dur >= MIN_DUR)]
print(f"of those, decode recovers a single episode for {len(resc)} "
      f"({resc.label.eq('mitotic_survived').sum()} without death)")
both = S[(S.old_peaks == 1) & (S.new_dur >= MIN_DUR)]
print(f"tracks scored by both methods: {len(both)}; median duration "
      f"old {both.old_dur.median():.1f} vs decoded {both.new_dur.median():.1f} frames")

# ---- pick five distinct exemplars ----
picks, used = [], set()
def take(sel, title):
    for pid in sel.pid:
        if pid not in used:
            used.add(pid)
            picks.append((pid, title))
            return

take(S[(S.path == 'IMP') & (S.new_dur.between(4, 12)) & (S.raw_flips <= 3)
       & (S.n > 60)].sort_values('n', ascending=False),
     'interphase -> mitosis -> survives (no death): both methods agree')
take(S[(S.path == 'IMP') & (S.old_peaks > 1) & (S.new_dur >= 4)]
     .sort_values('raw_flips', ascending=False),
     'survives mitosis, but flickering semantic makes the old method see several episodes')
take(S[(S.label == 'dead_in_mitosis') & (S.n > 40)].sort_values('new_dur', ascending=False),
     'death during mitosis')
take(S[(S.label == 'dead_post_mitosis') & (S.n > 50)].sort_values('new_dur', ascending=False),
     'death after mitosis')
take(S[(S.label == 'dead_no_mitosis') & (S.n > 60)].sort_values('n', ascending=False),
     'dead-like throughout: no true mitotic episode')

state_colors = ['#eef1f4', '#7fb3e0', '#a8d5a2', '#8c8c8c']
fig, axes = plt.subplots(len(picks), 1, figsize=(12, 2.0 * len(picks)), sharex=True)
for ax, (pid, title) in zip(np.atleast_1d(axes), picks):
    g, states, ev = decoded[pid]
    f = g.frame.to_numpy()
    for t, st in zip(f, states):
        ax.axvspan(t - 0.5, t + 0.5, color=state_colors[st], lw=0)
    ax.step(f, (g.semantic == 101) * 1.0, where='mid', color='k', lw=1.1)
    ax.step(f, g.semantic_smoothed * 1.0, where='mid', color='tab:orange', lw=1.8,
            alpha=0.85)
    ax.plot(f, g.dead_proba, '.', ms=4.5, color='tab:red')
    if ev['death_frame'] is not None:
        ax.axvline(ev['death_frame'], color='crimson', ls='--', lw=1.3)
        ax.text(ev['death_frame'] + 1.5, 0.45, f"death f{ev['death_frame']}",
                fontsize=7.5, rotation=90, va='center', color='crimson')
    row = S[S.pid == pid].iloc[0]
    ax.set_ylim(-0.08, 1.15)
    ax.set_ylabel(f'particle {pid}', fontsize=9)
    ax.set_title(f'{title}   |   old: {row.old_peaks} episode(s)   '
                 f'decoded: {row.new_dur} mitotic frames, {ev["fate_label"]}',
                 fontsize=9, loc='left')
axes[-1].set_xlabel('frame')
handles = [Patch(fc=c, label=f'decoded {n}') for c, n in zip(state_colors, STATE_NAMES)]
handles += [Line2D([], [], color='k', lw=1.1, label='raw semantic == 101'),
            Line2D([], [], color='tab:orange', lw=1.8, label='old smoothed (medfilt+closing)'),
            Line2D([], [], color='tab:red', marker='.', ls='', label='classifier dead_proba')]
fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=8, frameon=False)
fig.suptitle('Track-state assignment: symmetric smoothing vs constrained Viterbi decode',
             fontsize=12)
fig.tight_layout(rect=(0, 0.07, 1, 0.985))
fig.savefig(HERE / 'track_decode_comparison.png', dpi=140)
print('\nsaved figure')
