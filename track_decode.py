'''
Constrained state decoding for particle tracks.

Biology dictates a one-way state order over a track: a cell can enter and exit
mitosis at most once, and death is absorbing. The semantic segmentation and the
dead classifier both make per-frame errors, so the per-frame evidence is decoded
into the most likely monotone state sequence

    interphase (I) -> mitotic (M) -> post-mitotic (P) -> dead (D)

by Viterbi dynamic programming, where any segment may be empty (e.g. I->M->D is
death in mitosis, I->D is a cell that was dead-like without a mitotic episode).
This replaces symmetric smoothing (median filter + closing), which cannot
enforce single-episode or absorbing-death constraints and ignores the
classifier's confidence.

Per-frame evidence:
  mitotic_obs : whether the semantic segmentation labeled the cell mitotic
  dead_proba  : classifier P(dead-like); NaN (unclassified) counts as neutral
'''
import numpy as np
import numpy.typing as npt

STATE_NAMES = ('interphase', 'mitotic', 'post_mitotic', 'dead')
# transitions allowed into each state (no backward moves)
_ALLOWED = {0: (0,), 1: (0, 1), 2: (1, 2), 3: (0, 1, 2, 3)}


def decode_track(mitotic_obs: npt.NDArray, dead_proba: npt.NDArray,
                 flip_prob: float = 0.1, dead_sem_prob: float = 0.7,
                 switch_penalty: float = 2.5, dead_weight: float = 0.3,
                 proba_clip: float = 0.05) -> npt.NDArray:
    '''
    Decode one track into the monotone state sequence I -> M -> P -> D.

    Inputs:
    mitotic_obs    : bool array; per-frame semantic-mitotic observation
    dead_proba     : float array; classifier P(dead-like), NaN = unclassified
    flip_prob      : probability of a per-frame semantic mislabel
    dead_sem_prob  : probability that a dead-like cell still carries the
                     mitotic semantic label (they usually do - that is why
                     they contaminate label 101 in the first place)
    switch_penalty : -log prior for each state transition; larger values
                     suppress short spurious episodes
    dead_weight    : tempering on the classifier evidence. The classifier is
                     overconfident and flickers between adjacent frames, so
                     its log-likelihood is down-weighted relative to the
                     semantic evidence; death must be supported by a run of
                     frames rather than by any single confident one.
    proba_clip     : dead_proba is clipped to [proba_clip, 1-proba_clip]

    Returns int array of per-frame states (0=I, 1=M, 2=P, 3=D).
    '''
    T = len(mitotic_obs)
    s = np.asarray(mitotic_obs, dtype=bool)
    q = np.asarray(dead_proba, dtype=float)
    q = np.where(np.isnan(q), 0.5, np.clip(q, proba_clip, 1.0 - proba_clip))

    # Semantic evidence: I and P expect no mitotic label, M expects one, and
    # a dead cell carries it with probability dead_sem_prob.
    sem_off = -np.log(np.where(s, flip_prob, 1.0 - flip_prob))
    sem_on = -np.log(np.where(s, 1.0 - flip_prob, flip_prob))
    sem_dead = -np.log(np.where(s, dead_sem_prob, 1.0 - dead_sem_prob))

    # Classifier evidence: only the dead state expects a high dead_proba.
    alive = -dead_weight * np.log(1.0 - q)
    dead = -dead_weight * np.log(q)

    E = np.column_stack([sem_off + alive, sem_on + alive,
                         sem_off + alive, sem_dead + dead])

    dp = np.full((T, 4), np.inf)
    ptr = np.zeros((T, 4), dtype=int)
    dp[0] = E[0]                      # any start state, no prior cost
    for t in range(1, T):
        for k in range(4):
            best_cost, best_j = np.inf, k
            for j in _ALLOWED[k]:
                c = dp[t - 1, j] + (switch_penalty if j != k else 0.0)
                if c < best_cost:
                    best_cost, best_j = c, j
            dp[t, k] = best_cost + E[t, k]
            ptr[t, k] = best_j

    states = np.zeros(T, dtype=int)
    states[-1] = int(np.argmin(dp[-1]))
    for t in range(T - 2, -1, -1):
        states[t] = ptr[t + 1, states[t + 1]]
    return states


def track_events(states: npt.NDArray, frames: npt.NDArray) -> dict:
    '''
    Extract event summary from a decoded state sequence.

    Inputs:
    states : per-frame states from decode_track
    frames : movie frame number of each track row (same length)

    Returns dict with:
    mito_start, mito_end : movie frames bounding the mitotic episode (or None)
    death_frame          : movie frame of the first dead-classified state
    death_label          : none | dead_in_mitosis | dead_post_mitosis |
                           dead_no_mitosis
    '''
    frames = np.asarray(frames)
    m = np.flatnonzero(states == 1)
    d = np.flatnonzero(states == 3)
    events = {'mito_start': int(frames[m[0]]) if m.size else None,
              'mito_end': int(frames[m[-1]]) if m.size else None,
              'death_frame': int(frames[d[0]]) if d.size else None,
              'death_label': 'none'}
    if d.size:
        if not m.size:
            events['death_label'] = 'dead_no_mitosis'
        elif np.any(states == 2):
            events['death_label'] = 'dead_post_mitosis'
        else:
            events['death_label'] = 'dead_in_mitosis'
    return events
