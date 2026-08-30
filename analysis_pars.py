from pathlib import Path
import numpy as np
from skimage.morphology import disk
import joblib

class analysis_pars:

    def __init__(self, cell_type = "hela", frame_interval = None):
        self.cell_types = ["hela", "u2os", "rpe1", "ht1080"]

        try:
            if cell_type.lower() in self.cell_types:
                self.current_cell_type = cell_type
        except:
            raise ValueError(f"Choose one of {self.cell_types}")
        
        # Mask value for mitotic cells. Different cellaap versions write
        # different ones; mitotic_semantic_values lists every value that means
        # mitotic, and mitotic_mask_value is the one this dataset actually
        # uses (resolved against the segmentation in cellaap_analysis).
        self.mitotic_semantic_values = (100, 101)
        self.mitotic_mask_value = 101

        # Parameters for pre-processing inferred cell segmentations
        self.erode_footprint = disk(3)
        self.max_cell_size = 4000
        self.min_cell_size = 250

        # Shortest run of mitotic frames that counts as a mitosis. This is
        # biology, so it is held in MINUTES and converted to frames against the
        # acquisition's own interval by set_frame_interval. There is no default
        # interval: it is a property of the experiment, and the 10 that used to
        # sit here was silently wrong for datasets acquired at 4 min, which put
        # this threshold at 12 minutes rather than 30.
        self.min_mitotic_duration = 30 # minutes
        self.frame_interval = None                   # min, from the metadata
        self.min_mitotic_duration_in_frames = None   # set by set_frame_interval
        if frame_interval is not None:
            self.set_frame_interval(frame_interval)

        # Smoothing of the semantic label trace, in FRAMES rather than minutes,
        # and deliberately not converted by set_frame_interval. This closes
        # single-frame flicker in the segmentation's output - a classifier
        # error is one frame wide whether frames are 4 or 10 minutes apart - so
        # unlike min_mitotic_duration it does not scale with the interval.
        # Must be odd so it can be centered symmetrically on each pixel; otherwise the operation will translate the peak
        self.semantic_gap_closing = 3 # number of frames
        self.semantic_footprint = np.ones(self.semantic_gap_closing)

        # How far past each mitotic episode the dead/mitotic classifier is run.
        # A cell that dies on exiting mitosis loses the mitotic semantic label
        # while it is still rounded, so the frames that carry the death are
        # just past the peak. Scoring the whole rest of the track instead does
        # not work - by then the cell is flat, which this model reads as "dead"
        # whether or not it is (see the dead_classifier docstring).
        self.post_peak_frames = 10

        # Spurious daughter tracks. The segmentation labels anaphase mitotic,
        # so trackpy opens a fresh particle on a dividing cell: the track
        # begins late in the movie and is mitotic almost at once. Rejecting
        # tracks that are mitotic on their very first frame catches most, but
        # not the ones whose first frame or two are not yet labeled. A track is
        # also rejected when it starts after late_track_start_fraction of the
        # movie AND reaches mitosis within early_mitosis_frames of its start.
        # Both conditions are needed: a late start alone is ordinary, and so is
        # an early mitosis in a track that has been followed from the outset.
        self.late_track_start_fraction = 1 / 3
        self.early_mitosis_frames = 3

        # Spurious-track filter (summarize_data). Two independent causes put a
        # short "mitosis" in the summary, and they need different tests.
        #
        # A track that is only watched briefly after mitotic entry cannot show
        # a long mitosis whatever the cell does - at an observation window of
        # 50 frames, 96% of tracks report a short one regardless of how well
        # they were tracked. So the window after entry must be long enough to
        # contain a typical mitosis for this experiment. "Typical" is measured
        # from the data rather than set here: a fixed threshold cannot serve
        # both an arrest (median mitosis ~4 h) and an unperturbed control
        # (~40 min), and neither can one scaled only by the frame interval.
        # The reference is the median mitotic duration of tracks followed for
        # reference_track_fraction of the movie, and the window must be at
        # least min_window_factor of that.
        #
        # The other cause is trackpy opening a fresh particle late in the
        # movie on a fragment in a crowded area, which no window test catches:
        # holding the window fixed, tracks starting in the last two thirds of
        # the movie report a short mitosis 24-63% of the time against 2-9% for
        # the rest. That step is what max_track_start_fraction cuts, as a
        # fraction of the movie so it does not depend on its length.
        #
        # Note this threshold is the same as late_track_start_fraction but is
        # applied on its own rather than in conjunction with an early mitosis,
        # so it subsumes that test while the filter is enabled. The conjunction
        # is kept for when it is not.
        # summarize_data runs one position at a time, so the reference duration
        # is measured per position and the resulting threshold varies about
        # +/-20% across a plate - noise, since positions in a well share their
        # biology. Set reference_duration_frames to pin one value plate-wide
        # (the median over the whole plate's reference tracks) when that
        # variation matters; left None, each position calibrates itself.
        self.exclude_short_window_tracks = True
        self.reference_track_fraction = 0.75
        self.reference_duration_frames = None
        self.min_window_factor = 0.65
        self.max_track_start_fraction = 1 / 3

        # Death call. A cell is dead from the first frame of the first run of
        # death_run_frames consecutive frames with P(dead) above the
        # threshold. High threshold plus a required run keeps isolated
        # confident frames - classifier flicker - from ending a mitosis.
        self.death_proba_threshold = 0.8
        self.death_run_frames = 5

        # Border exclusion. The classifier needs a full CROP_SIZE (96 px, full
        # resolution) box around the centroid, so a cell within half of that of
        # the edge can never be scored and its mitotic frames carry no evidence
        # either way. Margin is in analysis-scale pixels, i.e. 96 / 2 / 2.
        self.border_margin = 24
        self.exclude_border_tracks = True

        # trackpy parameters
        self.max_pixel_movement = 20
        self.tracking_memory    = 1

        # Shortest track trackpy keeps (filter_stubs, in track_centroids).
        # This is a fragment filter, not a data filter. Since summarize_data
        # gained the episode-length, entry-observed and observation-window
        # rules it removes nothing from the summary: dropping it to 1 on two
        # positions of the 20260826 plate left the summary identical at 102
        # and 314 rows. What it does remove is the fragments every stage
        # between tracking and the summary would otherwise carry - the tracked
        # table goes 1043 -> 3200 and 2939 -> 7008 without it, which is +18-28%
        # mitotic tracks for measure_signal to scan the table for, and the same
        # again through the dead classifier, both workbooks and the particle
        # browser. It is also the unconditional floor for when
        # _filter_spurious_tracks switches itself off - too few reference
        # tracks to calibrate against, or exclude_short_window_tracks cleared.
        #
        # In FRAMES, and deliberately not scaled by the frame interval, like
        # semantic_gap_closing: a tracking fragment is a fragment whether
        # frames are 4 or 10 minutes apart. It sits below what the window test
        # already demands on the datasets checked - the shortest summarized
        # track is 11 frames at 4 min/frame and 15 at 10 - but it becomes the
        # binding constraint wherever min_window falls under it, which needs a
        # short mitosis AND a long interval (an unperturbed control at 10
        # min/frame). Compare it against the shortest track_length in the
        # summary before trusting it on such a dataset.
        self.min_track_length   = 10 # frames

        self.track_mode = "vanilla"

        _cell_type_overrides = {
            "ht1080": {"max_pixel_movement": 30, "max_cell_size": 9000, "track_mode": "predictive"},
            "u2os":   {"max_pixel_movement": 30, "max_cell_size": 9500, "track_mode": "predictive"},
            "rpe1":   {"max_pixel_movement": 30, "max_cell_size": 9500, "track_mode": "predictive"},
        }
        for attr, val in _cell_type_overrides.get(cell_type.lower(), {}).items():
            setattr(self, attr, val)

        #### Mitotic vs dead discrimination model
        # ResNet18-512 -> PCA(32) -> logistic regression, trained on 826 hand
        # labels pooled from the BUB1 (345) and CycB-oe 20576 (481) datasets
        # over 15 microscope positions. Class 1 = mitotic, 0 = dead.
        # Leave-one-position-out balanced accuracy 0.878, AUC 0.945, ECE 0.018.
        #
        # Keep the whole bundle: it declares which feature blocks the estimator
        # consumes (this model takes the 512-d embedding alone, the previous
        # one took embedding + 16 handcrafted), and classify_dead builds the
        # design matrix from that declaration.
        self.dead_classifier_bundle = joblib.load(
            Path(__file__).parent / "models" / "dead_classifier_pooled.joblib")
        self.dead_classifier = self.dead_classifier_bundle["model"]

    def set_frame_interval(self, minutes):
        '''Record the acquisition's frame interval and convert the durations.

        Kept separate from __init__ because the interval is discovered from the
        data - `analysis.files` reads it out of the acquisition metadata once
        it knows which folder the position lives in - while the parameters
        object is built before that. Calling it again with a different value
        re-derives cleanly, so a caller may override what the metadata claims.

        min_mitotic_duration_in_frames keeps the floor of the old `//`: the
        threshold is "at least this many frames", and rounding 30 min at 4
        min/frame up to 8 rather than down to 7 would tighten it beyond what
        changing the interval is meant to do. At 4 min/frame it becomes 7
        frames where the hard-coded 10 gave 3.
        '''
        minutes = float(minutes)
        if not np.isfinite(minutes) or minutes <= 0:
            raise ValueError(f"frame_interval must be a positive number of "
                             f"minutes, got {minutes}")
        self.frame_interval = minutes
        self.min_mitotic_duration_in_frames = max(
            1, int(self.min_mitotic_duration // minutes))
        return self

    def require_frame_interval(self):
        '''Fail with an actionable message rather than on a None comparison.

        Everything that reads min_mitotic_duration_in_frames goes through here
        first, so an unset interval is caught at the top of the stage that
        needs it instead of surfacing as a TypeError deep inside _episodes.
        '''
        if self.frame_interval is None or self.min_mitotic_duration_in_frames is None:
            raise ValueError(
                "frame_interval is not set, so the minimum mitotic duration "
                "cannot be converted to frames. Build the analysis through "
                "analysis.files()/from_analysis_file (which read it from the "
                "acquisition metadata), or call "
                "defaults.set_frame_interval(minutes) yourself.")
        return self.min_mitotic_duration_in_frames
