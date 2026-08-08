import os
from pathlib import Path
import pandas as pd
import numpy as np
from skimage.morphology import disk
import trackpy
import joblib

class analysis_pars:

    def __init__(self, cell_type = "hela"):
        self.cell_types = ["hela", "u2os", "rpe1", "ht1080"]

        try:
            if cell_type.lower() in self.cell_types:
                self.current_cell_type = cell_type
        except:
            raise ValueError(f"Choose one of {self.cell_types}")
        
        # Mask value for mitotic cells
        self.mitotic_mask_value = 101

        # Parameters for pre-processing inferred cell segmentations
        self.erode_footprint = disk(3)
        self.max_cell_size = 4000
        self.min_cell_size = 250

        # Median filter size for smoothing semantic label trace
        self.min_mitotic_duration = 30 # minutes
        self.frame_interval = 10 #  time step in min
        self.min_mitotic_duration_in_frames = self.min_mitotic_duration // self.frame_interval

        # Must be odd so it can be centered symmetrically on each pixel; otherwise the operation will translate the peak
        self.semantic_gap_closing = 3 # number of frames
        self.semantic_footprint = np.ones(self.semantic_gap_closing)

        # Constrained track decoding (interphase -> mitotic -> post-mitotic -> dead)
        self.decode_flip_prob      = 0.1 # per-frame semantic mislabel probability
        self.decode_dead_sem_prob  = 0.7 # dead cell still carries the mitotic label
        self.decode_switch_penalty = 2.5 # -log prior per state transition; larger
                                         # values suppress short spurious episodes
        self.decode_dead_weight    = 0.3 # tempering on classifier evidence; death
                                         # must be supported by a run of frames

        # trackpy parameters
        self.max_pixel_movement = 20
        self.tracking_memory    = 1
        self.min_track_length   = 10 # min track length

        self.track_mode = "vanilla"

        _cell_type_overrides = {
            "ht1080": {"max_pixel_movement": 30, "max_cell_size": 9000, "track_mode": "predictive"},
            "u2os":   {"max_pixel_movement": 30, "max_cell_size": 9500, "track_mode": "predictive"},
            "rpe1":   {"max_pixel_movement": 30, "max_cell_size": 9500, "track_mode": "predictive"},
        }
        for attr, val in _cell_type_overrides.get(cell_type.lower(), {}).items():
            setattr(self, attr, val)

        #### Dead cell discrimination model
        # HistGradientBoosting on handcrafted features; 1 = live mitotic, 0 = dead-like
        self.dead_classifier = joblib.load(
            Path(__file__).parent / "models" / "dead_classifier_hgb.joblib")["model"]