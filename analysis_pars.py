import os
import pandas as pd
import numpy as np
from skimage.morphology import disk
import trackpy

class analysis_pars:

    def __init__(self, cell_type = "hela"):
        self.cell_types = ["hela", "u2os", "rpe1", "ht1080"]

        try:
            if cell_type.lower() in self.cell_types:
                self.current_cell_type = cell_type
        except:
            raise ValueError(f"Choose one of {self.cell_types}")
        
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

        # trackpy parameters
        self.max_pixel_movement = 20
        self.tracking_memory    = 1
        self.min_track_length   = 10 # min track length

        self.track_mode = "vanilla"

        if cell_type.lower() == "ht1080":
            self.max_pixel_movement = 30
            self.max_cell_size = 9000
            self.track_mode = "predictive"

        if cell_type.lower() == "u2os":
            self.max_pixel_movement = 30
            self.max_cell_size = 9500
            self.track_mode = "predictive"

        if cell_type.lower() == "rpe1":
            self.max_pixel_movement = 30
            self.max_cell_size = 9500
            self.track_mode = "predictive"