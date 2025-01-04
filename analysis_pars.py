import os
import pandas as pd
from skimage.morphology import disk

class analysis_pars:

    def __init__(self, cell_type = "hela"):
        self.cell_types = ["hela", "u2os", "rpe1", "ht1080"]

        try:
            if cell_type.lower() in self.cell_types:
                self.current_cell_type = cell_type
        except:
            raise ValueError(f"Choose one of {self.cell_types}")
        
        self.erode_footprint = disk(7)
        self.max_cell_size = 4000
        self.min_cell_size = 500
        self.min_mitotic_duration = 2 # Number of frames
        self.delta_t = 10 # default time step in min

        # 5 - max. centroid movement between frames, memory - if the segmentation flickers
        self.max_pixel_movement = 20
        self.tracking_memory    = 2
        self.min_track_length   = 10 # min track length

        if cell_type.lower() == "ht1080":
            self.max_pixel_movement = 60