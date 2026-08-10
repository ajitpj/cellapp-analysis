import os, re, tifffile, warnings
from pathlib import Path, PurePosixPath, PureWindowsPath
from skimage.io import imread # type: ignore
from skimage.morphology import erosion, closing
from skimage.measure import regionprops_table
import numpy as np
import pandas as pd
import trackpy as tp
import scipy.ndimage as ndi
from scipy.signal import medfilt
from analysis_pars import analysis_pars
from cellaap_utils import *
from dead_classifier import classify_dead, rows_to_classify

class analysis:
    
    def __init__(self, root_folder: Path, plotting_only: False):
        '''
        Object initializes with default parameter values and definitions of 
        root and inference folders. It also reads in tif files with either 
        intensity or background in their names as the corresponding, channel-
        specific correction maps.
        '''
        # self.__dict__.update((key, False) for key in self.suffixes)
        # self.__dict__.update((key, value) for key, value in file_dict.items() if key in self.suffixes)

        ## Default parameters
        self.paths = {}  # Dictionary stores stack paths
        self.stacks = {} # Dictionary stores image stacks
        self.quality = {} # Dictionary to store quality checks
        ##
        try:
            root_folder.exists()
        except:
            raise ValueError(f"{root_folder} not a valid path")
        
        self.root_folder = root_folder
        self.inference_folders = [] # list of folders with _inference in their names
        for directory in os.scandir(root_folder):
            if "_inference" in directory.name:
                self.inference_folders.append(directory)

        ## Check for and set path names for correction maps
        # These maps are loaded with the root directory because they apply to all 
        # segmentations in the root directory.
        self.intensity_map_present  = False
        self.background_map_present = False
        
        if not plotting_only:
            self._load_maps()
        else:
            print(f"Opening {root_folder} in plotting only mode.")

    def _load_maps(self, ):
        '''
        Function to detect and load background and intensity correction maps from the root folder
        '''
        paths = [os.path.join(dirpath,f) for (dirpath, _, filenames) in os.walk(self.root_folder) for f in filenames]
        maps_types_chnls = [
            (
                name, re.search(r"background|intensity", name).group(), re.search(r"GFP|Texas Red|Cy5|phs", name).group()
                ) for name in paths if re.search(r"background|intensity", name)
            ]
        for name, type, channel_name in maps_types_chnls:
            # if channel_name  == "Texas Red":
            #     channel_name = "Texas_Red"
            match type:
                case "intensity":
                    self.paths[channel_name + "_intensity_map"]  = Path(name)
                    self.stacks[channel_name + "_intensity_map"] = tifffile.imread(Path(name))
                    self.intensity_map_present = True
                    print(f"{name} used as the {channel_name} intensity map")

                case "background":
                    self.paths[channel_name + "_background_map"]  = Path(name)
                    self.stacks[channel_name + "_background_map"] = imread(Path(name))
                    self.background_map_present = True
                    print(f"{name} used as the {channel_name} background map")

    def files(self, cellaap_dir: Path, cell_type: str):
        '''
        Inputs:
        cellaap_dir: directory containing cellapp inference; must contain "instance" and "semantic" tif files
        cell_type: specify the cell type so appropriate default pars are set
        '''
        # Process the path objects to retrieve the parent directories and suffixes
        try:
            self.cellaap_dir = cellaap_dir
            self.paths["instance"] = Path([name for name in cellaap_dir.glob('*.tif') if "instance" in name.name][0])
            self.paths["semantic"] = Path([name for name in cellaap_dir.glob('*.tif') if "semantic" in name.name][0])
            # Keep the name stub to infer other file names
            self.name_stub = re.search(r"[A-H]([1-9]|[0][1-9]|[1][0-2])_s(\d{2}|\d{1})", str(self.paths["semantic"].name)).group()
            self.expt_name = self.paths["instance"].name.split(f"{self.name_stub}")[0]
            self.defaults = analysis_pars(cell_type=cell_type)
        except:
            raise ValueError("Instance and/or semantic segmentations not found!")

        # Set path names for existing channel files
        self.data_dir = self.cellaap_dir.parent
        image_paths = list(self.data_dir.glob('*.tif')) + list(self.data_dir.glob('*.tiff'))
        image_paths = map(str, image_paths)
        valid_paths = [
            name for name in image_paths if re.search(r"(.tif|.tiff)", name) and self.name_stub+"_" in name
        ]
        for name in valid_paths:
            channel = re.search(r"GFP|Texas Red|Cy5|phs", name)
            if channel:
                print(f"{str(name)} can be used for analysis or display")
                match channel.group():
                    case "phs":
                        self.paths["phase"] = Path(name)
                    case "GFP":
                        self.paths["GFP"] = Path(name)
                    case "Texas Red":
                        self.paths["Texas Red"] = Path(name)
                    case 'Cy5':
                        self.paths["Cy5"] = Path(name)
                
        # Only read instance and semantic stacks
        self.stacks["phase"]    = imread(self.paths["phase"]) #For mitotic/dead discrimination
        self.stacks["semantic"] = imread(self.paths["semantic"])
        self.stacks["instance"] = imread(self.paths["instance"])
        
        # Apply erosion before zooming to avoid duplicate computations
        instance_shape = self.stacks["instance"].shape
        # Saving the number of planes for track filtering (summarize_data)
        self.max_timepoints = instance_shape[0]

        # int16 to match the erosion below, which already casts to int16. The
        # default float64 costs 8 bytes per pixel of a 4x-upsampled stack -
        # 15 GB for a 450-frame 1024x1024 segmentation, enough to swap the
        # machine - and buys nothing, since these are integer instance labels
        # only ever used in the equality test in measure_signal.
        instance_zoomed = np.zeros((instance_shape[0], instance_shape[1]*2, instance_shape[2]*2),
                                   dtype=np.int16)
        print(f"Computing zoomed and eroded instance mask...")
        for i in np.arange(instance_shape[0]):
            pre_zoom = erosion(self.stacks["instance"][i,:,:].astype(np.int16),
                               self.defaults.erode_footprint)
            instance_zoomed[i, :, :] = ndi.zoom(pre_zoom, 2, order=0)
            
        self.stacks["instance_zoomed"] = instance_zoomed
        print(f"Finished computing zoomed and eroded instance mask!")
        return
    

    def create_correction_maps(self, type_file_dict: dict):
        '''
        Function uses the assigned background and intensity mapping stacks
        to create the corresponding background and intensity correction maps.
        the position_name_channel dictionary should include key-value pairs as:
        key: intensity or background , value: Path_to_file
        Currently, has to be used for each channel individually.
        '''

        for key, value in type_file_dict.items():
            if type(value) is PurePosixPath or PureWindowsPath or Path:
                if "intensity" in key:
                    intensity_map_name = value.parent / Path(value.stem + "_intensity_map.tif")
                    channel_map = value.stem.split('_')[-1] + "_intensity_map"
                    self.stacks[channel_map] = gen_intensity_correction_map(imread(str(value)))
                    tifffile.imsave(intensity_map_name, self.stacks[channel_map].astype(np.float16))
                    self.paths[channel_map] = intensity_map_name
                    print(f"Intensity map saved in the data dir. as {intensity_map_name}")

                elif "background" in key:
                    background_map_name = value.parent / Path(value.stem + "_background_map.tif")
                    channel_map = value.stem.split('_')[-1] + "_background_map"
                    self.stacks[channel_map] = gen_background_correction_map(imread(str(value)))
                    tifffile.imsave(background_map_name, self.stacks[channel_map].astype(np.int16))
                    self.paths[channel_map] = background_map_name
                    print(f"Background map saved in the data dir. as {background_map_name}")
            
            else:
                print(f"Dictionary values must be Path objects; try again.")
        self._load_maps()
        return


    def track_centroids(self, mode: str, save_flag: False, memory = None, max_pixel_movement = None) -> pd.DataFrame:
        '''
        Function uses the instance segmentation file generated by cell_aap and trackpy
        to track cell centroids.
        Inputs:
        saveflag  : Save the pandas dataframe as a csv file.
        Optionally, track memory and min. track length can be changed.
        Outputs:
        csv_name  : name of the csv file saved
        dataframe : dataframe with tracking information
        '''
        
        frames = self.stacks["instance"].shape[0] # type: ignore

        property_list  = ['centroid','area','eccentricity','bbox', 'label']

        df_list = []

        for i in np.arange(frames):
            props = regionprops_table(self.stacks["instance"][i,:,:], properties=property_list)
            props_table = pd.DataFrame(props)
            props_table["centroid-0"] = props_table["centroid-0"].apply(lambda x: int(x))
            props_table["centroid-1"] = props_table["centroid-1"].apply(lambda x: int(x))
            props_table["frame"] = i #trackpy needs this to track.
            df_list.append(props_table)

        track_table = pd.concat(df_list)

        # trackpy requirements for the dataframe it will use to link tracks
        track_table.rename(columns={"centroid-0":"x", "centroid-1":"y"}, inplace=True)
        # remove large and small cells
        track_table = track_table[(track_table.area < self.defaults.max_cell_size) & 
                                  (track_table.area > self.defaults.min_cell_size)]

        # Check if default values have been changed
        if memory is None:
            memory = self.defaults.tracking_memory
        if max_pixel_movement is None:
            max_pixel_movement = self.defaults.max_pixel_movement

        match mode:

            case "predictive":
                # tp.linking.Linker.MAX_SUB_NET_SIZE = 40
                track_pred = tp.predict.NearestVelocityPredict(span=3)
                self.tracked = track_pred.link_df(track_table,
                                                  max_pixel_movement,
                                                  adaptive_stop=10, 
                                                  adaptive_step=0.9, 
                                                  memory=memory)

            case "adaptive":
                self.tracked = tp.link(track_table, 
                                       max_pixel_movement, 
                                       adaptive_stop=10, 
                                       adaptive_step=0.9,
                                       memory=memory)
            
            case "vanilla":
                self.tracked = tp.link(track_table, 
                                       max_pixel_movement, 
                                       memory=memory)
        
        
        self.tracked = tp.filter_stubs(self.tracked, self.defaults.min_track_length) 
        
        # tracked annoying use the frame number as the row index.
        # So drop the frame column and reset the index to recreate it.
        self.tracked.drop(columns=["frame"], inplace=True)
        self.tracked.reset_index(inplace=True)
        # For displaying tracks in naapri
        self.tracked.sort_values(by=["particle","frame"],inplace=True)

        # Obtain the semantic label
        # This drops the old index and uses serial numbers
        self.tracked.reset_index(inplace=True)

        semantic_label = []
        for i in np.arange(len(self.tracked)):
            semantic_label.append(self.stacks["semantic"][self.tracked.loc[i,"frame"],
                                                      self.tracked.loc[i,"x"],  
                                                      self.tracked.loc[i,"y"]])

        self.tracked["semantic"] = semantic_label

        # Different cellaap versions write different mitotic mask values (101 in
        # this pipeline, 100 in others). If the configured value is absent the
        # analysis silently finds no mitosis anywhere, so resolve it against
        # what the segmentation actually contains, and fail loudly when that is
        # ambiguous rather than returning an empty result.
        observed = set(np.unique(semantic_label).tolist())
        if self.defaults.mitotic_mask_value not in observed:
            present = [v for v in self.defaults.mitotic_semantic_values
                       if v in observed]
            if len(present) == 1:
                print(f"mitotic_mask_value={self.defaults.mitotic_mask_value} is "
                      f"absent from this segmentation; using {present[0]} instead "
                      f"(observed values {sorted(observed)}).")
                self.defaults.mitotic_mask_value = present[0]
            else:
                raise ValueError(
                    f"mitotic_mask_value={self.defaults.mitotic_mask_value} does not "
                    f"occur in the semantic segmentation, which contains {sorted(observed)}. "
                    f"Set analysis_pars.mitotic_mask_value to the mitotic value for "
                    f"this dataset, otherwise no mitotic events will be detected.")

        # remove 0's and 2's, and fill gaps in the semantic vector.
        self.tracked.loc[self.tracked.semantic != self.defaults.mitotic_mask_value, "semantic"] = 1

        self.tracked.loc[:, "semantic"] = medfilt(self.tracked.semantic,
                                                  self.defaults.semantic_gap_closing)

        # Turn into Boolean. Compare against the mask value directly; the former
        # (semantic - 1)//99 only happened to work for values of 100 or 101.
        semantic_smoothed = (self.tracked.semantic
                             == self.defaults.mitotic_mask_value).astype(int)
        semantic_smoothed = closing(semantic_smoothed,
                                    self.defaults.semantic_footprint)
        self.tracked["semantic_smoothed"] = semantic_smoothed

        
        ###########################################################################
        # Mitotic/dead discrimination
        # Detections whose semantic label is mitotic, plus post_peak_frames
        # after each episode, are classified from the instance-masked phase
        # crop. The tail catches a cell dying as it leaves mitosis, while it is
        # still rounded. Every other row keeps NaN probabilities: the model was
        # never asked about it, and a 0 there would read as a confident "not
        # dead". Do not widen the tail much - far from mitosis the cell is flat
        # and this model reads flat as dead (see the dead_classifier docstring).
        # remembered for the border test in summarize_data, which may run later
        # from a saved analysis file with no stack in memory
        self.frame_shape = tuple(self.stacks["instance"].shape[-2:])

        n_to_classify = len(rows_to_classify(
            self.tracked, self.defaults.mitotic_semantic_values,
            self.defaults.post_peak_frames))
        print(f"Classifying {n_to_classify} detections (mitotic + "
              f"{self.defaults.post_peak_frames} frames after each episode)...")
        label_df = classify_dead(self.stacks["phase"],
                                 self.stacks["instance"],
                                 self.tracked,
                                 self.defaults.dead_classifier_bundle,
                                 mitotic_values=self.defaults.mitotic_semantic_values,
                                 post_peak_frames=self.defaults.post_peak_frames)
        self.tracked["mitotic_proba"] = np.nan
        self.tracked["dead_proba"] = np.nan
        self.tracked["dead_flag"] = 0
        self.tracked.loc[label_df.index, "mitotic_proba"] = label_df.mitotic_proba
        self.tracked.loc[label_df.index, "dead_proba"] = label_df.dead_proba
        self.tracked.loc[label_df.index, "dead_flag"] = label_df.dead_flag

        # The classifier is the only consumer of the full-resolution phase
        # stack, and at 2048x2048x450 it is ~3.8 GB that would otherwise sit
        # alongside the channel stack measure_signal loads next. Drop it here;
        # measure_signal and _display_tracks both re-read from self.paths.
        self.stacks.pop("phase", None)
        ###########################################################################

        
        # classify the cells as dividing or non-dividing
        # observed division = 1; no division = 0
        for id in list(set(self.tracked.particle)):
            index  = self.tracked[self.tracked.particle==id].index

            if np.isin(self.defaults.mitotic_mask_value, self.tracked[self.tracked.particle==id].semantic):
                self.tracked.loc[index, "mitotic"] = 1
            else:
                self.tracked.loc[index, "mitotic"] = 0
        
        if save_flag:
            self.tracked.to_excel(self.cellaap_dir / Path(self.expt_name+self.name_stub+"_tracks.xlsx"))

        return self.tracked
    
    def measure_signal(self, channel: str, save_flag: False, id = -1,):
        '''
        Measures the average cell signal over the eroded cell masks. Also
        calculates the position-dependent correction factors for background fluorescence and 
        excitation intensity variation. 
        Inputs:
        channel   - Channel Name; must be one of: (phase, GFP, Texas_Red, Cy5)
        save_flag - Whether to export dataframe as xlsx
        id        - ID of the cell (assigned by trackpy); -1 will analyze all cells
        '''
        # try:
        #     if channel in []:
        #         pass
        # except:
        #     raise ValueError(f"")
        
        # Read image stack
        channel_stack = imread(self.paths[channel])

        ##
        # if particle id is set to -1 - measure all particles. 
        # otherwise, just the specified particle.
        if id > -1:
            if type(id) == int:
                try:
                    if id in set(self.tracked.particle):
                        id_list = []
                        id_list.append(id)
                except:
                    raise ValueError(f"id not in the tracking list")
            
            elif type(id) == list:
                id_list = id
        else:
            id_list = sorted(list(set(self.tracked[self.tracked.mitotic==1].particle)))
            print(f"{len(id_list)} tracks to process")
        
        # Default values for all entries
        self.tracked[channel] = np.nan
        self.tracked[channel+"_int_corr"] = 1.
        self.tracked[channel+"_bkg_corr"] = 0. 
        
        
        for id in id_list:
            
            # Measurement decision
            measure_cell = True
            print(f"Processing cell #{id}...")
            
            if measure_cell:
                frames = self.tracked[self.tracked.particle==id].frame.tolist()
                labels = self.tracked[self.tracked.particle==id].label.tolist()
                index  = self.tracked[self.tracked.particle==id].index
                
                signal = np.zeros(len(frames))
                background_correction = np.zeros_like(signal)
                intensity_correction = np.ones_like(signal)
                counter = 0
                for f,l in zip(frames, labels):
                    # print("Processing frame # {f}...")
                    mask = self.stacks["instance_zoomed"][f,:,:]==l

                    signal[counter] = mean_signal_from_mask(channel_stack[f,:,:], mask)
                    
                    if self.background_map_present:
                        map_name = channel + '_background_map'
                        background_correction[counter] = mean_signal_from_mask(self.stacks[map_name][f,:,:].astype(float), mask)
                    
                    if self.intensity_map_present:
                        map_name = channel + '_intensity_map'
                        intensity_correction[counter] = mean_signal_from_mask(self.stacks[map_name], mask)

                    counter = counter + 1
                
                self.tracked.loc[index, channel] = signal
                self.tracked.loc[index, channel+"_bkg_corr"] = background_correction
                self.tracked.loc[index, channel+"_int_corr"] = intensity_correction


        if save_flag:
            with pd.ExcelWriter(self.cellaap_dir / Path(self.expt_name+self.name_stub+'_analysis.xlsx')) as writer:  
                self.tracked.to_excel(writer, sheet_name='cell_data')
                # There are no scalars - so turn into list; transform
                pd.DataFrame([self.paths]).T.to_excel(writer,   sheet_name='file_data')
                pd.DataFrame([self.defaults.__dict__]).T.to_excel(writer,sheet_name='parameters')
            # self.tracked.to_excel()

        return self.tracked
    
    
    def _analysis_frame_shape(self):
        '''(rows, cols) of the analysis-scale frame, for the border test.

        Taken from the instance stack when the object was built by the normal
        pipeline. When summarize_data is driven from a saved analysis file
        there is no stack in memory, so fall back to a caller-set frame_shape,
        and finally to the largest bounding box in the table - across hundreds
        of frames some cell always touches the edge, so that maximum is the
        frame size to within a pixel or two.
        '''
        shape = getattr(self, 'frame_shape', None)
        if shape is None:
            stacks = getattr(self, 'stacks', None)
            if stacks is not None and stacks.get('instance') is not None:
                shape = tuple(stacks['instance'].shape[-2:])
        if shape is None:
            cols = self.tracked.columns
            if {'bbox-2', 'bbox-3'}.issubset(cols):
                shape = (int(self.tracked['bbox-2'].max()),
                         int(self.tracked['bbox-3'].max()))
                print(f'frame shape not available; inferred {shape} from bounding '
                      f'boxes for the border test')
            else:
                raise ValueError('cannot determine frame shape for the border '
                                 'test; set self.frame_shape = (rows, cols)')
        return shape

    @staticmethod
    def _runs(flags):
        '''(start, stop) half-open row ranges of each True run in a 1-D mask.'''
        flags = np.asarray(flags, dtype=bool)
        padded = np.concatenate(([False], flags, [False]))
        edges = np.flatnonzero(padded[1:] != padded[:-1])
        return list(zip(edges[::2], edges[1::2]))

    def _episodes(self, sem):
        '''Mitotic episodes of a track, as half-open row ranges.

        A run of mitotic-labeled frames shorter than min_mitotic_duration_in_frames
        is segmentation noise, not a mitosis. Runs are read directly off the
        trace rather than through find_peaks, which reports peak *bases* - the
        frame before the episode - and needed an off-by-one correction.
        '''
        return [(s, e) for s, e in self._runs(sem)
                if e - s >= self.defaults.min_mitotic_duration_in_frames]

    def _death_row(self, dead_proba):
        '''Row of the first frame at which the cell is called dead, or None.

        The rule: P(dead) above death_proba_threshold for death_run_frames
        consecutive frames, and the death is dated to the first frame of that
        run. Unscored frames (NaN) are not evidence of death and break a run -
        the model was never asked about them.

        Death is irreversible, so a run the cell visibly recovers from was a
        transient burst, not a death, and is skipped: the classifier can read a
        cell as dead for a few frames as it rounds up, then correct itself.
        "Recovers" means death_run_frames consecutive scored frames back under
        1 - death_proba_threshold after the run. Without this, particle 87 of
        E10_s7 died on the first frame of its mitosis at P(dead) > 0.85, was
        back at 0.01 five frames later, ran a second mitosis 300 frames on at
        P(dead) = 0.00, and still reported mitosis = 0.
        '''
        p = np.asarray(dead_proba, dtype=float)
        over = np.where(np.isnan(p), False, p > self.defaults.death_proba_threshold)
        alive = np.where(np.isnan(p), False,
                         p < 1.0 - self.defaults.death_proba_threshold)
        need = max(int(self.defaults.death_run_frames), 1)
        for s, e in self._runs(over):
            if e - s < need:
                continue
            if any(b - a >= need for a, b in self._runs(alive[e:])):
                continue                      # recovered: not a death
            return int(s)
        return None

    def summarize_data(self, save_flag: True):
        '''
        Summarizes data stored in the tracked dataframe; operates on all measured channels.

        Two rules do the work, both per track:

        * the mitotic episodes are the runs of mitotic semantic label at least
          min_mitotic_duration_in_frames long. The FIRST episode is the mitosis
          that gets reported; n_peaks says how many there were.
        * the cell is dead from the first frame of the first run of
          death_run_frames consecutive frames with P(dead) above
          death_proba_threshold, discarding runs the cell visibly recovers
          from - death is irreversible (_death_row). Where that frame falls
          relative to the first episode gives the fate:

          no qualifying run        -> mitotic_survived
          before the cell has been -> interphase death. The cell never had a
          mitotic for                 mitosis - it rounded up because it was
          min_mitotic_duration_       dying - so there is no mitotic duration
          in_frames                   and no mitotic signal to report, and the
                                      track is left OUT of the summary.
          inside the first episode -> dead_in_mitosis. Reports the frame of
                                      death, the time in mitosis before it
                                      (time_to_death), and signal averaged over
                                      those pre-death mitotic frames.
          after the first episode  -> dead_post_mitosis. The mitosis completed,
                                      so the full mitotic duration and its
                                      signal are reported, plus the frame of
                                      death.

        Death is only visible where the classifier ran: the mitotic frames plus
        defaults.post_peak_frames AFTER each episode. Nothing before mitotic
        entry is ever scored - this model was built to separate mitotic from
        dead cells, both rounded, and has no meaning on an interphase cell.
        The tail is what catches a cell dying on mitotic exit, and is
        deliberately short for the same reason.

        A cell is summarized only if it has at least one episode and (default)
        none of its mitotic-labeled detections sit within defaults.border_margin
        of the frame edge, where the classifier cannot see the whole cell.

        Columns: mito_start and death_frame are absolute movie frames, and
        directly comparable. mitosis is a count of frames - the length of the
        first episode as observed, never truncated, so a cell that rounds up
        and dies still reports how long it was seen rounded. time_to_death is
        the frames from mitotic entry to the death call (NaN if the cell never
        dies). Fluorescence is averaged over the episode up to the death call
        only, so it is NaN for a cell dead from its first mitotic frame.

        Note what never reaches this function: tracks shorter than
        defaults.min_track_length are removed by trackpy in track_centroids,
        so a cell that rounds up, dies and loses its track inside that window
        is gone before any of this runs.

        Inputs -
        save_flag : whether to export the data as an xlsx file

        Outputs -
        None
        '''
        # Select only those tracks where mitosis was observed
        idlist    = list(set(self.tracked[self.tracked.mitotic==1].particle))

        # Border test: any mitotic-labeled detection closer than the margin
        # means the classifier could not score the cell there.
        margin = self.defaults.border_margin
        rows, cols = self._analysis_frame_shape()
        near_border_ids = set()
        if margin:
            mit_rows = self.tracked[self.tracked.semantic_smoothed == 1]
            near = ((mit_rows.x < margin) | (mit_rows.x > rows - margin) |
                    (mit_rows.y < margin) | (mit_rows.y > cols - margin))
            near_border_ids = set(mit_rows.loc[near, 'particle'].unique())
        if self.defaults.exclude_border_tracks and near_border_ids:
            dropped = len(near_border_ids & set(idlist))
            idlist = [i for i in idlist if i not in near_border_ids]
            print(f'excluded {dropped} tracks with mitotic detections within '
                  f'{margin} px of the frame edge ({rows}x{cols})')
        
        # A list to store the number of peaks
        # Multiple peaks will reveal either tracking errors or segmentation issues
        peaks_per_track = np.zeros(len(idlist)) 
        # Fluctuations in the mask size will indicate segmentation quality
        cell_area_std  = np.zeros_like(peaks_per_track)

        mitosis          = []
        time_to_death    = [] # frames from mitotic entry to the death call
        mito_start       = []
        cell_area        = []
        particle         = []
        track_length     = []
        channels         = []
        max_displacement = []
        dead_cell_score  = [] # keep track of "dead" flags
        fate_label       = [] # mitotic_survived / dead_in_mitosis / dead_post_mitosis
        death_frame      = [] # movie frame of death (NaN if the cell never dies)
        n_peaks          = [] # mitotic episodes in the track
        n_sem_mitotic    = [] # frames the segmentation called mitotic
        n_scored         = [] # frames the classifier actually scored
        n_interphase_death = 0 # died without ever having a mitosis; excluded

        # Check which channels have been measured. If none, return only "mitotic duration"
        # Need to find a better way to code this.

        if "GFP" in self.tracked.columns:
            channels.append('GFP')
        if "Texas Red" in self.tracked.columns:
            channels.append("Texas Red")
        if "Cy5" in self.tracked.columns:
            channels.append("Cy5")
        signal_storage = {}
        for channel in channels:
            signal_storage[f'{channel}'] = []
            signal_storage[f'{channel}_std'] = []
            signal_storage[f'{channel}_bkg_corr'] = []
            signal_storage[f'{channel}_bkg_corr_std'] = []
            signal_storage[f'{channel}_int_corr'] = []
            signal_storage[f'{channel}_int_corr_std'] = []

        for index, id in enumerate(idlist):

            track_rows = self.tracked[self.tracked.particle==id]
            frames    = track_rows.frame.to_numpy()
            sem_raw   = (track_rows.semantic_smoothed == 1).to_numpy()
            dead_flag = track_rows.dead_flag.to_numpy()
            proba     = track_rows.dead_proba.to_numpy(dtype=float)

            episodes = self._episodes(sem_raw)
            peaks_per_track[index] = len(episodes)
            cell_area_std[index]   = track_rows.area.std()
            if not episodes:
                continue

            # The first episode is the mitosis being reported; any later one is
            # a re-rounding, most often the cell dying after it divided.
            first_start, first_stop = episodes[0]
            death = self._death_row(proba)

            if death is None:
                fate = 'mitotic_survived'
            elif death - first_start < self.defaults.min_mitotic_duration_in_frames:
                # The cell was called dead before it had been mitotic for the
                # minimum duration, so it never had a mitosis - it rounded up
                # because it was dying. That is an interphase death, not a
                # mitotic event: no mitotic duration and no mitotic signal to
                # report, so the track is left out of the summary rather than
                # carried as a row of NaNs. Same threshold that decides what
                # counts as an episode in the first place.
                n_interphase_death += 1
                continue
            elif death < first_stop:
                fate = 'dead_in_mitosis'
            else:
                fate = 'dead_post_mitosis'

            # Fluorescence window: the first episode, ending at death if the
            # cell died during it, so no signal is measured from a cell already
            # called dead. A cell dead from its first mitotic frame leaves this
            # empty and its channel means are NaN - there was no live mitotic
            # frame to measure.
            window = np.zeros(len(frames), dtype=int)
            window[first_start:first_stop] = 1
            if death is not None:
                window[death:] = 0

            # Duration is the observed episode, NOT truncated at death: a cell
            # that rounds up and dies was still seen rounded for that long, and
            # zeroing it made those rows unusable. time_to_death carries the
            # truncated quantity - frames from mitotic entry to the death call,
            # which for a post-mitotic death runs past the end of the episode.
            mitosis.append(int(first_stop - first_start))
            time_to_death.append(int(death - first_start) if death is not None
                                 else np.nan)
            mito_start.append(int(frames[first_start]))
            death_frame.append(int(frames[death]) if death is not None else np.nan)
            fate_label.append(fate)
            n_peaks.append(len(episodes))
            dead_cell_score.append(np.sum(sem_raw*dead_flag))
            n_sem_mitotic.append(int(sem_raw.sum()))
            n_scored.append(int(np.isfinite(proba).sum()))
            cell_area.append(track_rows.area.mean())
            particle.append(id)
            track_length.append(len(frames))
            coords = self.tracked.loc[self.tracked.particle==id, ['x', 'y']]
            disp_vector = calculate_displacement(coords)
            max_displacement.append(np.max(disp_vector))

            # A cell that dies on the first frame of its mitosis has an empty
            # averaging window, so the channel means are legitimately NaN;
            # numpy's empty-slice warnings would otherwise flood the log.
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                for channel in channels:
                    signal, bkg_corr, int_corr, area, signal_std, bkg_std, int_std, area_std = calculate_signal(
                                                                window,
                                                                self.tracked[self.tracked.particle==id][f'{channel}'].to_numpy(),
                                                                self.tracked[self.tracked.particle==id][f'{channel}_bkg_corr'].to_numpy(),
                                                                self.tracked[self.tracked.particle==id][f'{channel}_int_corr'].to_numpy(),
                                                                self.tracked[self.tracked.particle==id].area.to_numpy(),
                                                                self.defaults.semantic_footprint
                                                                )
                    signal_storage[f'{channel}'].append(signal)
                    signal_storage[f'{channel}_std'].append(signal_std)
                    signal_storage[f'{channel}_bkg_corr'].append(bkg_corr)
                    signal_storage[f'{channel}_bkg_corr_std'].append(bkg_std)
                    signal_storage[f'{channel}_int_corr'].append(int_corr)
                    signal_storage[f'{channel}_int_corr_std'].append(int_std)
                    # signal_storage[f'{channel}_area_mean'].append(area)
                    # signal_storage[f'{channel}_area_std'].append(area_std)


        n_dead = sum(1 for f in fate_label if f != 'mitotic_survived')
        print(f'summarized {len(particle)} tracks; {n_dead} died during or '
              f'after mitosis (P(dead) > {self.defaults.death_proba_threshold} '
              f'for {self.defaults.death_run_frames} consecutive frames)')
        if n_interphase_death:
            print(f'excluded {n_interphase_death} tracks dead within '
                  f'{self.defaults.min_mitotic_duration_in_frames} frames of '
                  f'mitotic entry (interphase death - no mitosis to measure)')

        # Quality metrics
        n_obs, bins = np.histogram(peaks_per_track, np.arange(0,15))
        peaks_per_track_df = pd.DataFrame({"n_peaks"     : bins[:-1],
                                           "track number" : n_obs})
        
        self.quality["peaks_per_track"] = peaks_per_track_df
        self.quality["cell_area_std"]   = pd.DataFrame(cell_area_std, columns = ["Area std."])

        # Construct summary DF
        other_storage = {
                        "particle"         : particle,
                        "track_length"     : track_length,
                        "max_displacement" : max_displacement,
                        "mito_start"       : mito_start,
                        "cell_area"        : cell_area,
                        "mitosis"          : mitosis,
                        "time_to_death"    : time_to_death,
                        "dead_cell_score"  : dead_cell_score,
                        "fate_label"       : fate_label,
                        "death_frame"      : death_frame,
                        "n_peaks"          : n_peaks,
                        "n_sem_mitotic"    : n_sem_mitotic,
                        "n_scored"         : n_scored
                        }
        
        summary_storage = other_storage | signal_storage
        self.summaryDF = pd.DataFrame(summary_storage)


        if save_flag:
            with pd.ExcelWriter(self.cellaap_dir / Path(self.expt_name+self.name_stub+"_summary.xlsx")) as writer: 
                self.summaryDF.to_excel(writer,sheet_name = "Summary", index=False)
                pd.DataFrame([self.paths]).T.to_excel(writer, sheet_name='file_data')
                pd.DataFrame([self.defaults.__dict__]).T.to_excel(writer,sheet_name='parameters')
                self.quality["cell_area_std"].join(self.quality["peaks_per_track"]).to_excel(writer, 
                                                                           sheet_name="quality", 
                                                                           index=False)

        return self.summaryDF
    