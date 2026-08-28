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

# Fraction of a cell's own mask that must carry the mitotic value for the
# detection to be called mitotic. The measurement is close to binary - a cell's
# mask is typically >99% one class or the other - so anything from 0.1 to 0.7
# gives the same answer; 0.5 is the honest middle.
MITOTIC_MASK_FRACTION = 0.5


def _resolve_semantic_classes(occurs, defaults):
    """Which values in this segmentation are classes, and which one is mitotic.

    cellaap writes one value per class - 0 background, 1 or 99 interphase, 100
    or 101 mitotic, depending on version - but where two masks overlap it
    writes their SUM: 1 + 100 = 101, 100 + 100 = 200, 100 + 101 = 201. A sum is
    not a class. The pixel belongs to two cells at once and says nothing about
    either, so everything outside the class set is treated as background: the
    sums, and any stray value a deeper pile-up produces (1 + 1 + 100 = 102).

    There are only ever three class values, so they are identified rather than
    filtered: background, the mitotic value, and the commonest value left,
    which is interphase - a field is mostly interphase, and a contact rim is
    thinner than the cells it lies between.

    The mitotic value comes from analysis_pars.mitotic_semantic_values, which
    lists what the different cellaap versions write. When more than one of
    those occurs, a sum is told from a class by arithmetic and by rarity: a sum
    is the total of two other present values that are each more common than it
    is. That separates the two datasets that look alike from the outside - 101
    alongside a much larger 1 and 100 is a contact rim, while 101 alongside 100
    with no 1 present cannot be built from anything else and is a class in its
    own right.

    Takes the value histogram of the whole stack. Returns (mitotic_value, keep)
    where keep is a Boolean mask over the value axis, True for the three class
    values and False for everything else.
    """
    present = np.flatnonzero(occurs)
    present = present[present > 0]

    candidates = [v for v in defaults.mitotic_semantic_values
                  if v < len(occurs) and occurs[v] > 0]
    if not candidates:
        raise ValueError(
            f"none of {defaults.mitotic_semantic_values} occurs in this "
            f"semantic segmentation, whose values are "
            f"{sorted(int(v) for v in present)}. Set "
            f"analysis_pars.mitotic_semantic_values for this dataset, "
            f"otherwise no mitotic events will be detected.")

    def is_sum(v):
        parts = present[present < v]
        parts = parts[occurs[parts] > occurs[v]]
        return bool(np.isin(v - parts, parts).any())

    real = [v for v in candidates if not is_sum(v)]
    if not real:
        raise ValueError(
            f"every candidate mitotic value {candidates} looks like the sum of "
            f"two overlapping masks rather than a class of its own (pixel "
            f"counts {[int(occurs[v]) for v in candidates]}). Set "
            f"analysis_pars.mitotic_semantic_values for this dataset.")
    if len(real) == 1:
        mitotic_value = real[0]
    else:
        # Both candidates are classes in their own right, so this dataset uses
        # both. The mitotic one is the minority - a field is mostly interphase
        # - unless analysis_pars names it outright.
        if defaults.mitotic_mask_value in real:
            mitotic_value = defaults.mitotic_mask_value
        else:
            mitotic_value = min(real, key=lambda v: occurs[v])
        print(f"semantic values {real} are all listed as mitotic (pixel counts "
              f"{[int(occurs[v]) for v in real]}); reading {mitotic_value} as "
              f"the mitotic class. Set analysis_pars.mitotic_mask_value if "
              f"that is wrong.")

    keep = np.zeros(len(occurs), dtype=bool)
    keep[0] = True                                  # background
    keep[mitotic_value] = True
    # Interphase is the commonest value left that is not itself a sum. The
    # guard matters where mitotic cells crowd together: the region two of them
    # share can outweigh what is left of either, and without it that sum would
    # be promoted to a class and both cells scored on it.
    others = [v for v in present if v != mitotic_value and not is_sum(v)]
    if others:
        keep[max(others, key=lambda v: occurs[v])] = True
    return int(mitotic_value), keep


def mask_semantic_values(instance, semantic, frames, labels, defaults,
                         threshold: float = MITOTIC_MASK_FRACTION):
    """Each detection's semantic label, read over its whole instance mask.

    The obvious rule - read the semantic frame at the centroid - is wrong for
    two kinds of cell. A concave or crescent-shaped mask does not contain its
    own centroid, so the lookup lands on background or on the neighbour; and
    where two cells touch, the segmentation writes their SUM (199 = 99 + 100
    for an interphase cell overlapping a mitotic one), which no lookup can
    decode. Both show up as isolated wrong frames in a track.

    So each detection is scored over its own mask instead: the fraction of the
    cell's pixels carrying the mitotic value. On real data that fraction is
    close to binary - mean 0.998 for the detections a centroid calls mitotic
    and 0.000 for the rest - so the threshold does no real work.

    Overlap sums count as background (see _resolve_semantic_classes), which
    keeps them out of the denominator as well as the numerator. Left in the
    denominator they would dilute an overlapped cell below the threshold: two
    touching mitotic cells write 100 + 100 = 200 across the whole region they
    share, and both would read as interphase. Each cell is judged on the
    pixels that belong to it alone; a detection with no such pixels left is
    the only one dropped, and it is counted in the diagnostics.

    Returns (values, mask_px, diagnostics). `values` goes into
    `tracked.semantic` and is then smoothed by `_label_semantic` exactly as a
    centroid-derived column was; `mask_px` is each detection's mask size.
    """
    frames = np.asarray(frames, dtype=np.int64)
    labels = np.asarray(labels, dtype=np.int64)

    # Resolve which values are classes and which one means mitotic against the
    # segmentation itself, the same way _label_semantic does against the table
    # - different cellaap versions write 100 or 101, and guessing wrong finds
    # no mitosis at all - while telling a class apart from the sum two
    # overlapping masks write.
    top = int(semantic.max())
    if top > 4095:
        raise ValueError(f"semantic segmentation holds values up to {top}; this "
                         f"is not a class map and cannot be scored over masks")
    width = top + 1
    occurs = np.zeros(width, dtype=np.int64)
    for f in range(len(semantic)):
        occurs += np.bincount(semantic[f].reshape(-1), minlength=width)
    mitotic_value, keep = _resolve_semantic_classes(occurs, defaults)

    # Per frame, the histogram of semantic values under each instance label, in
    # one bincount over label*width + value. Doing it per detection instead
    # would compare a full 1024x1024 frame per row - hundreds of thousands of
    # passes over the stack for a single position.
    counts = np.zeros((len(frames), width), dtype=np.int64)
    for f in np.unique(frames):
        lab = instance[f].reshape(-1).astype(np.int64)
        val = semantic[f].reshape(-1).astype(np.int64)
        n = int(lab.max()) + 1
        hist = np.bincount(lab * width + val,
                           minlength=n * width).reshape(n, width)
        rows = np.flatnonzero(frames == f)
        inside = labels[rows] < n
        counts[rows[inside]] = hist[labels[rows[inside]]]

    # Overlap sums are background as far as the score is concerned. Zeroing
    # their columns here is the same as zeroing those pixels in the stack, but
    # without copying it - a value the mask never counts is a value the mask
    # never saw.
    before = counts.sum(axis=1)
    counts[:, ~keep] = 0
    dropped_row = before - counts.sum(axis=1)
    dropped_px = int(dropped_row.sum())
    rows_with_dropped = int((dropped_row > 0).sum())

    mask_px = counts.sum(axis=1) - counts[:, 0]      # value 0 is background
    mitotic_px = counts[:, mitotic_value]
    fraction = mitotic_px / np.maximum(mask_px, 1)
    is_mitotic = (mask_px > 0) & (fraction > threshold)

    # Non-mitotic rows keep the most common real value under the mask rather
    # than a sentinel, so tracked.semantic still reads as segmentation output.
    # _label_semantic collapses all of them to 1 anyway.
    other = counts.copy()
    other[:, 0] = 0
    other[:, mitotic_value] = 0
    values = np.where(other.any(axis=1), other.argmax(axis=1), 1)
    values[is_mitotic] = mitotic_value

    diagnostics = {
        "mitotic_value": int(mitotic_value),
        "mask_fraction_threshold": float(threshold),
        "rows": int(len(frames)),
        "rows_mitotic": int(is_mitotic.sum()),
        "rows_empty_mask": int((mask_px == 0).sum()),
        "overlap_values": {int(v): int(occurs[v])
                           for v in np.flatnonzero(~keep) if occurs[v]},
        "overlap_px_in_masks": dropped_px,
        "rows_with_overlap_px": rows_with_dropped,
    }
    return values, mask_px, diagnostics


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
                    self.stacks[channel_name + "_background_map"] = tifffile.imread(Path(name))
                    self.background_map_present = True
                    print(f"{name} used as the {channel_name} background map")

    def files(self, cellaap_dir: Path, cell_type: str, frame_interval = None):
        '''
        Inputs:
        cellaap_dir: directory containing cellapp inference; must contain "instance" and "semantic" tif files
        cell_type: specify the cell type so appropriate default pars are set
        frame_interval: minutes between frames. Read from the acquisition
                        metadata in the data directory when omitted; pass it
                        explicitly when that metadata is wrong or absent.
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

        # The frame interval comes from the acquisition, not from the analysis
        # defaults, so it is resolved here - once data_dir is known - rather
        # than carried as a constant in analysis_pars. It decides what counts
        # as a mitotic episode at all, so an explicit argument wins over the
        # metadata, and unreadable metadata is fatal rather than silently
        # replaced by a guess.
        if frame_interval is None:
            frame_interval = read_frame_interval(self.data_dir)
            print(f"frame interval {frame_interval:g} min, from the "
                  f"acquisition metadata in {self.data_dir.name}")
        self.defaults.set_frame_interval(frame_interval)
        print(f"minimum mitotic duration {self.defaults.min_mitotic_duration} "
              f"min = {self.defaults.min_mitotic_duration_in_frames} frames")

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
                    self.stacks[channel_map] = gen_intensity_correction_map(tifffile.imread(str(value)))
                    tifffile.imwrite(intensity_map_name, self.stacks[channel_map].astype(np.float16))
                    self.paths[channel_map] = intensity_map_name
                    print(f"Intensity map saved in the data dir. as {intensity_map_name}")

                elif "background" in key:
                    background_map_name = value.parent / Path(value.stem + "_background_map.tif")
                    channel_map = value.stem.split('_')[-1] + "_background_map"
                    self.stacks[channel_map] = gen_background_correction_map(tifffile.imread(str(value)))
                    tifffile.imwrite(background_map_name, self.stacks[channel_map].astype(np.int16))
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

        # Scored over each detection's whole instance mask rather than at its
        # centroid: a concave mask need not contain its own centroid, and
        # touching cells are written as the SUM of their classes (199 = 99 +
        # 100), which no single-pixel lookup can decode. This costs about a
        # second more per position than the lookup it replaces - its cost is
        # per pixel where the lookup's was per detection - which is nothing
        # against the rest of the stage. It is for robustness, not speed.
        values, _, sem_diag = mask_semantic_values(
            self.stacks["instance"], self.stacks["semantic"],
            self.tracked["frame"].to_numpy(), self.tracked["label"].to_numpy(),
            self.defaults)
        self.tracked["semantic"] = values
        self.quality["semantic_from_masks"] = sem_diag
        # The value resolved against the segmentation is the one the rest of
        # the analysis must use, so _label_semantic has nothing left to guess.
        self.defaults.mitotic_mask_value = sem_diag["mitotic_value"]
        if sem_diag["overlap_values"]:
            print(f'semantic values {sorted(sem_diag["overlap_values"])} are '
                  f'sums written where masks overlap; treated as background '
                  f'({sem_diag["overlap_px_in_masks"]} px inside masks, '
                  f'affecting {sem_diag["rows_with_overlap_px"]} detection(s))')
        if sem_diag["rows_empty_mask"]:
            print(f'{sem_diag["rows_empty_mask"]} detection(s) had no pixels '
                  f'under their instance label; semantic left non-mitotic')

        self._label_semantic()


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

        if save_flag:
            self.tracked.to_excel(self.cellaap_dir / Path(self.expt_name+self.name_stub+"_tracks.xlsx"))

        return self.tracked
    
    def _label_semantic(self):
        '''Turn the raw semantic label into the columns the summary needs.

        Adds semantic_smoothed (gaps closed, 1 = mitotic) and mitotic (1 for a
        track that was ever mitotic), and collapses semantic to the mitotic
        mask value against 1. Shared by track_centroids and
        from_analysis_file, since an analysis file written before these
        columns existed needs exactly the same treatment.
        '''
        # Different cellaap versions write different mitotic mask values (101 in
        # this pipeline, 100 in others). If the configured value is absent the
        # analysis silently finds no mitosis anywhere, so resolve it against
        # what the segmentation actually contains, and fail loudly when that is
        # ambiguous rather than returning an empty result.
        observed = set(np.unique(self.tracked.semantic).tolist())
        if self.defaults.mitotic_mask_value not in observed:
            present = [v for v in self.defaults.mitotic_semantic_values
                       if v in observed]
            if not present:
                raise ValueError(
                    f"mitotic_mask_value={self.defaults.mitotic_mask_value} does not "
                    f"occur in the semantic segmentation, which contains {sorted(observed)}. "
                    f"Set analysis_pars.mitotic_mask_value to the mitotic value for "
                    f"this dataset, otherwise no mitotic events will be detected.")
            # On the tracking path mask_semantic_values has already resolved
            # this and left only class values in the column. An analysis file
            # written before it did can still hold both - a real class and the
            # sum of two overlapping masks - so take the rarer, which is the
            # mitotic class either way: a field is mostly interphase, and a
            # contact rim is thinner still.
            chosen = min(present,
                         key=lambda v: int((self.tracked.semantic == v).sum()))
            print(f"mitotic_mask_value={self.defaults.mitotic_mask_value} is "
                  f"absent from this segmentation; using {chosen} instead "
                  f"(observed values {sorted(observed)}).")
            self.defaults.mitotic_mask_value = chosen

        # remove 0's and 2's, and fill gaps in the semantic vector.
        mask_value = self.defaults.mitotic_mask_value
        self.tracked.loc[self.tracked.semantic != mask_value, "semantic"] = 1

        # Smoothing runs per particle, over that track's frames in order. Run
        # over the whole column it would bleed across track boundaries, letting
        # the end of one cell's trace close a gap at the start of the next -
        # the rows are merely adjacent in the table, not in the movie.
        semantic = self.tracked.semantic.to_numpy().copy()
        smoothed = np.zeros(len(semantic), dtype=int)
        frames = self.tracked.frame.to_numpy()
        for rows in self.tracked.groupby("particle", sort=False).indices.values():
            rows = rows[np.argsort(frames[rows])]      # frame order within track
            trace = medfilt(semantic[rows], self.defaults.semantic_gap_closing)
            semantic[rows] = trace
            # Turn into Boolean. Compare against the mask value directly; the
            # former (semantic - 1)//99 only happened to work for 100 or 101.
            smoothed[rows] = closing((trace == mask_value).astype(int),
                                     self.defaults.semantic_footprint)
        self.tracked["semantic"] = semantic
        self.tracked["semantic_smoothed"] = smoothed

        # dividing (1) vs non-dividing (0), per track
        is_mitotic = self.tracked.semantic == self.defaults.mitotic_mask_value
        self.tracked["mitotic"] = self.tracked.particle.map(
            is_mitotic.groupby(self.tracked.particle).any()).astype(int)
        return self.tracked

    @classmethod
    def from_analysis_file(cls, analysis_xlsx: Path, cell_type: str = "hela",
                           frame_shape=None, frame_interval=None):
        '''Build an object that can summarize a saved *_analysis.xlsx.

        Only what summarize_data reads is populated - no image stacks are
        loaded - so this is the cheap path for re-summarizing an existing
        analysis after its dead/mitotic labels change, as augment_dead_label
        does. Columns the file predates (semantic_smoothed, mitotic) are
        derived here.

        Inputs:
        analysis_xlsx : path to a *_analysis.xlsx written by measure_signal
        cell_type     : selects the analysis_pars defaults. Only the tracking
                        parameters differ by cell type and summarize_data uses
                        none of them, so this rarely matters here.
        frame_shape   : (rows, cols) at analysis scale, for the border test.
                        Inferred from the bounding boxes when omitted.
        frame_interval: minutes between frames. Read from the acquisition
                        metadata beside the image stacks when omitted.

        The interval is deliberately NOT taken from the file's own parameters
        sheet. Files written before this change carry the hard-coded 10 that
        was never a measurement, and reading it back would quietly restore the
        wrong minimum mitotic duration on exactly the re-summarizing path this
        constructor exists for. The metadata is the source; when it is
        unreadable this raises and asks for the value rather than guessing.

        semantic_smoothed and mitotic are derived only when absent. Deriving
        them from a file that already has them would median-filter an already
        filtered trace, which is not idempotent, so a file that carries them is
        left alone.
        '''
        analysis_xlsx = Path(analysis_xlsx)
        stub = re.search(r"[A-H]([1-9]|[0][1-9]|[1][0-2])_s(\d{2}|\d{1})",
                         analysis_xlsx.name)
        if stub is None:
            raise ValueError(f"cannot parse a well/position stub from "
                             f"{analysis_xlsx.name}; expected something like A12_s2")

        sheets = pd.read_excel(analysis_xlsx, sheet_name=None, index_col=0)
        # the cell table is 'cell_data' in newer files and 'Sheet1' in older ones
        cell = sheets["cell_data" if "cell_data" in sheets else list(sheets)[0]]

        self = cls.__new__(cls)
        self.defaults = analysis_pars(cell_type=cell_type)
        self.paths = {"analysis": analysis_xlsx}
        self.stacks = {}
        self.quality = {}
        self.tracked = cell
        self.cellaap_dir = analysis_xlsx.parent
        self.root_folder = analysis_xlsx.parent.parent
        self.name_stub = stub.group()
        self.expt_name = analysis_xlsx.name.split(self.name_stub)[0]
        if frame_shape is not None:
            self.frame_shape = frame_shape

        # data_dir is where the image stacks and their metadata live: the
        # inference folder's parent, the same directory files() reads.
        self.data_dir = self.cellaap_dir.parent
        if frame_interval is None:
            frame_interval = read_frame_interval(self.data_dir)
        self.defaults.set_frame_interval(frame_interval)

        missing = {"semantic", "particle", "frame", "x", "y", "area"} - set(cell.columns)
        if missing:
            raise ValueError(f"{analysis_xlsx.name} is missing {sorted(missing)}, "
                             f"which summarize_data needs")
        if not {"semantic_smoothed", "mitotic"}.issubset(cell.columns):
            self._label_semantic()
        return self

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
            self.write_analysis_file()

        return self.tracked


    def analysis_file_path(self) -> Path:
        '''Where this position's `*_analysis.xlsx` lives.'''
        return self.cellaap_dir / Path(self.expt_name + self.name_stub + '_analysis.xlsx')


    def write_analysis_file(self):
        '''
        Write `self.tracked` and the run's provenance to `*_analysis.xlsx`.

        Split out of measure_signal because the table is written more than
        once: measure_signal saves it per channel, and anything that adds
        columns afterwards - signal_correction's per-frame corrected signal -
        has to put them on disk before summarize_data runs. Rewrites the whole
        workbook each time, so the sheets never hold a stale mix of columns.
        '''
        with pd.ExcelWriter(self.analysis_file_path()) as writer:
            self.tracked.to_excel(writer, sheet_name='cell_data')
            # There are no scalars - so turn into list; transform
            pd.DataFrame([self.paths]).T.to_excel(writer,   sheet_name='file_data')
            pd.DataFrame([self.defaults.__dict__]).T.to_excel(writer,sheet_name='parameters')
        return self.analysis_file_path()
    
    
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
                if e - s >= self.defaults.require_frame_interval()]

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
        1 - death_proba_threshold after the run. Without this guard such a cell
        is killed on the first frame of its mitosis and reports a zero-length
        one, even when it goes on to divide again later.
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

    def _filter_spurious_tracks(self, summary, movie_last, censored):
        '''Drop tracks whose reported mitotic duration cannot be trusted.

        Runs on the finished summary rather than inside the per-track loop,
        because the calibration needs the durations the loop produces. Adds
        two columns to every row, kept whether or not the filter is enabled:

        obs_window_after_entry  frames from mitotic entry to the END OF THE
                                TRACK - how long the cell could have been
                                watched in mitosis, irrespective of what it
                                did. It is fixed before the outcome is known,
                                which is what makes it safe to filter on: it
                                is exposure time, not a measurement.
        duration_censored       the track ends while the cell is still mitotic,
                                somewhere other than the last frame of the
                                movie, so its duration is a LOWER BOUND. This
                                is a flag only. Excluding on it costs a fifth
                                of the data and does not improve the
                                distribution (median moves 37.0 -> 37.2 h and
                                the short-mitosis rate gets slightly worse),
                                because a censored duration is still a real
                                mitosis seen for a while - unlike the tracks
                                the window test removes, which were barely
                                seen at all.

        Two tests, each anchored to the scale it is actually about:

        * the window must be at least min_window_factor of the reference
          mitotic duration, itself the median over tracks followed for
          reference_track_fraction of the movie. Anchoring to the data rather
          than to a constant is what lets one rule serve a mitotic arrest and
          an unperturbed control - on the datasets this was built against the
          same factor yields a 155-frame threshold for the first and 7 for the
          second.
        * the track must start within max_track_start_fraction of the movie.

        Returns the filtered summary. The reference set is reported in
        self.quality["track_filter"] so a later reader can see what the
        thresholds were calibrated against, which matters because they move
        with the data.
        '''
        n_frames = int(movie_last) + 1
        summary = summary.copy()
        summary["obs_window_after_entry"] = (
            summary.track_start_frame + summary.track_length - 1
            - summary.mitotic_start_frame + 1)
        summary["duration_censored"] = summary.particle.map(censored).fillna(False)
        if not self.defaults.exclude_short_window_tracks:
            return summary          # columns still added, nothing excluded

        reference = summary[summary.track_length
                            >= self.defaults.reference_track_fraction * n_frames]
        if self.defaults.reference_duration_frames is not None:
            duration = float(self.defaults.reference_duration_frames)
        elif len(reference) < 20:
            print(f"only {len(reference)} tracks reach "
                  f"{self.defaults.reference_track_fraction:.0%} of the movie; "
                  f"too few to calibrate the spurious-track filter, so it is "
                  f"skipped for this position. Set "
                  f"defaults.reference_duration_frames to filter it anyway.")
            return summary
        else:
            duration = float(reference.corrected_frames_in_mitosis.median())
        min_window = max(int(np.ceil(self.defaults.min_window_factor * duration)),
                         self.defaults.require_frame_interval())
        max_start = int(self.defaults.max_track_start_fraction * n_frames)

        keep = ((summary.obs_window_after_entry >= min_window)
                & (summary.track_start_frame <= max_start))
        self.quality["track_filter"] = pd.DataFrame([{
            "n_frames": n_frames,
            "n_reference_tracks": len(reference),
            "reference_pinned": self.defaults.reference_duration_frames is not None,
            "reference_duration_frames": duration,
            "min_window_frames": min_window,
            "max_track_start_frame": max_start,
            "n_before": len(summary),
            "n_after": int(keep.sum()),
        }])

        interval = self.defaults.frame_interval
        print(f"spurious-track filter: reference duration "
              f"{duration:.0f} frames ({duration * interval:.0f} min) from "
              f"{len(reference)} tracks followed >= "
              f"{self.defaults.reference_track_fraction * n_frames:.0f} frames")
        print(f"  excluded {int((~keep).sum())} of {len(summary)} tracks "
              f"(observation window < {min_window} frames after mitotic entry, "
              f"or track starting after frame {max_start}); "
              f"{int(keep.sum())} remain")
        return summary[keep].reset_index(drop=True)

    def summarize_data(self, save_flag: True, suffix: str = ""):
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
                                      death and the frames from entry to it.
          after the first episode  -> dead_post_mitosis. The mitosis completed,
                                      so its full duration is reported, plus
                                      the frame of death.

        Death is only visible where the classifier ran: the mitotic frames plus
        defaults.post_peak_frames AFTER each episode. Nothing before mitotic
        entry is ever scored - this model was built to separate mitotic from
        dead cells, both rounded, and has no meaning on an interphase cell.
        The tail is what catches a cell dying on mitotic exit, and is
        deliberately short for the same reason.

        A cell is summarized only if it has at least one episode and (default)
        none of its mitotic-labeled detections sit within defaults.border_margin
        of the frame edge, where the classifier cannot see the whole cell.

        Columns. track_start_frame, mitotic_start_frame and death_frame are
        absolute movie frames and compare directly; the rest of the frame
        counts are durations.

        sem_frames_in_mitosis        every mitotic-labeled frame in the track.
        corrected_frames_in_mitosis  of those, the ones before the death call:
                                     the time the cell spent in mitosis while
                                     still alive, and the window the channel
                                     means are taken over. Equal to
                                     sem_frames_in_mitosis when the cell never
                                     dies, and the difference between the two
                                     is the time it lay dead but still labeled
                                     mitotic.
        frames_to_death              frames from mitotic entry to the death
                                     call; NaN if the cell never dies. It can
                                     exceed corrected_frames_in_mitosis, since
                                     a post-mitotic death happens after the
                                     cell has left mitosis.

        Fluorescence is averaged over the corrected window only, so nothing is
        measured from a cell already called dead.

        A track is excluded outright when

          * it is already mitotic on its own first frame - no entry was
            observed, and the segmentation's habit of labeling anaphase mitotic
            means most of these are a daughter cell trackpy picked up
            mid-division;
          * it begins after defaults.late_track_start_fraction of the movie AND
            reaches mitosis within defaults.early_mitosis_frames of its start -
            the same daughter-cell artifact, for the ones whose first frame or
            two are not yet labeled;
          * it is still mitotic on the last frame of the movie, where the
            acquisition cut the exit off;
          * a mitotic-labeled detection lies within defaults.border_margin of
            the frame edge, where the classifier cannot see a full crop box.

        Note what never reaches this function: tracks shorter than
        defaults.min_track_length are removed by trackpy in track_centroids,
        so a cell that rounds up, dies and loses its track inside that window
        is gone before any of this runs.

        Inputs -
        save_flag : whether to export the data as an xlsx file
        suffix    : appended to the summary file name, before the extension,
                    to match a suffixed analysis file (augment_dead_label)

        Outputs -
        None
        '''
        # Fail here rather than inside the per-track loop: every duration this
        # function reports is scaled by the frame interval, and the episode
        # test below depends on it.
        self.defaults.require_frame_interval()

        # Select only those tracks where mitosis was observed
        idlist    = sorted(set(self.tracked[self.tracked.mitotic==1].particle))

        # Border test: any mitotic-labeled detection closer than the margin
        # means the classifier could not score the cell there.
        margin = self.defaults.border_margin
        rows, cols = self._analysis_frame_shape()
        if margin and self.defaults.exclude_border_tracks:
            mit_rows = self.tracked[self.tracked.semantic_smoothed == 1]
            near = ((mit_rows.x < margin) | (mit_rows.x > rows - margin) |
                    (mit_rows.y < margin) | (mit_rows.y > cols - margin))
            near_border_ids = set(mit_rows.loc[near, 'particle'].unique())
            dropped = len(near_border_ids & set(idlist))
            idlist = [i for i in idlist if i not in near_border_ids]
            print(f'excluded {dropped} tracks with mitotic detections within '
                  f'{margin} px of the frame edge ({rows}x{cols})')

        # Two ways a track can carry an episode that was never fully observed,
        # and they need different tests.
        #
        # Starting mitotic disqualifies a track wherever in the movie it
        # begins. The segmentation labels anaphase mitotic, so when a cell
        # divides trackpy commonly opens a fresh particle on a daughter that is
        # still carrying the mitotic label: its "mitosis" is the tail of the
        # mother's division, with no entry of its own.
        #
        # Ending mitotic only disqualifies a track that runs to the LAST frame
        # of the movie, where the acquisition cut the episode short. A track
        # that merely stops mid-movie is trackpy losing the cell, which says
        # nothing about the mitosis and would cost most of the data to exclude.
        ordered = self.tracked.sort_values(['particle', 'frame'])
        movie_last = self.tracked.frame.max()
        edge = ordered.groupby('particle').agg(
            last_frame=('frame', 'last'),
            first_sem=('semantic_smoothed', 'first'),
            last_sem=('semantic_smoothed', 'last'))
        starts_mitotic = set(edge.index[edge.first_sem == 1])
        ends_clipped = set(edge.index[(edge.last_frame == movie_last)
                                      & (edge.last_sem == 1)])
        # Same test one frame short of the end: the track stops while the cell
        # is still mitotic, so the duration is a lower bound. Reported as a
        # column rather than excluded - see _filter_spurious_tracks.
        censored = ((edge.last_sem == 1) & (edge.last_frame < movie_last))
        for ids, why in ((starts_mitotic, 'starting in mitosis (no entry '
                                          'observed; usually a daughter cell '
                                          'picked up mid-division)'),
                         (ends_clipped, 'still mitotic on the last frame of '
                                        'the movie (no exit observed)')):
            dropped = len(ids & set(idlist))
            if dropped:
                idlist = [i for i in idlist if i not in ids]
                print(f'excluded {dropped} tracks {why}')


        # A list to store the number of peaks
        # Multiple peaks will reveal either tracking errors or segmentation issues
        peaks_per_track = np.zeros(len(idlist)) 
        # Fluctuations in the mask size will indicate segmentation quality
        cell_area_std  = np.zeros_like(peaks_per_track)

        corrected_frames_in_mitosis = [] # mitotic frames before the death call
        sem_frames_in_mitosis = []  # every mitotic-labeled frame in the track
        frames_to_death  = [] # frames from mitotic entry to the death call
        track_start_frame   = [] # movie frame the track begins on
        mitotic_start_frame = [] # movie frame of mitotic entry
        particle         = []
        track_length     = []
        channels         = []
        dead_cell_score  = [] # keep track of "dead" flags
        fate_label       = [] # mitotic_survived / dead_in_mitosis / dead_post_mitosis
        death_frame      = [] # movie frame of death (NaN if the cell never dies)
        n_peaks          = [] # mitotic episodes in the track
        n_scored         = [] # frames the classifier actually scored
        n_interphase_death = 0 # died without ever having a mitosis; excluded
        n_late_daughter    = 0 # late-starting track that is mitotic at once

        # A track starting after this frame is late enough in the movie to be a
        # daughter picked up mid-division rather than a cell followed from the
        # start; frames are 0-indexed, hence the +1 for the movie length.
        late_start_after = (self.defaults.late_track_start_fraction
                            * (int(movie_last) + 1))

        # Check which channels have been measured. If none, return only "mitotic duration"
        # Need to find a better way to code this.

        if "GFP" in self.tracked.columns:
            channels.append('GFP')
        if "Texas Red" in self.tracked.columns:
            channels.append("Texas Red")
        if "Cy5" in self.tracked.columns:
            channels.append("Cy5")
        # Which per-channel columns to average over each track's window. The
        # list is discovered rather than fixed, so a column added upstream -
        # `<ch>_corrected` from signal_correction is the reason this is not
        # hard-coded any more - reaches the summary without another edit here.
        # The second element is what an all-zero trace should report: 0 for a
        # quantity that is added or subtracted, 1 for one that divides.
        CHANNEL_COLUMNS = (("", 0), ("_bkg_corr", 0), ("_int_corr", 1),
                           ("_corrected", 0))
        # NB `column_suffix`, not `suffix`: `suffix` is this function's
        # file-name argument, and a loop over it here would leak into the
        # output name at the bottom.
        measured = {channel: [(cs, empty) for cs, empty in CHANNEL_COLUMNS
                              if f'{channel}{cs}' in self.tracked.columns]
                    for channel in channels}

        signal_storage = {}
        for channel in channels:
            for column_suffix, _ in measured[channel]:
                signal_storage[f'{channel}{column_suffix}'] = []
                signal_storage[f'{channel}{column_suffix}_std'] = []

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

            # The rest of the daughter-cell artifact. Rejecting tracks mitotic
            # on their very first frame misses the ones whose first frame or
            # two are not yet labeled, so a track that both begins late in the
            # movie and reaches mitosis almost immediately is rejected too.
            # Neither condition alone is suspicious.
            if (frames[0] > late_start_after
                    and frames[first_start] - frames[0]
                    < self.defaults.early_mitosis_frames):
                n_late_daughter += 1
                continue

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

            # The measured window: every mitotic-labeled frame before the death
            # call. Nothing is measured from a cell already called dead, so the
            # channel means describe the cell while it was both mitotic and
            # alive, and its size is reported as corrected_frames_in_mitosis -
            # sem_frames_in_mitosis minus whatever fell at or after the death.
            window = sem_raw.astype(int)
            if death is not None:
                window = window * (np.arange(len(frames)) < death)
            corrected_frames_in_mitosis.append(int(window.sum()))
            sem_frames_in_mitosis.append(int(sem_raw.sum()))
            frames_to_death.append(int(death - first_start) if death is not None
                                   else np.nan)
            track_start_frame.append(int(frames[0]))
            mitotic_start_frame.append(int(frames[first_start]))
            death_frame.append(int(frames[death]) if death is not None else np.nan)
            fate_label.append(fate)
            n_peaks.append(len(episodes))
            dead_cell_score.append(np.sum(sem_raw*dead_flag))
            n_scored.append(int(np.isfinite(proba).sum()))
            particle.append(id)
            track_length.append(len(frames))

            # A cell that dies on the first frame of its mitosis has an empty
            # averaging window, so the channel means are legitimately NaN;
            # numpy's empty-slice warnings would otherwise flood the log.
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                track_rows_id = self.tracked[self.tracked.particle==id]
                for channel in channels:
                    for column_suffix, empty in measured[channel]:
                        column = f'{channel}{column_suffix}'
                        mean, std = window_stats(
                            window, track_rows_id[column].to_numpy(), empty)
                        signal_storage[column].append(mean)
                        signal_storage[f'{column}_std'].append(std)


        n_dead = sum(1 for f in fate_label if f != 'mitotic_survived')
        print(f'summarized {len(particle)} tracks; {n_dead} died during or '
              f'after mitosis (P(dead) > {self.defaults.death_proba_threshold} '
              f'for {self.defaults.death_run_frames} consecutive frames)')
        if n_late_daughter:
            print(f'excluded {n_late_daughter} tracks starting after frame '
                  f'{late_start_after:.0f} that reached mitosis within '
                  f'{self.defaults.early_mitosis_frames} frames (daughter '
                  f'cells picked up mid-division)')
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
        # Column order is deliberate: identity, then when the track and the
        # mitosis start, then the durations, then the classifier's evidence and
        # verdict, with the fluorescence columns last.
        other_storage = {
                        "particle"                    : particle,
                        "track_length"                : track_length,
                        "n_peaks"                     : n_peaks,
                        "track_start_frame"           : track_start_frame,
                        "mitotic_start_frame"         : mitotic_start_frame,
                        "frames_to_death"             : frames_to_death,
                        "sem_frames_in_mitosis"       : sem_frames_in_mitosis,
                        "corrected_frames_in_mitosis" : corrected_frames_in_mitosis,
                        "n_scored"                    : n_scored,
                        "dead_cell_score"             : dead_cell_score,
                        "fate_label"                  : fate_label,
                        "death_frame"                 : death_frame,
                        }
        
        summary_storage = other_storage | signal_storage
        self.summaryDF = pd.DataFrame(summary_storage)

        # The window and start tests run last: the first needs the durations
        # computed above to calibrate itself against, and both add columns
        # that describe rows the loop has already built.
        self.summaryDF = self._filter_spurious_tracks(
            self.summaryDF, movie_last, censored)

        if save_flag:
            out = self.cellaap_dir / Path(
                self.expt_name + self.name_stub + "_summary" + suffix + ".xlsx")
            with pd.ExcelWriter(out) as writer:
                self.summaryDF.to_excel(writer,sheet_name = "Summary", index=False)
                pd.DataFrame([self.paths]).T.to_excel(writer, sheet_name='file_data')
                pd.DataFrame([self.defaults.__dict__]).T.to_excel(writer,sheet_name='parameters')
                self.quality["cell_area_std"].join(self.quality["peaks_per_track"]).to_excel(writer, 
                                                                           sheet_name="quality", 
                                                                           index=False)

        return self.summaryDF
    