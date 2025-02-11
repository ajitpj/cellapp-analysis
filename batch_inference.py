import sys, tifffile
from pathlib import Path
import inference as inf
from skimage.io import imread
import numpy as np

model_name = 'HeLa' # can be on of ['HeLa', 'U2OS']
confluency_est = 1800 # can be in the interval (0, 2000]
conf_threshold = .25 # can be in the interval (0, 1)

# folder definition
root_folder = Path('/nfs/turbo/umms-ajitj/anishjv/HeLa/for_inference/')
save_dir = Path('/nfs/turbo/umms-ajitj/anishjv/HeLa/complete/')
filter_str  = '*.tif'

file_list = []
for phs_file_name in root_folder.glob(filter_str):
    file_list.append(phs_file_name)
num_files = len(file_list)

# Following code is modified from Ajit P. Joglekar

def main():

    container = inf.configure(model_name, confluency_est, conf_threshold, save_dir = save_dir)
    for i in np.arange(num_files):
        phs_file = tifffile.TiffFile(file_list[i])
        interval = [0, len(phs_file.pages)-1]
        result = inf.run_inference(container, file_list[i], interval)
        inf.save(container, result)
        print(f"{file_list[i]} written!")
        print(f"{i} out of {num_files} processed")

    return


if __name__ == "__main__":
    main()



