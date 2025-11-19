"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-11-19

This should load the CSV that contains all of the points to use to interpolate to get the 
correction for the dip. It is possible that this will change over time.


Copilot help with pathing

"""
# read_relative_numpy.py
from pathlib import Path
import numpy as np
import sys
import scipy.interpolate as scinterp

import matplotlib.pyplot as plt


input_file='dip_correction_points.csv'
# Get the directory that contains this script
try:
    script_dir = Path(__file__).resolve().parent
except NameError:
    # __file__ may not exist in some interactive environments; fall back to argv[0]
    script_dir = Path(sys.argv[0]).resolve().parent

file_path = script_dir / input_file   # change to your filename

# Example: text file with whitespace-separated columns
input_array = np.genfromtxt(file_path)
print(input_array.shape)

corrector=scinterp.Akima1DInterpolator(input_array[0],input_array[1])

#waves=np.arange(3700,8000., 1.)
#plt.plot(waves,corrector(waves))
#plt.show()





