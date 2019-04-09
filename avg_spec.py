"""
Created by Ben Kaiser 2019-03-27 (UNC-Chapel Hill)

@author: Ben Kaiser

Should average together the counts and uncertainties of matching extracted spectra.

"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const
from astropy import convolution as conv
import scipy.interpolate as scinterp
import time
start = time.time()

import spec_plot_tools as spt


filenames= glob('wctb*')
low_index=10
high_index=-5

filename_matches=[]
prev_core=''
filename_set=[]
for filename in filenames:
    core_name= filename[low_index:high_index]
    if core_name == prev_core:
        filename_set.append(filename)
        print('match')
        print(core_name, prev_core)
    else:
        filename_matches.append([filename_set])
        filename_set=[]
    prev_core= core_name
    print('======')
    print(filename)
    print(core_name)
    print('+++++++')
filename_matches.append([filename_set])
for sets in filename_matches:
    print('======')
    print(sets)
    print('+++++++')
