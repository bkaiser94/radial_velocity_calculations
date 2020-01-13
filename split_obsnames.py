"""
Created by Ben Kaiser (UNC-Chapel Hill) 2020-01-13


Basically take some time-series spectroscopic observations and slightly alter the names of the files so that they'll 
be interepreted by avg_spec.py to be different objects


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
import shutil
start = time.time()

import spec_plot_tools as spt

input_filenames= glob("fwctb.*SDSS*.fits")


n_divisions=4

index_num=0

for name in input_filenames:
    big_parts=name.split('.')
    small_parts=big_parts[1].split('_')
    obj_name= small_parts[1]
    obj_name= obj_name+"n"+str(index_num%n_divisions)
    small_parts[1]=obj_name
    big_parts[1]="_".join(small_parts)
    output_name=".".join(big_parts)
    
    print("copying", name, "to", output_name)
    shutil.copyfile(name, output_name)
    index_num+=1
    
