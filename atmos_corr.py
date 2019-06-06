"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-06-06, D-Day

This should be a test-bed for the more complicated atmospheric correction processes that I'll be undertaking
I'm not sure if this is going to evolve into a separate script for atmospheric corrections or if it will be copied and 
pasted into flux_calibration.py and calibrate_flux.py; I just don't want to clutter those even further with
unnecessary plotting as I'm about to do here.

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
