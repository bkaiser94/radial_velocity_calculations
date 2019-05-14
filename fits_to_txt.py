"""
Created by Ben Kaiser 2019-05-14 (UNC-Chapel Hill)

Read in spectra in the FITS format typical of my own 'radial_velocity_calcuations/' spectra and then output
them as .txt files that can be used by PyHammer...hopefully.

The output filenames should be identical in everyway except that they'll be .txt files instead



"""
from __future__ import print_function
import numpy as np
#import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
#from astropy.time import Time
#from astropy import coordinates as coords
#from astropy import units as u
#from astropy import constants as const
#from astropy import convolution as conv
#import scipy.interpolate as scinterp
#import time
#start = time.time()



#print start
import spec_plot_tools as spt


input_string= 'avg_fwctb*fits'

input_filenames= glob(input_string)
