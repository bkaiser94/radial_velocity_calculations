"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-08-05


Random script for testing random things without having to get a whole thing working in the actual reduction 
process.



"""


import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
import cosmics
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const
from astropy.modeling import models as asmodels
from astropy.modeling import fitting as asfitting
from astropy.table import Table, Column

import get_cal_params as gcp
import cal_params as cp
import spec_plot_tools as spt


input_file= 'ctb.0312_EG274_400m2.fits'

i= fits.open(input_file)
header = fits.getheader(input_file)
img_data= np.copy(i[0].data)

spt.rebin_image(img_data, rebin_axis=1, rebin_num=10, plot_all=True)
