"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-07-22


This is supposed to take a line list (I'm using one from NIST at the moment, but there's no reason it has to be
from there other than consistency of headers), and it uses that line list on a lamp image that is wavelength 
calibrated already using existing line lists and a different lamp image (one that doesn't have the new lines to be 
identified or actually might as a check).

This doesn't do any actual polynomial fitting. It should literally produce an intensity plot, and that's it. Well, it 
should also plot the line labels over top of it...





"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import balmer_line_ranges as blr
from astropy import units as u
from astropy import constants as const
from astropy.table import Table
import scipy.interpolate as scinterp

import cal_params as cp
import spec_plot_tools as spt



linelist_file= ''



