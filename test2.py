"""


"""


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
from astropy.table import Table, Column
import scipy.interpolate as scinterp
import time
start = time.time()

plt.rc('lines',linewidth=0.5)
#plt.rc('font', size =18)

#print start
#import wdatmos
import spec_plot_tools as spt
import cal_params as cp
import plot_spec as ps



fitsfile='ravg_fwctb.EC01578m1743_930_blue_045asec.fits'

textfile=''
