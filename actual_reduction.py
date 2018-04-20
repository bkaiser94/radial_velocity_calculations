"""
This should read in the raw image files for the spectra and the iron lamp that corresponds
"""

#I need something to indicate that this is supposed to be a standard star or the actual target in order to fix the 


import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
import cosmics
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u


c= 2.998E5  #km/s
zerolistname= 'listZero'
#flatlistname = 'listFlat'
#speclistname= 'listSpec'

#output_filename = 'radial_velocities.txt'


linefilename = 'JJ_FeAr_lines.txt'
#linefilename = 'FeAr_3650to5250_lines_GOODMAN.txt'
#flatlist= np.genfromtxt(flatlistname,dtype = 'str' )
zerolist = np.genfromtxt(zerolistname, dtype ='str')
#speclist= np.genfromtxt(speclistname, dtype = 'str')

#########3
slit_ystart = 1   #The beginning of the image that has light from outside
slit_yend= 199     #The end of the image with same
trace_xstart = 9
trace_xend = 2055
#bkg_width= 10   #How many pixel rows should be sampled on each edge of the slit
trace_band_mid= 105   #y-pixel that's about the center of the bulge of the galaxy
trace_band_width = 20 #pixel width to determine the centroid of the galaxy
poly_degree= 3
core_sides=  5
bkg_width= core_sides
bkg_shift= 50


lamp_sigma_guess = 2
line_search_width= 3
balmer_sigma_guess= 14
lamp_p0 = [100, 500,  lamp_sigma_guess, 0]
#balmer_p0= [-1, 500, balmer_sigma_guess,balmer_line_sides[0],0]
balmer_p0= [-100, 500, balmer_sigma_guess,0]


fear_array= np.genfromtxt(linefilename, names = True)
line_x_checks = np.copy(fear_array['Pixel']) +90
lamp_lines = np.copy(fear_array['User'])
line_sides = np.ones(line_x_checks.shape[0])*line_search_width


