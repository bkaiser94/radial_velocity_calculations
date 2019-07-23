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
import scipy.optimize as sciop

import cal_params as cp
import spec_plot_tools as spt



linelist_file= ''
wave_soln_file= ''

wave_sol_binning= 2 #binning of the wavelength solution


trace_width=10
trace_mid=100


############################

def get_binning(header):
    """
    get the dispersion direction binning of pixels to figure out how to map a polynomial from the wavelength calibrations.
    """
    binning= header['CCDSUM'].split(' ')
    xbinning= binning[0]
    return xbinning



def wavelength_to_pixel(lambda_val, in_wave_coeffs, lamp_poly_degree=5, bounds=[0,2100]):
    """
    input wave_coeffs should already have an offset subtracted from the x-values everywhere.... you can't really 
    do that...
    
    bounds: the pixel boundaries at which the lambda_val might be located (it should be 
    """
    wave_coeffs= np.copy(in_wave_coeffs)
    wave_coeffs[-1]= wave_coeffs[-1]-lambda_val
    
    def func_to_solve(x):
        if lamp_poly_degree==5:
            return wave_coeffs[0]*x**5+ wave_coeffs[1]*x**4 +wave_coeffs[2]*x**3+ wave_coeffs[3]*x**2+wave_coeffs[4]*x + wave_coeffs[5]
        else:
            print("don't have function to solve for inversion of wavelengths for that lamp_poly_degree:", lamp_poly_degree)
            return np.polyval(wave_coeffs, x)
    #plt.plot(np.polyval(wave_coeffs, np.linspace(0,2000,2000)), label='changed wave_coeffs')
    #plt.plot(np.polyval(in_wave_coeffs, np.linspace(0,2000,2000)), label='og wave_coeffs')
    #plt.plot(func_to_solve(np.linspace(0,2000,2000)),label='func_to_solve')
    #plt.legend(loc='best')
    #plt.show()
    pixel= sciop.brentq(func_to_solve, bounds[0],bounds[1])
    return pixel




######################

