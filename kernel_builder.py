
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const
from astropy import convolution as conv


import spec_plot_tools as spt

slit_width = 1.0 #arcseconds
wave_poly_deg =5


    
#filename = glob('fwctb')[0]
def get_waves_to_pix_poly(wavelengths):
    indices = np.indices(wavelengths)
    pixel_values = indices[0]
    waves_to_pix_poly = np.polyfit(wavelengths, pixel_values, wave_poly_deg)
    pixel_values_calc = np.polyval(waves_to_pix_poly, wavelengths)
    plt.title('Residuals of Wavelength to pixel polynomial fit')
    plt.plot(wavelengths, pixel_values- pixel_values_calc, linestyle  = 'none', marker = 'o')
    
    
    
    return waves_to_pix_poly


def convolve_model(model_spec, target_spec, header):
    """
    receive the fits file input of the target because you need a number of things from the header.
    """
    angs_per_pix = target_spec[0][target_spec[0].shape[0]/2]-target_spec[0][target_spec[0].shape[0]/2-1] #angstroms per pixel at the midpoint of the spectrum
    pixel_scale = float(header['PIX_SCAL'])
    see_sig = float(header['SEE_SIG']) #sigma value of gaussian fit to do the kernel construction
    slit_width = slit_width/pixel_scale #eventually should read this out of the
    
    
    angs_per_ind = model_spec[0]
    
    
    see_kernel = conv.Gaussian1DKernel(see_sig)
    #slit_kernel = conv.
    #make the gaussian kernel
    final_kernel =1
    
    
    return final_kernel
