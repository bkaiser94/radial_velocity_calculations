"""
Created 2019-01-29 by Ben Kaiser (UNC-Chapel Hill)

@author: Ben Kaiser

This should just be able to plot a given spectrum and various things about it.

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



#print start
import wdatmos
import spec_plot_tools as spt

slit_width = 1.0 #arcseconds
pixel_scale = 0.3 #arcseconds per pixel_scale
slit_width = slit_width/pixel_scale #slit width in pixels
test_wavelength = 4686
test_width = 40
test_side = test_width/2

#filenames = glob.(sys.argv[1])
filenames= glob('wctb*')
#print(filenames)

def convolve_spectrum(target_spec, header, kernel_type='gaussian'):
    pix_width =3
    see_sig = float(header['SEE_SIG']) #sigma value of gaussian fit to do the 
    fluxes= np.copy(target_spec[1])
    wavelengths = np.copy(target_spec[0])
    if kernel_type=='gaussian':
        see_kernel = conv.Gaussian1DKernel(see_sig, x_size = int(slit_width), mode = 'oversample')
        see_kernel.normalize()
        spec_conv = conv.convolve(fluxes, see_kernel)
    elif kernel_type== 'box':
        pix_kernel = conv.Box1DKernel(width = int(pix_width), mode = 'oversample')
        pix_kernel.normalize()
        spec_conv = conv.convolve(fluxes, pix_kernel)
    else:
        pass
    spec_out = np.vstack([wavelengths, spec_conv])
    return spec_out


def plot_spectrum(spec, filename, header, smooth=False, kernel_type='gaussian'):
    if smooth:
        spec= convolve_spectrum(spec, header, kernel_type=kernel_type)
    else:
        pass
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux')
    plt.title(filename)
    plt.plot(spec[0], spec[1])
    plt.show()
    return

def plot_dwavelength(spec):
    plt.ylabel(r'delta Wavelength ($\AA$)')
    plt.xlabel('Wavelength ($\AA$)')
    plt.title(filename)
    dlambda= spec[0][1:]-spec[0][:-1]
    plt.plot(spec[0][:-1], dlambda)
    plt.show()
    
    #plt.ylabel(r'Wavelength ($\AA$)')
    #plt.xlabel('pixel')
    #plt.title(filename)
    #plt.plot(spec[0])
    #plt.show()
    #return


def plot_SNR(spec, noise, filename):
    center_pixel = np.argmin(np.abs(spec[0]-test_wavelength))
    measured_std = np.std(spec[1][center_pixel-test_side:center_pixel+test_side])
    print("sigma in " + str(test_width) + " pixel range around " + str(test_wavelength)+ " angstroms", measured_std)
    sigma_range = noise[1][center_pixel-test_side:center_pixel+test_side]
    print("Predicted sigmas of " + str(test_width) +" pixel range around " + str(test_wavelength)+ " angstroms", "min:" + str(np.min(sigma_range)), "mean:" + str(np.mean(sigma_range)), "max:" + str(np.max(sigma_range)))
    print("Mean S/N:", np.mean(spec[1]/noise[1]))
    plt.xlabel('Noise')
    plt.title(filename)
    plt.hist(sigma_range)
    plt.show()
    
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Signal/Noise')
    plt.vlines([spec[0][center_pixel-test_side],spec[0][center_pixel+test_side]], np.min(spec[1]/noise[1]), np.max(spec[1]/noise[1]))
    plt.title(filename)
    plt.plot(spec[0], spec[1]/noise[1], color = 'r')
    plt.show()
    return
    return

for filename in filenames:
    target_spec, header, target_noise= spt.retrieve_spec(filename)
    #plot_spectrum(target_spec, filename, header, smooth=True)
    #plot_spectrum(target_spec, filename, header, smooth=True, kernel_type='box')
    plot_spectrum(target_spec, filename, header)
    plot_SNR(target_spec, target_noise, filename)
    plot_dwavelength(target_spec)

    
    
