"""
Created by Ben Kaiser (UNC - Chapel Hill) 06-30-2018

This should open all of the spectra that are supposed to be combined to produce the final one and then they
should interpolate the spectra to all have the same wavelength values, and the noise and flux should be 
interpolated in that same idea. I believe we would want to interpolate the normalized noise, so that's what my plan 
will be.

I think I'm supposed to do a weighted mean of the values.

I'm also going to have to truncate the spectra to a desired wavelength range since they all cover slightly different
ones due to the radial velocity differences that were corrected.

"""
import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
import astropy.constants as const
from astropy.units import cds
import scipy.stats as scistats
import scipy.optimize as sciop
import scipy.interpolate as scinterp
cds.enable()

import spec_plot_tools as spt


listfile = 'listZFWCTB'
target_list = np.genfromtxt(listfile, dtype ='str')

wavelength_range= [3600, 5290]
flux_list = []
noise_list = []
first=fits.open(target_list[-4])
#header = fits.getheader(target_file)
standard_waves= first[0].data
standard_flux = first[1].data
first_spec = np.vstack([standard_waves, standard_flux])

trimmed_first = spt.trim_spec(first_spec, wavelength_range[0], wavelength_range[1])
trimmed_waves = trimmed_first[0] #These are the wavelength values that all of the other spectra should be interpolated to



for target_file in target_list:
    target_spec, header, target_noise = spt.retrieve_spec(target_file)
    interpolator = scinterp.CubicSpline(target_spec[0], target_spec[1])
    interp_flux = interpolator(trimmed_waves)
    print "max",  np.max(target_noise[1])
    bad_noise = np.isnan(target_noise[1])
    target_noise[1][bad_noise]= 1e6
    other_noise = np.isinf(target_noise[1])
    target_noise[1][other_noise]= 1e6
    interpolator_noise = scinterp.CubicSpline(target_noise[0], target_noise[1])
    interp_noise = interpolator_noise(trimmed_waves)
    flux_list.append(interp_flux)
    noise_list.append(interp_noise)
    #plt.plot(target_spec[0], target_spec[1], label = target_file)
    #plt.plot(trimmed_waves, interp_flux, label = 'interpolated')
    #plt.legend()
    #plt.show()
    
flux_array = np.array(flux_list)
noise_array = np.array(noise_list)
print flux_array.shape
plt.title('median_combined')
med_comb = np.nanmedian(flux_array, axis = 0)
avg_comb = np.nanmean(flux_array, axis=0)
noise_comb = np.nanmedian(noise_list, axis= 0)
weights = 1./noise_array
weight_comb= np.average(flux_array, axis= 0, weights = weights)
combined_noise= np.sqrt(np.average(noise_array**2, axis=0))
plt.plot(trimmed_waves, med_comb, label = 'median')
plt.plot(trimmed_waves, avg_comb, label ='mean')
plt.plot(trimmed_waves, weight_comb, label= 'weighted')

plt.xlabel(r'Wavelength $(\AA)$')
plt.ylabel(r'Flux (ergs/cm/cm/s/A 10**-16)')
plt.legend()
plt.show()

plt.plot(trimmed_waves, noise_comb, label = 'med')
plt.plot(trimmed_waves, combined_noise, label = 'quad')
plt.legend()
plt.show()

hdu=fits.PrimaryHDU(trimmed_waves, header= header)
hdu1= fits.ImageHDU(weight_comb)
hdu2= fits.ImageHDU(np.ones(trimmed_waves.shape))
hdu3 = fits.ImageHDU(noise_comb/weight_comb)
hdulist= fits.HDUList([hdu,hdu1,hdu2,hdu3])
hdulist.writeto('combined_PSRJ1431m4715.fits', overwrite= True)
