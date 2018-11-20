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

#output_filename='combined_PSRJ1431m4715_20181017.fits'
#output_filename='combined_PSRJ1431m4715_20181031_2spec.fits'
#output_filename='combined_PSRJ1431m4715_20181031_cube_interp.fits'
#output_filename='combined_PSRJ1431m4715_20181105.fits'
output_filename='combined_PSRJ1431m4715_20181114_new.fits'

#output_filename = "wcmtb.GD108930blue.fits"
#listfile= 'listGD108'

#output_filename = 'wcmtb.feige67930blue.fits'
#listfile = 'listFeige67'

#output_filename  = 'wcmtb.ltt6248930blue.fits'
#listfile= 'listLtt6248'

#output_filename = 'wcmtb.eg274930blue.fits'
#listfile= 'listEg274'


target_list = np.genfromtxt(listfile, dtype ='str')

wavelength_range= [3600, 5290]
rms_range= [4600,4650]


flux_list = []
noise_list = []
#first=fits.open(target_list[-4])
first=fits.open(target_list[0])
#header = fits.getheader(target_file)
standard_waves= first[0].data
standard_flux = first[1].data
first_spec = np.vstack([standard_waves, standard_flux])

trimmed_first = spt.trim_spec(first_spec, wavelength_range[0], wavelength_range[1])
trimmed_waves = trimmed_first[0] #These are the wavelength values that all of the other spectra should be interpolated to

def deal_with_bad_vals(input_spec, noise_spec):
    bad_noise = np.isnan(noise_spec[1])
    print "sum bad_noise:", np.sum(bad_noise)
    print "Bad Pixel:" ,input_spec[0][bad_noise], input_spec[1][bad_noise], noise_spec[1][bad_noise]
    #noise_spec[1][bad_noise]= 1.
    input_spec[1][bad_noise]= 0. #setting the flux of the pixels that have messed up things equal to zero
    noise_spec[1][bad_noise]=1.
    other_noise = np.isinf(noise_spec[1])
    print "sum other_noise:", np.sum(other_noise)
    #noise_spec[1][other_noise]= 1.
    input_spec[1][other_noise]=0. #setting the flux to zero
    return input_spec, noise_spec

for target_file in target_list:
    target_spec, header, target_noise = spt.retrieve_spec(target_file)
    #interpolator = scinterp.CubicSpline(target_spec[0], target_spec[1])
    #interp_flux = interpolator(trimmed_waves)
    #interp_flux= np.interp(trimmed_waves, target_spec[0], target_spec[1])
    print "max",  np.max(target_noise[1])
    #bad_noise = np.isnan(target_noise[1])
    #target_noise[1][bad_noise]= 1e6
    #other_noise = np.isinf(target_noise[1])
    #target_noise[1][other_noise]= 1e6
    mask_list=[]
    rms_target= np.copy(spt.clean_spectrum(target_spec, rms_range[0],rms_range[1],mask_list))
    rms_noise_spec= np.copy(spt.clean_spectrum(target_noise, rms_range[0], rms_range[1],mask_list))
    mean_noise = np.nanmean(rms_noise_spec[1])
    rms_about_mean= np.sqrt(np.nanmean((rms_target[1]-np.nanmean(rms_target[1]))**2))
    rms_about_median= np.sqrt(np.nanmean((rms_target[1]-np.nanmedian(rms_target[1]))**2))
    print "============="
    print "Pre-interpolation"
    print "Mean sigma of", rms_range,":", mean_noise
    print "RMS in", rms_range,":", rms_about_mean, "compared to mean:", np.nanmean(rms_target[1])
    print "RMS in",rms_range,":", rms_about_median, "compared to median:", np.nanmedian(rms_target[1])
    print "max noise in", rms_range, ":", np.nanmax(rms_noise_spec[1])
    print "statistical std dev:", np.std(rms_target[1])
    print "-------------------"
    print "Post-interpolation"
    target_spec, target_noise= deal_with_bad_vals(target_spec, target_noise)
    interpolator_noise = scinterp.CubicSpline(target_noise[0], target_noise[1])
    interp_noise = interpolator_noise(trimmed_waves)
    
    print "---------------"
    #interp_flux= np.interp(trimmed_waves, target_spec[0], target_spec[1])
    interpolator= scinterp.CubicSpline(target_spec[0], target_spec[1])
    interp_flux = interpolator(trimmed_waves)
    #interp_flux= np.interp(trimmed_waves, target_spec[0], target_spec[1])
    #interp_noise= np.interp(trimmed_waves, target_noise[0], target_noise[1])
    rms_target=np.copy(np.vstack([trimmed_waves,interp_flux]))
    rms_noise_spec= np.copy(np.vstack([trimmed_waves, interp_noise]))
    rms_target= spt.clean_spectrum(rms_target, rms_range[0], rms_range[1], mask_list)
    rms_noise_spec= spt.clean_spectrum(rms_noise_spec, rms_range[0],rms_range[1], mask_list)
    mean_noise = np.nanmean(rms_noise_spec[1])
    rms_about_mean= np.sqrt(np.nanmean((rms_target[1]-np.nanmean(rms_target[1]))**2))
    rms_about_median= np.sqrt(np.nanmean((rms_target[1]-np.nanmedian(rms_target[1]))**2))
    print "Mean sigma of", rms_range,":", mean_noise
    print "RMS in", rms_range,":", rms_about_mean, "compared to mean:", np.nanmean(rms_target[1])
    print "RMS in",rms_range,":", rms_about_median, "compared to median:", np.nanmedian(rms_target[1])
    print "max noise in", rms_range, ":", np.nanmax(rms_noise_spec[1])
    print "statistical std dev:", np.std(rms_target[1])
    print "new max", np.max(target_noise[1])
    print "================"
    flux_list.append(interp_flux)
    noise_list.append(interp_noise)
    #noise_list.append(np.copy(target_noise[1]))
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
print "noise_array", noise_array
print "avg_comb", avg_comb
#noise_comb = np.nanmedian(noise_list, axis= 0) #median combining the noises for each spectrum
#noise_comb = np.sqrt(np.sum(noise_array**2, axis=0)/noise_array.shape[0]) #sum in quadrature of the noises for each spectrum
noise_comb = np.sqrt(np.sum(noise_array**2, axis=0))/noise_array.shape[0] #sum in quadrature of the noises for each spectrum'
print "noise_array.shape[0]:", noise_array.shape[0]
for value in noise_array:
    print "sigma", value[0]
print "combined sigma in quad for that:", noise_comb[0]
print "np.nanmean(noise_array, axis=0)[0]", np.nanmean(noise_array, axis=0)[0]
print "array**2", (noise_array**2)[1][2]
print "array[0][0]**2", noise_array[1][2]**2
print "noise_comb", noise_comb
#noise_comb= noise_comb/noise_array.shape[0] #per David, this should be done inside the square root
#weights = 1./noise_array
weights = 1./(noise_array**2) #weight should be 1 over square of sigma
#for weight in weights:
    #print "weight", weight
weight_comb= np.average(flux_array, axis= 0, weights = weights)
#combined_noise= np.sqrt(np.average(noise_array**2, axis=0))
plt.plot(trimmed_waves, med_comb, label = 'median')
plt.plot(trimmed_waves, avg_comb, label ='mean')
plt.plot(trimmed_waves, weight_comb, label= 'weighted')

plt.xlabel(r'Wavelength $(\AA)$')
plt.ylabel(r'Flux (ergs/cm/cm/s/A 10**-16)')
plt.legend()
plt.show()

plt.plot(trimmed_waves, weight_comb/avg_comb, label= "weight_comb/avg_comb")
plt.xlabel(r'Wavelength $(\AA)$')
plt.ylabel(r'Flux (ergs/cm/cm/s/A 10**-16)')
plt.legend()
plt.show()

plt.plot(trimmed_waves, noise_comb, label = 'quad')
#plt.plot(trimmed_waves, combined_noise, label = 'quad')
plt.plot(trimmed_waves, np.median(noise_array, axis=0), label='median')
plt.legend()
plt.show()

hdu=fits.PrimaryHDU(trimmed_waves, header= header)
hdu1= fits.ImageHDU(weight_comb)
hdu2= fits.ImageHDU(np.ones(trimmed_waves.shape))
hdu3 = fits.ImageHDU(noise_comb/weight_comb)
hdulist= fits.HDUList([hdu,hdu1,hdu2,hdu3])
#hdulist.writeto('combined_PSRJ1431m4715_new.fits', overwrite= True)
#hdulist.writeto('combined_PSRJ1431m4715_new.fits', overwrite= True)
hdulist.writeto(output_filename, overwrite= True)
