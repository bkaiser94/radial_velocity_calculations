"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-08-28


take a "flux-calibrated" standard spectrum and correct any wavelength offset between the balmer line of the observed spectrum and the standard model


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
from astropy.table import Table, Column

import get_cal_params as gcp
import cal_params as cp
import spec_plot_tools as spt

wave_shift_range= [-20, 20]
wave_targ_center = 6562 #approximate wavelength center to be targeted
search_width= 40 #in angstroms


use_fnu= False

#gauss_p0= [-3.e-14, 6562. , 2., 1e-14]
gauss_p0= [-4000., 6562. , 2., 1e4]



##standard_name = "GD108"
#standard_name = 'Feige67'
#standard_name = 'LTT6248'
standard_name='EG274'
#standard_name = 'GD153'
#standard_name= 'LTT3218'
#standard_name='Feige110'
#standard_name= 'LTT7987'


#observed_file='avg_fwctb.LTT7987_400m2.fits'
#observed_file='avg_fwctb.LTT7987_400m1.fits'


#observed_file='avg_fwctb.EG274_400m1.fits'
observed_file='avg_fwctb.EG274_400m2.fits'
#observed_file='fwctb.0291_EG274_400m2.fits'

#observed_file='avg_fwctb.Feige110_400m2.fits'
#observed_file='avg_fwctb.Feige110_400m1.fits'

#observed_file='avg_fwctb.GD153_400m2.fits'
#observed_file='avg_fwctb.GD153_400m1.fits'


def get_star_info(starname):
    standard_dict= cp.standard_dict[starname.lower()]
    standard_dict['filename']=cp.standard_dir+standard_dict['filename']
    return standard_dict

def get_output_header(header):
    airmass= header['AIRMASS']
    obs_time = header['OPENTIME']
    obs_date = header['OPENDATE']
    obs_time = obs_date+'T'+obs_time
    obs_time = Time(obs_time, format = 'isot', scale = 'utc').mjd
    output_header_list= ['Airmass'+header_char+str(airmass), 'MJD'+header_char+str(obs_time)]
    for in_header, out_header in zip(cp.in_headers, cp.out_headers):
        value= header[in_header]
        new_entry= out_header+header_char+str(value)
        output_header_list.append(new_entry)
    output_header= header_delim.join(output_header_list)
    return output_header

def degrade_model(model_vals, obs_vals, header):
    """
    Convolve model spectrum with the seeing and rebin it (using flux-conservative method) to the pixel-scale of
    the observation.
    
    INPUTS:
        model_vals - [wavelengths, fluxes, wavelength bin widths]
        obs_vals - [wavelengths, fluxes, wavelength bin widths]
        header - observation header, which should include the pixel width, slit width, seeing, etc.
        
    OUTPUTS:
        model_spec 
    
    """
    rebinned_model_spec= spt.rebin_generic_spec(model_vals[:2], model_vals[2], obs_vals[0], obs_vals[2])
    slit_width = spt.get_slit_width(header)
    rebinned_model_spec= mm.convolve_model_new(rebinned_model_spec, header, slit_width=slit_width)
    output_spec= rebinned_model_spec
    return output_spec


def gaussian_curve(x, a, x0, sigma,b):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+b


def fit_gaussian_curve(x_pixels, light_values, p0_list, search_width, plot_all = False, bounds = (-np.inf, np.inf), fixed_width=True):
    """
    Those bounds are the default for scipy.optimize.curve_fit(), so now changing them changes the bounds
    """
    cut_region = np.where(x_pixels> (p0_list[1]-search_width ))
    high_x_pixels= np.copy(x_pixels[cut_region])
    high_light_values= np.copy(light_values[cut_region])
    upper_cut = np.where(high_x_pixels < (p0_list[1]+search_width))
    cut_x_pixels = high_x_pixels[upper_cut]
    cut_light_values= high_light_values[upper_cut]
    #print('p0_list it lets you use:', p0_list)
    popt, pcov = sciop.curve_fit(gaussian_curve, cut_x_pixels, cut_light_values, p0= p0_list, bounds = bounds)
    
    return popt, pcov


def rescale_flux(stand_flux1, standard_info):
    """
    Make sure the flux is in the correct units for your purposes, which should be 1e-16
    
    """
    standard_type= standard_info['filename'].split('/')[-2]
    print('standard_type:', standard_type)
    if standard_type == 'xshooter_standards':
        print('should be resetting fluxes')
        stand_flux1= stand_flux1*1e16 #converts to 10**-16 flux vals hopefully
    else:
        stand_flux1= stand_flux1
    return stand_flux1



####################################



obs_fits = fits.open(observed_file)
header = fits.getheader(observed_file)
obs_waves1= obs_fits[0].data
obs_flux1 = obs_fits[1].data
obs_spec= np.vstack([obs_waves1, obs_flux1])
obs_dlambda= obs_fits[4].data
airmass = header['AIRMASS']
obs_time = header['OPENTIME']
obs_date = header['OPENDATE']
obs_time = obs_date+'T'+obs_time
obs_time = Time(obs_time, format = 'isot', scale = 'utc').mjd
exptime = header['EXPTIME']

obs_spec= np.vstack([obs_waves1, obs_flux1])

    

setup_dict= gcp.get_cal_params(header)
setup_name=setup_dict['setupname']
sens_fit_method= cp.flux_cal_dict['sens_fit_method'][setup_name]

#standard_file = standard_directory+standard_file
standard_info = get_star_info(standard_name)
print type(standard_info['balmer_masks'])
print type(standard_info['other_masks'])
wavelength_masks=standard_info['balmer_masks']+standard_info['other_masks']
print "wavelength_masks:", wavelength_masks

stand_array = np.genfromtxt(glob(standard_info['filename'])[0]).T
#output_filename= standard_dict['sens_filename']

stand_waves1 = stand_array[0]


stand_flux1 = stand_array[1]  #ergs/cm/cm/s/A (That's exactly how it's written in the README for X-shooter)

stand_flux1= rescale_flux(stand_flux1, standard_info)

unshifted_waves= np.copy(stand_waves1)
unshifted_flux= np.copy(stand_flux1)
#plt.plot(stand_waves1, stand_flux1, label='unshifted')
stand_waves1= spt.barycentric_vel_uncorr(header, stand_waves1, sys_vel=standard_info['sys_vel'])
#stand_waves1= spt.barycentric_vel_uncorr(header, stand_waves1, sys_vel=0)

    


if use_fnu:
    #stand_flux1= rescale_flux(stand_flux1, standard_info)
    plt.plot(stand_waves1, stand_flux1, label='flambda')
    stand_spec= np.vstack([stand_waves1, stand_flux1])
    stand_spec=spt.flambda_to_fnu(stand_spec)
    stand_flux1=stand_spec[1]
    stand_waves1=stand_spec[0]
    plt.plot(stand_waves1, stand_flux1, label='fnu')
    plt.legend()
    plt.show()
else:
    pass



standard_type= standard_info['filename'].split('/')[-2]
print('standard_type:', standard_type)
if standard_type == 'xshooter_standards':
    print('no bin widths provided for model because X-Shooter')
    model_bin_widths= np.copy(np.roll(stand_waves1, -1) - stand_waves1)
    print(stand_waves1.shape, stand_flux1.shape, model_bin_widths.shape)
    stand_waves1= stand_waves1[:-1]
    stand_flux1=stand_flux1[:-1]
    model_bin_widths= model_bin_widths[:-1] #need to remove the first two pixels because they're going to be weird... or mayb it's the last two...
    plt.plot(stand_waves1, model_bin_widths)
    plt.title('model bin widths')
    plt.show()

plt.plot(stand_waves1, stand_flux1, label='shifted')
plt.plot(unshifted_waves, unshifted_flux, label='from file')
plt.plot(obs_spec[0], obs_spec[1], label='observation')
plt.legend()
plt.show()

model_fit_gauss, discard_pcov= fit_gaussian_curve(stand_waves1, stand_flux1, gauss_p0, search_width)
print(model_fit_gauss)
print('central wavelength of gaussian for model: ', model_fit_gauss[1])

gauss_values= gaussian_curve(stand_waves1, model_fit_gauss[0], model_fit_gauss[1], model_fit_gauss[2], model_fit_gauss[3])
plt.plot(stand_waves1, gauss_values, label='gaussian')
plt.plot(stand_waves1, stand_flux1, label='model')
plt.legend()
plt.show()




obs_fit_gauss, discard_pcov= fit_gaussian_curve(obs_spec[0], obs_spec[1], gauss_p0, search_width)
print(obs_fit_gauss)
print('central wavelength of gaussian for obs: ', obs_fit_gauss[1])

gauss_values= gaussian_curve(obs_spec[0], obs_fit_gauss[0], obs_fit_gauss[1], obs_fit_gauss[2], obs_fit_gauss[3])
plt.plot(obs_spec[0], gauss_values, label='gaussian')
plt.plot(obs_spec[0], obs_spec[1], label='obs')
plt.legend()
plt.show()


