"""
using a model-template already fitted to a median-combined group of observations, we now find the radial velocity of each individual observation and zero it out to the rest frame while appending a radial velocity value to the header in addition to the model that was used for the 
"""

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
from astropy import constants as const
from astropy import convolution as conv
import scipy.interpolate as scinterp


import wdatmos
import spec_plot_tools as spt

wd=wdatmos.wdmodel(filename='ELM.hdf5')
model = wd(Teff = teff, logg = logg)
model_waves = model['w']
model_flux = model['flux']
model_spec = np.vstack([model_waves,model_flux])
teff = 7250
logg = 4.75

target_list_name = 'listFWCTB'
target_list = np.genfromtxt(target_list_name, dtype = 'str')
slit_width = 1.0 #arcseconds
pixel_scale = 0.3 #arcseconds per pixel_scale
slit_width = slit_width/pixel_scale #slit width in pixels

poly_degree = 5

first_conv_bin = 0.1 #width in angstroms of the first interpolation of the model to then be used in the convolution.
test_loc = 1200 #pixel location in the target spectrum to look to get a pixel to wavelength value to use for the seeing


low_wave_cut= 3800
high_wave_cut= 5200

#low_wave_cut = 4000
#high_wave_cut = 5200


#####
#####
#continuum_list = [[3597,3603],
                  #[3670,3678],
                  #[3782,3785],
                  #[3861,3864],
                  #[4014,4034],
                  #[4183, 4214],
                  #[4589,4608],
                  #[4645,4680],
                  #[4740,4760],
                  #[4930,4935],
                  #[5045,5070],
                  #[5110,5130],
                  #[5220,5240],
                  #[5275,5290]]
                  
                  
continuum_list = [[3809,3812],
                  [3861,3864],
                  [3907,3911],
                  [4014,4017],
                  [4036,4040],
                  [4183, 4214],
                  [4422,4427],
                  [4427,4432],
                  [4432,4437],
                  [4589,4608],
                  [4645, 4650],
                  [4655, 4660],
                  [4665, 4670],
                  [4675,4680],
                  [4720,4725],
                  [4730,4735],
                  [4740,4745],
                  [4750,4755],
                  [4760,4765],
                  [4770,4775],
                  [4970,4975],
                  [5045, 5050],
                  [5055,5060],
                  [5065,5070],
                  [5110,5130],
                  [5190,5195]]
                  
#continuum_list = [[4014,4034],
                  #[4183, 4214],
                  #[4589,4608],
                  #[4645,4680],
                  #[4740,4760],
                  #[4930,4935],
                  #[5045,5070],
                  #[5110,5130]]

target_continuum_list = [[3861,3864],
                  [3900,3915],
                  [4014,4034],
                  [4183, 4214],
                  [4589,4608],
                  [4645,4680],
                  [4740,4760],
                  [4930,4935],
                  [5045,5070],
                  [5110,5130],
                  [5187,5192]]




######



def get_doppler_shifted(wavelengths, radial_velocity):
    #print "doppler shifting by ", radial_velocity
    lambda_obs = wavelengths * (radial_velocity*u.km/u.s + const.c.to(u.km/u.s)) / const.c.to(u.km/u.s)
    return lambda_obs.value


def dopp_shift_continuum_list(radial_velocity):
    dopp_cont_list = []
    for waves in continuum_list:
        shift_waves = get_doppler_shifted(waves, radial_velocity)
        dopp_cont_list.append(shift_waves)
    return dopp_cont_list



def chi_squared(observed, actual):
    return (observed - actual)**2/actual

#def rebin_model(target_spec, model_spec)

def get_med_val(input_spec, wave_range):
    sub_spec = spt.trim_spec(input_spec, wave_range[0], wave_range[1])
    length = sub_spec[0].shape[0]
    if length%2 == 0:
        #even number
        sub_spec = sub_spec[:, :-1] #trim off the last point
    
    #med_val = np.nanmedian(sub_spec, axis =1)
    med_flux = np.nanmedian(sub_spec[1])
    med_index = np.where(sub_spec[1] == med_flux)[0]
    med_wave = sub_spec[0, med_index][0]
    med_val =[med_wave, med_flux]
    #print med_val
    return med_val

def make_continuum(input_spec, continuum_list= continuum_list):
    waves= []
    flux = []
    for ranges in continuum_list:
        new_vals = get_med_val(input_spec, ranges)
        waves.append(new_vals[0])
        flux.append(new_vals[1])
    wave_array = np.array(waves)
    flux_array = np.array(flux)
    continuum_spec = np.vstack([wave_array, flux_array])
    #plt.plot(input_spec[0], input_spec[1], label = 'input_spec')
    #plt.plot(continuum_spec[0], continuum_spec[1], linestyle = 'none', marker = 'o', label = 'continuum')
    #plt.legend()
    #plt.show()
    return continuum_spec

def get_norm_polynomial(input_spec, continuum_list = continuum_list):
    continuum_spec = make_continuum(input_spec, continuum_list = continuum_list)
    poly_coeffs= np.polyfit(continuum_spec[0], continuum_spec[1], poly_degree)
    #plt.plot(input_spec[0], input_spec[1], label = 'input_spec')
    #plt.plot(continuum_spec[0], continuum_spec[1], linestyle = 'none', marker = 'o', label = 'continuum', color = 'r')
    #plt.plot(input_spec[0], np.polyval(poly_coeffs, input_spec[0]), label = 'fit')
    #plt.title(continuum_list[0])
    #plt.legend()
    #plt.show()
    return poly_coeffs

def poly_norm_spec(input_spec, continuum_list = continuum_list):
    poly_coeffs = get_norm_polynomial(input_spec, continuum_list = continuum_list)
    poly_vals = np.polyval(poly_coeffs, input_spec[0])
    input_spec[1]= np.float_(input_spec[1])/poly_vals
    return input_spec

def calc_sq_dist(target_spec, model_spec, error_spec = np.array([])):
    interp_model_flux = np.interp(target_spec[0], model_spec[0], model_spec[1])
    interp_model= np.vstack([np.copy(target_spec[0]),interp_model_flux])
    #print "interp_model.shape", interp_model.shape
    if error_spec.shape[0] != 0:
        #norm_difs = np.abs(interp_model[1]-target_spec[1])/np.float_(error_spec[1])
        norm_difs = (interp_model[1]-target_spec[1])**2/np.float_(error_spec[1])**2
        #norm_difs = np.abs(interp_model[1]-target_spec[1])/np.float_(interp_model[1])
    else:
        #print "no uncertainties provided"
        norm_difs =(interp_model[1]-target_spec[1])**2/np.float_(interp_model[1])
    #norm_difs = np.abs(interp_model[1]-target_spec[1])

    nan_remove = np.isinf(norm_difs)
    norm_difs= norm_difs[~nan_remove]
    dif = np.sum(norm_difs)/norm_difs.shape[0]
    
    return dif


def convolve_model(model_spec, target_spec, header):
    """
    receive the fits file input of the target because you need a number of things from the header.
    """
    wavelengths = np.arange(np.nanmin(model_spec[0]), np.nanmax(model_spec[0]), first_conv_bin)
    #fluxes = scinterp.interp1d(wavelengths)
    fluxes = np.interp(wavelengths, model_spec[0], model_spec[1])
    dlam = target_spec[0][test_loc+1]-target_spec[0][test_loc] #angstroms per pixel at this location in the target
    see_sig = float(header['SEE_SIG']) #sigma value of gaussian fit to do the 
    see_sig = see_sig*dlam/first_conv_bin #seeing value in units of indices of the model
    pix_slit_width = slit_width*dlam/first_conv_bin  #slit width value in units of indices of the model
    #print "pix_slit_width", pix_slit_width, int(pix_slit_width)
    try:
        see_kernel = conv.Gaussian1DKernel(see_sig, x_size = int(pix_slit_width), mode = 'oversample')
        see_kernel.normalize()
        model_conv = conv.convolve(fluxes, see_kernel)
    except ValueError as error:
        #print error
        #print "so making it odd"
        pix_slit_width= pix_slit_width+1
        see_kernel = conv.Gaussian1DKernel(see_sig, x_size = int(pix_slit_width), mode = 'oversample')
        see_kernel.normalize()
        model_conv = conv.convolve(fluxes, see_kernel)
    pix_width = dlam/first_conv_bin #width in pixels of model of a pixel from the original spectrum
    #print "pix_width", pix_width
    pix_kernel = conv.Box1DKernel(width = int(pix_width), mode = 'oversample')
    pix_kernel.normalize()
    model_conv = conv.convolve(model_conv, pix_kernel)
    model_out = np.vstack([wavelengths, model_conv])
    return model_out


#===========================
##########################
#===========================


model_spec = poly_norm_spec(model_spec)



def fit_rv(target_file):
    return
