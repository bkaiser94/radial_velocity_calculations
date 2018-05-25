"""
this script should open a model file (or all of them I suppose more accurately, and step through them.

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


import wdatmos
import spec_plot_tools as spt

target_list_name = 'listFWCTB'
target_list = np.genfromtxt(target_list_name, dtype = 'str')
scaling_range = [4600,4650]
slit_width = 1.0 #arcseconds
pixel_scale = 0.3 #arcseconds per pixel_scale
slit_width = slit_width/pixel_scale #slit width in pixels


flux_stack = []
for index in range(23,27):
    filename = target_list[index]
    print filename
    i=fits.open(filename)
    file_waves= i[0].data
    file_flux = i[1].data
    flux_stack.append([file_flux])
    
target_waves = file_waves
target_flux= np.nanmedian(flux_stack, axis=0)[0]
print target_waves.shape
print target_flux.shape

target_file = target_list[0]
print target_file
#i= fits.open(target_file)
#header = fits.getheader(target_file)
#target_waves= i[0].data
#target_flux = i[1].data

#######


def chi_squared(observed, actual):
    return (observed - actual)**2/actual

#def rebin_model(target_spec, model_spec)

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


def chi_square_countours(teff_array, logg_array, dist_array):
    min_index = np.argmin(dist_array)
    print "Teff and logg min chi-squared values: ", teff_array[min_index],logg_array[min_index], "|chi-sq:", dist_array[min_index]
    contour_array = np.vstack([teff_array,logg_array, dist_array])
    #plt.imshow(contour_array, aspect= 100)
    #plt.contour(teff_array, logg_array, dist_array)
    marker_scale = 1/dist_array* dist_array.min() *40.
    #plt.scatter(teff_array, logg_array, s= 1./dist_array*30, c = 1./dist_array*20)
    plt.scatter(teff_array, logg_array, s=marker_scale, c = marker_scale)
    plt.plot(teff_array[min_index],logg_array[min_index], marker = '*', markersize = 14)
    plt.xlabel('T_eff')
    plt.ylabel('logg')
    plt.show()

def get_scale_factor(target_spec, model_spec, scaling_range=scaling_range):
    model_scale_region = spt.trim_spec(model_spec,scaling_range[0],scaling_range[1])
    target_scale_region = spt.trim_spec(target_spec, scaling_range[0], scaling_range[1])
    scale_factor= np.mean(target_scale_region[1, :])/np.mean(model_scale_region[1, :])
    return scale_factor

def plot_overlays(spec1, spec2, model_string = 'model'):
    plt.plot(spec1[0], spec1[1], label = 'observed')
    #plt.errorbar(spec1[0],spec1[1], yerr = errors[1], label='observed')
    plt.plot(spec2[0], spec2[1], label= model_string, color = 'r')
    plt.legend(numpoints=1, fontsize=14, loc='best' )
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux (cgs units)')
    plt.title(target_file )
    plt.show()
    return ''

def plot_overlays_convolve(spec1, spec2, model_string = 'model'):
    #plt.plot(spec1[0], spec1[1], label = 'observed')
    plt.plot(spec1[0], conv.convolve(spec1[1], conv.Gaussian1DKernel(3)), label = 'observed convolved')
    #plt.errorbar(spec1[0],spec1[1], yerr = errors[1], label='observed')
    #spec2conv = conv.convolve(spec2[1], conv.convolve(conv.Gaussian1DKernel(2.2), conv.Gaussian1DKernel(3)))
    spec2conv =conv.convolve( conv.convolve(spec2[1], conv.Gaussian1DKernel(2.2)), conv.Gaussian1DKernel(5))
    #plt.plot(spec2[0], spec2[1], label= model_string, color = 'r')
    plt.plot(spec2[0], spec2conv, label = model_string)
    plt.legend(numpoints=1, fontsize=14, loc='best' )
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux (cgs units)')
    plt.title(target_file )
    plt.show()
    return ''

def convolve_model(model_spec, target_spec, header):
    """
    receive the fits file input of the target because you need a number of things from the header.
    """
    see_sig = float(header['SEE_SIG']) #sigma value of gaussian fit to do the 
    see_kernel = conv.Gaussian1DKernel(see_sig)
    #slit_kernel = conv.
    return

######

#David's instructions for loading the model
wd=wdatmos.wdmodel(filename='ELM.hdf5')
teff = 8000
logg = 6.25
####3

#teff_array = np.arange(6000, 15000, 250)
#logg_array = np.arange(3.75, 6.5, 0.25)
teff_array=wd.Teffs
logg_array = wd.loggs
model = wd(Teff = teff, logg = logg)
print wd.Teffs[0]
print wd.loggs
model_num =0

#print model
#for teff,logg in zip(teff_array, logg_array):
    #model = wd(Teff = teff , logg = logg)
    ##print model['w'][0]
    #if model != None:
        #model_num +=1

#####

model_waves = model['w']
model_flux = model['flux'] #since we'll be arbitrarily-ish scaling this it won't work.

model_spec  = np.vstack([model_waves, model_flux])
target_spec = np.vstack([target_waves, target_flux])

scale_factor= get_scale_factor(target_spec, model_spec)
scale_model_flux = model_flux* scale_factor
print scale_model_flux.mean()
print target_flux.mean()


plt.plot(model_waves, scale_model_flux, label = 'model'+str(teff) + ' ' + str(logg))
plt.plot(target_waves, target_flux, label = 'Target')
plt.legend()
plt.xlabel('Angstroms')
plt.ylabel('Flux in cgs 10**-16')
plt.show()

def run_model_grid(target_spec):
    mask_list = []
    #target_spec = spt.clean_spectrum(target_spec, min_wave, max_wave, mask_list)
    dist_list = []
    for teff,logg in zip(teff_array, logg_array):
        model = wd(Teff = teff , logg = logg)
        model_spec = np.vstack([model['w'], model['flux']])
        scaling_coefficient= get_scale_factor(target_spec, model_spec)
        model_spec[1]=model_spec[1]*scaling_coefficient
        new_dist = calc_sq_dist(target_spec, model_spec)
        dist_list.append(new_dist)
    dist_array = np.array(dist_list)
    for teff,logg, dist_mod in zip(teff_array, logg_array, dist_list):
        print "Teff:", teff, "logg:", logg, "chi-squared:", dist_mod
    min_index = np.argmin(dist_list)
    #min_model = model_file_list[min_index]
    min_teff = teff_array[min_index]
    min_logg = logg_array[min_index]
    min_dist = dist_array[min_index]
    print "best fit model:", "Teff", min_teff, "logg", min_logg, "chi-squared", min_dist
    #model_spec= get_model_fromfile(min_model)
    min_model = wd(Teff= min_teff, logg = min_logg)
    model_spec = np.vstack([min_model['w'], min_model['flux']])
    #model_spec = spt.trim_spec(model_spec, np.min(target_spec[0]), np.max(target_spec[0]))
    model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
    scaling_coefficient= get_scale_factor(target_spec, model_spec)
    model_spec[1]= model_spec[1]*scaling_coefficient
    #calc_rdist(scaling_coefficient)
    plot_overlays(target_spec, model_spec, model_string = 'Teff ' + str(min_teff) + ' logg ' +str(min_logg))
    interp_model_flux = np.interp(target_spec[0], model_spec[0], model_spec[1])
    interp_model= np.vstack([np.copy(target_spec[0]),interp_model_flux])
    plot_overlays(target_spec,interp_model, model_string = 'interp Teff ' + str(min_teff) + ' logg ' +str(min_logg))
    plot_overlays(model_spec,interp_model, model_string = 'interp Teff ' + str(min_teff) + ' logg ' +str(min_logg))
    plot_overlays_convolve(target_spec, model_spec, model_string = 'Teff ' + str(min_teff) + ' logg ' +str(min_logg))
    chi_square_countours(teff_array,logg_array, dist_array)
    #output_array = np.vstack([target_spec[0], target_spec[1], target_err[1]]).T
    #np.savetxt( 'output_spectrum.csv',output_array, header = 'Wavelength, Flux (cgs units), Error', delimiter = ',')

    #print model['w'][0]
    
run_model_grid(target_spec)
