import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

percentile= 80

def trim_spec(input_spec, min_wave, max_wave):
    lower_indices = np.where(input_spec[0]< max_wave)
    trimmed_waves= input_spec[0][lower_indices]
    trimmed_flux= input_spec[1][lower_indices]
    upper_indices= np.where(trimmed_waves > min_wave)
    trimmed_waves= trimmed_waves[upper_indices]
    trimmed_flux = trimmed_flux[upper_indices]
    trimmed_spec= np.vstack([trimmed_waves, trimmed_flux])
    #print trimmed_spec.shape
    return trimmed_spec

def remove_range(input_spec, bound_list):
    """
    Removes wavelengths and flux values from the array that fall in the range specified by bound_list
    """
    wave_array= input_spec[0]
    other_array= input_spec[1]
    lower_bound = bound_list[0]
    upper_bound= bound_list[1]
    low_mask = np.where(wave_array < lower_bound)
    high_mask= np.where(wave_array > upper_bound)
    low_waves= wave_array[low_mask]
    high_waves= wave_array[high_mask]
    low_other = other_array[low_mask]
    high_other= other_array[high_mask]
    merge_waves= np.append(low_waves, high_waves)
    merge_other= np.append(low_other, high_other)
    return np.vstack([merge_waves, merge_other])

def sort_spectrum(input_spec):
    """
    Fixes the ordering of the wavelengths (and associated fluxes) so that the go from least to greatest, thereby removing errant lines through the plotted spectra and allowing for correct interpolation in model-fitting.
    """
    sort_indices = np.argsort(input_spec[0])
    sorted_waves= input_spec[0][sort_indices]
    sorted_flux= input_spec[1][sort_indices]
    sorted_spectrum= np.vstack([sorted_waves, sorted_flux])
    
    return sorted_spectrum



def clean_spectrum(input_spec, min_wave, max_wave, mask_list):
    """
    input_spec should be a vstack of wavelengths and the flux (or error)
    """
    clean_spec= trim_spec(input_spec, min_wave, max_wave)
    for mask in mask_list:
        clean_spec= remove_range(clean_spec, mask)
    clean_spec= sort_spectrum(clean_spec)
    return clean_spec


def get_med_val(input_spec, wave_range):
    sub_spec = trim_spec(input_spec, wave_range[0], wave_range[1])
    length = sub_spec[0].shape[0]
    if length%2 == 0:
        #even number
        sub_spec = sub_spec[:, :-1] #trim off the last point
    
    #med_val = np.nanmedian(sub_spec, axis =1)
    #med_flux = np.nanmedian(sub_spec[1])
    med_flux = np.percentile(sub_spec[1], percentile)
    med_index = np.where(sub_spec[1] == med_flux)[0]
    try:
        med_wave = sub_spec[0, med_index][0]
    except IndexError:
        min_index = np.argmin(np.abs(med_flux-sub_spec[1]))
        med_wave = sub_spec[0][min_index]
    med_val =[med_wave, med_flux]
    #print med_val
    return med_val

def make_continuum(input_spec, continuum_list= []):
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

def get_norm_polynomial(input_spec, continuum_list = [], poly_degree = 3, plot_all = False):
    continuum_spec = make_continuum(input_spec, continuum_list = continuum_list)
    poly_coeffs= np.polyfit(continuum_spec[0], continuum_spec[1], poly_degree)
    if plot_all:
        plt.plot(input_spec[0], input_spec[1], label = 'input_spec')
        plt.plot(continuum_spec[0], continuum_spec[1], linestyle = 'none', marker = 'o', label = 'continuum', color = 'r')
        plt.plot(input_spec[0], np.polyval(poly_coeffs, input_spec[0]), label = 'fit')
        plt.title(continuum_list[0])
        plt.legend()
        plt.show()
    else:
        pass
    return poly_coeffs

def poly_norm_spec(input_spec, continuum_list = [], poly_degree = 3, plot_all  = False):
    poly_coeffs = get_norm_polynomial(input_spec, continuum_list = continuum_list, poly_degree = poly_degree, plot_all = plot_all)
    poly_vals = np.polyval(poly_coeffs, input_spec[0])
    input_spec[1]= np.float_(input_spec[1])/poly_vals
    return input_spec



def retrieve_spec(filename):
    """
    Input: filename for the target spectrum you want to get
    
    Output: Spectrum made of a 2xN numpy array, header of the fits file you loaded it from
    """
    #print filename
    i=fits.open(filename)
    header = fits.getheader(filename)
    file_waves= i[0].data
    file_flux = i[1].data
    file_noise = i[3].data
    file_spec = np.vstack([file_waves, file_flux])
    file_noise_spec = np.vstack([file_waves, file_noise])
    file_noise_spec[1] = file_spec[1]*file_noise_spec[1]
    return file_spec, header, file_noise_spec

def rescale_spectrum(input_spec, reference_spec, scale_range):
    input_value = get_med_val(input_spec, scale_range)[1]
    reference_value = get_med_val(reference_spec, scale_range)[1]
    scale_factor = reference_value/np.float_(input_value)
    input_spec[1] = input_spec[1]*scale_factor
    return input_spec
    
    
    
    
    
    
    
    
    
    
    
    
