import numpy as np

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

def remove_range( input_spec, bound_list):
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
    clean_spec= trim_spec(input_spec, min_wave, max_wave)
    for mask in mask_list:
        clean_spec= remove_range(clean_spec, mask)
    clean_spec= sort_spectrum(clean_spec)
    return clean_spec
