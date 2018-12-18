"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-12-18

This should contain the chi-square functions and model convolution functions that are needed by both model_fitting.py and model_rv_fit.py, so that we guarantee they are using the same methods to calculate the chi-square values and do fits.


"""
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
import scipy.stats as scistats
#import cosmics
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const
from astropy import convolution as conv
import scipy.interpolate as scinterp



#import wdatmos
import spec_plot_tools as spt
import kernel_builder



def interpolate_model(target_spec, model_spec):
    interp_model_flux = np.interp(target_spec[0], model_spec[0], model_spec[1])
    #interpolator3= scinterp.CubicSpline(model_spec[0], model_spec[1])
    #interp_model_flux= interpolator3(target_spec[0])
    interp_model= np.vstack([np.copy(target_spec[0]),interp_model_flux])
    #print "interp_model.shape", interp_model.shape
    return interp_model

def calc_residuals(target_spec, model_spec):
    """
    model_spec should NOT already be convolved an interpolated!!!
    returns a residual spectrum with the wavelengths of target_spec[0]
    """
    interp_model= interpolate_model(target_spec, model_spec)
    residuals= target_spec[1]-interp_model[1]
    return np.vstack([target_spec[0], residuals])

def calc_sq_dist(target_spec, model_spec, error_spec = np.array([]), free_parameters= free_parameters, norm=chi_norm, raw_chi= raw_chi):
    """
    Return the reduced chi-square value if provided using the error spectrum (already rescaled to the spectrum values) if provided; otherwise it will return the reduced chi-square values using the model values for the denominator.
    
    """
    #target_spec= np.copy(in_target_spec)
    #model_spec= np.copy(in_model_spec)
    #error_spec= np.copy(in_error_spec)
    #interp_model_flux = np.interp(target_spec[0], model_spec[0], model_spec[1])
    #interpolator3= scinterp.CubicSpline(model_spec[0], model_spec[1])
    #interp_model_flux= interpolator3(target_spec[0])
    #interp_model= np.vstack([np.copy(target_spec[0]),interp_model_flux])
    #print "interp_model.shape", interp_model.shape
    interp_model= interpolate_model(target_spec, model_spec)
    if norm:
        if error_spec.shape[0] != 0:
            #have to rescale the errors to the normalized values
            error_spec[1]= np.copy(error_spec[1]/np.sum(target_spec[1]))
        target_spec[1]= np.float_(target_spec[1])/np.sum(target_spec[1])
        interp_model[1]= np.float_(interp_model[1])/np.sum(interp_model[1])
    else:
        pass
    if error_spec.shape[0] != 0:
        #norm_difs = np.abs(interp_model[1]-target_spec[1])/np.float_(error_spec[1])
        norm_difs = (interp_model[1]-target_spec[1])**2/np.float_(error_spec[1])**2
        #norm_difs = np.abs(interp_model[1]-target_spec[1])/np.float_(interp_model[1])
    else:
        print "no uncertainties provided"
        norm_difs =(interp_model[1]-target_spec[1])**2/np.float_(interp_model[1])
    #norm_difs = np.abs(interp_model[1]-target_spec[1])

    #nan_remove = np.isinf(norm_difs)
    #norm_difs= norm_difs[~nan_remove]
    #dif = np.sum(norm_difs)/norm_difs.shape[0]
    #print "norm_difs.shape[0]:", norm_difs.shape[0]
    if raw_chi:
        dif = np.sum(norm_difs)/norm_difs.shape[0]
        #dif =np.sum(norm_difs)
    else:
        dif = np.sum(norm_difs)/(norm_difs.shape[0]-1-free_parameters) #based on Numerical Recipes in C page 621. (Section 14.3)

    return dif
