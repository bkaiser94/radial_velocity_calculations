"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-12-18

This should contain the chi-square functions and model convolution functions that are needed by both model_fitting.py and model_rv_fit.py, so that we guarantee they are using the same methods to calculate the chi-square values and do fits. Hopefully they are the correct methods of doing those calculations, but they will certainly be the same methods.


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
#import kernel_builder

fixed_minimum=False #Determines if delta-chi2 parabola fitting is done with the minimum of the parabola fixed to a gridpoint or if it's allowed to be a free parameter (technically two free parameters since it's the x-position and then the chi2 value of that xpostion)

first_conv_bin = 0.1 #width in angstroms of the first interpolation of the model to then be used in the convolution.
test_loc = 1200 #pixel location in the target spectrum to look to get a pixel to wavelength value to use for the seeing
#slit_width = 1.0 #arcseconds
pixel_scale = 0.3 #arcseconds per pixel_scale
#slit_width = slit_width/pixel_scale #slit width in pixels
slit_width=3.0


delta_chi2= np.array(
    [np.nan,
     1.,
     2.30,
     3.53,
     4.72,
     5.82])
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

def calc_sq_dist(target_spec, model_spec, error_spec = np.array([]), free_parameters= 2, norm=False, raw_chi= False):
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
        print("no uncertainties provided")
        norm_difs =(interp_model[1]-target_spec[1])**2/np.float_(interp_model[1])
    #norm_difs = np.abs(interp_model[1]-target_spec[1])

    #nan_remove = np.isinf(norm_difs)
    #norm_difs= norm_difs[~nan_remove]
    #dif = np.sum(norm_difs)/norm_difs.shape[0]
    #print "norm_difs.shape[0]:", norm_difs.shape[0]
    if raw_chi:
        #dif = np.sum(norm_difs)/norm_difs.shape[0]
        dif =np.sum(norm_difs)
    else:
        dif = np.sum(norm_difs)/(norm_difs.shape[0]-1-free_parameters) #based on Numerical Recipes in C page 621. (Section 14.3)

    return dif


def plot_overlays(spec1, spec2, model_string = 'model'):
    plt.plot(spec1[0], spec1[1], label = 'observed')
    #plt.errorbar(spec1[0],spec1[1], yerr = errors[1], label='observed')
    plt.plot(spec2[0], spec2[1], label= model_string, color = 'r')
    plt.axhline(y=1, label = 'y=1', color = 'cyan')
    #plt.plot(spec2[0], spec2[1], label = model_string, linestyle ='none', marker = 'o')
    plt.legend(numpoints=1, fontsize=14, loc='best' )
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux (Arbitrary Units)')
    #plt.title(target_file )
    plt.show()
    return ''

def plot_overlays_convolve(spec1, spec2, model_string = 'model'):
    #plt.plot(spec1[0], spec1[1], label = 'observed')
    plt.plot(spec1[0], conv.convolve(spec1[1], conv.Gaussian1DKernel(3)), label = 'observed convolved')
    #plt.errorbar(spec1[0],spec1[1], yerr = errors[1], label='observed')
    #spec2conv = conv.convolve(spec2[1], conv.convolve(conv.Gaussian1DKernel(2.2), conv.Gaussian1DKernel(3)))
    #spec2conv =conv.convolve( conv.convolve(spec2[1], conv.Gaussian1DKernel(2.2)), conv.Gaussian1DKernel(5))
    spec2conv =conv.convolve(spec2[1], conv.Gaussian1DKernel(3))
    #plt.plot(spec2[0], spec2[1], label= model_string, color = 'r')
    plt.plot(spec2[0], spec2conv, label = model_string)
    #plt.plot(spec2[0], spec2conv, label = model_string, linestyle ='none', marker = 'o')

    plt.legend(numpoints=1, fontsize=14, loc='best' )
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux (Arbitrary Units)')
    #plt.title(target_file )
    plt.show()
    return ''

def convolve_model(model_spec, target_spec, header):
    """
    receive the fits file input of the target because you need a number of things from the header.
    """
    wavelengths = np.arange(np.nanmin(model_spec[0]), np.nanmax(model_spec[0]), first_conv_bin)
    #fluxes = scinterp.interp1d(wavelengths)
    fluxes = np.interp(wavelengths, model_spec[0], model_spec[1])
    #interpolator= scinterp.CubicSpline(model_spec[0], model_spec[1])
    #fluxes = interpolator(wavelengths)
    dlam = target_spec[0][test_loc+1]-target_spec[0][test_loc] #angstroms per pixel at this location in the target
    see_sig = float(header['SEE_SIG']) #sigma value of gaussian fit to do the 
    see_sig = see_sig*dlam/first_conv_bin #seeing value in units of indices of the model
    slit_width = slit_width/pixel_scale
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

def convolve_model_new(model_spec, header, slit_width=slit_width):
    """
    convolve a model spectrum in the new way in which we don't interpolate because it's already done the
    rebinned average flux-conservative method. Also, the model spectrum should already be in the scale of the
    target by this point, so there's no need for dlambda values... I think.
    """
    fluxes=model_spec[1]
    wavelengths=model_spec[0]
    pix_slit_width = slit_width/float(header['pix_scal']) #slit width in pixels
    see_sig = float(header['SEE_SIG']) #sigma value of gaussian fit to do the 
    try:
        see_kernel = conv.Gaussian1DKernel(see_sig, x_size = int(pix_slit_width), mode = 'oversample')
        see_kernel.normalize()
        model_conv = conv.convolve(fluxes, see_kernel,boundary='extend')
    except (ValueError,conv.utils.KernelSizeError) as error:
        #print error
        #print "so making it odd"
        pix_slit_width= pix_slit_width+1
        see_kernel = conv.Gaussian1DKernel(see_sig, x_size = int(pix_slit_width), mode = 'oversample')
        see_kernel.normalize()
        model_conv = conv.convolve(fluxes, see_kernel,boundary='extend')
    except TypeError as error:
        #print error
        #print "so making it odd"
        pix_slit_width= pix_slit_width+1
        see_kernel = conv.Gaussian1DKernel(see_sig, x_size = int(pix_slit_width), mode = 'oversample')
        see_kernel.normalize()
        model_conv = conv.convolve(fluxes, see_kernel,boundary='extend')
    model_out = np.vstack([wavelengths, model_conv])
    return model_out


def degrade_model(model_spec, target_spec, header):
    """
    receive the fits file input of the target because you need a number of things from the header.
    """
    #wavelengths = np.arange(np.nanmin(model_spec[0]), np.nanmax(model_spec[0]), first_conv_bin)
    wavelengths=model_spec[0]
    model_dlam=model_spec[0][1]-model_spec[0][0]
    #fluxes = scinterp.interp1d(wavelengths)
    #fluxes = np.interp(wavelengths, model_spec[0], model_spec[1])
    fluxes=model_spec[1]
    #interpolator= scinterp.CubicSpline(model_spec[0], model_spec[1])
    #fluxes = interpolator(wavelengths)
    dlam = target_spec[0][test_loc+1]-target_spec[0][test_loc] #angstroms per pixel at this location in the target
    see_sig = float(header['SEE_SIG']) #sigma value of gaussian fit to do the 
    #see_sig = see_sig*dlam/first_conv_bin #seeing value in units of indices of the model
    see_sig=see_sig*dlam/model_dlam
    #slit_width = slit_width/pixel_scale
    #pix_slit_width = slit_width*dlam/first_conv_bin  #slit width value in units of indices of the model
    pix_slit_width=slit_width/pixel_scale*dlam/model_dlam
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
    #pix_width = dlam/first_conv_bin #width in pixels of model of a pixel from the original spectrum
    #pix_width = dlam/model_dlam #width in pixels of model of a pixel from the original spectrum
    ##print "pix_width", pix_width
    #pix_kernel = conv.Box1DKernel(width = int(pix_width), mode = 'oversample')
    #pix_kernel.normalize()
    #model_conv = conv.convolve(model_conv, pix_kernel)
    new_model=np.vstack([model_spec[0],model_conv])
    model_out= spt.rebin_generic_spec(new_model, model_dlam, target_spec[0], dlam)
    #model_out = np.vstack([wavelengths, model_conv])
    return model_out


#def fit_fixed_parabola(xvals, yvals, dof=1, plot_fit= False, fixed_minimum=fixed_minimum):
def fit_parabola(xvals, yvals, dof=1, plot_fit= False, fixed_minimum=fixed_minimum):
    """
    Fit a parabola to the data after first affixing its minimum to be at the same location as the minimum of whatever data you're looking at.
    
    Should return the uncertainty values...
    """
    #first find the actual minimum of the data you're looking at...
    min_index= np.argmin(yvals)
    minx= np.copy(xvals[min_index])
    miny= np.copy(yvals[min_index])
    if fixed_minimum:
        def parabola_func(xvals, a): #version with fixed minimum at measured values
            """
            A parabola to be fitted with another function that calls curve_fit to fit this function to data and its chi2 vals
            """
            return a*(xvals-minx)**2+miny #version with fixed minimum at measured values
        popt, pcov= sciop.curve_fit(parabola_func, xvals, yvals) #version for fixed minimum at measured values
    else:
        def parabola_func(xvals, a, x_min, y_min):
            """
            A parabola to be fitted with another function that calls curve_fit to fit this function to data and its chi2 vals.
            Guesses for the fit are the minimum grid point, but it's allowed to vary
            """
            return a*(xvals-x_min)**2 + y_min
        popt, pcov= sciop.curve_fit(parabola_func, xvals, yvals, p0=[1,minx, miny])
        minx= popt[1]
        miny=popt[2]
    #popt, pcov= sciop.curve_fit(parabola_func, xvals, yvals) #version for fixed minimum at measured values
    sigma= np.sqrt(delta_chi2[dof]/popt[0])
    #Need to change the minimum location to be the minimum of the freely fitted parabola
    #just redefine minx and miny to be the minimum x and y values from the fitted parabola
    if plot_fit:
        plt.scatter(xvals, yvals)
        x_line= np.linspace(np.min(xvals), np.max(xvals), 1000)
        if fixed_minimum:
            y_line= parabola_func(x_line, popt)
        else:
            y_line= parabola_func(x_line, popt[0], popt[1], popt[2])
        plt.plot(x_line, y_line, color= 'r')
        plt.plot(minx,miny, marker='o', color='r')
        plt.show()
        #plt.scatter(xvals, yvals)
        x_line= np.linspace(np.min(xvals), np.max(xvals), 1000)
        #y_line= parabola_func(x_line, popt)-miny-1
        if fixed_minimum:
            y_line= parabola_func(x_line, popt)-miny-delta_chi2[dof]
        else:
            y_line= parabola_func(x_line, popt[0], popt[1], popt[2])-miny-delta_chi2[dof]
        plt.plot(x_line, y_line, color= 'r')
        plt.axhline(y=0, color = 'k', linestyle='--')
        bound_points= np.array([minx-sigma, minx, minx+sigma])
        #calc_bounds= parabola_func(bound_points, popt)-miny-1
        if fixed_minimum:
            calc_bounds= parabola_func(bound_points, popt)-miny-delta_chi2[dof]
        else:
            calc_bounds= parabola_func(bound_points, popt[0],popt[1],popt[2])-miny-delta_chi2[dof]
        print( "calc_bounds", calc_bounds)
        plt.plot(bound_points, calc_bounds, color='b')
        plt.scatter(bound_points, [0, -1*delta_chi2[dof], 0], color='b')
        #plt.plot(minx,miny, marker='o', color='r')
        plt.show()
    else:
        pass
    return minx, sigma
