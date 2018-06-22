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


teff = 7250
logg = 5.25
mask_list = []
wd=wdatmos.wdmodel(filename='ELM.hdf5')
model = wd(Teff = teff, logg = logg)
model_waves = model['w']
model_flux = model['flux']
model_spec = np.vstack([model_waves,model_flux])

target_list_name = 'listFWCTB'
target_list = np.genfromtxt(target_list_name, dtype = 'str')
slit_width = 1.0 #arcseconds
pixel_scale = 0.3 #arcseconds per pixel_scale
slit_width = slit_width/pixel_scale #slit width in pixels

plot_fit = False

poly_degree = 5

first_conv_bin = 0.1 #width in angstroms of the first interpolation of the model to then be used in the convolution.
test_loc = 1200 #pixel location in the target spectrum to look to get a pixel to wavelength value to use for the seeing

#velocity_bound = 400 #km/s
velocity_step  = 100 #km
first_prev_velocity_step = 200
velocity_step_list = [200., 100., 10., 1., 0.1] #km/s (the first one doesn't actually get used except to set the outer bounds of the grid)
velocity_center = -100 #km/s
velocity_grid_radius = 8 #number of gridpoints away from the central one to include
overlap_radius = 4 #was 18
velocity_low_bound = -500 #km/s
velocity_high_bound = 300 #km/s
velocity_tests = np.arange(velocity_low_bound, velocity_high_bound+velocity_step, velocity_step)

low_wave_cut= 3800
high_wave_cut= 5200

#low_wave_cut = 4000
#high_wave_cut = 5200


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
def make_velocity_grid(velocity_center, velocity_step, prev_velocity_step, overlap_radius= overlap_radius):
#def make_velocity_grid(velocity_center, velocity_step, velocity_grid_radius= velocity_grid_radius):
    """
    Produce a np.arange() that covers the desired velocity range.
    """
    low_bound = velocity_center - prev_velocity_step*overlap_radius
    high_bound= velocity_center +prev_velocity_step *(overlap_radius +1)
    grid = np.arange(low_bound, high_bound, velocity_step)
    return grid

def get_doppler_shifted(wavelengths, radial_velocity):
    #print "doppler shifting by ", radial_velocity
    lambda_obs = wavelengths * (radial_velocity*u.km/u.s + const.c.to(u.km/u.s)) / const.c.to(u.km/u.s)
    return lambda_obs.value

#def wave2doppler(w, w0):
    #w0_equiv = u.doppler_optical(w0)
    #w_equiv = w.to(u.km/u.s, equivalencies=w0_equiv)
    #return w_equiv.to(u.km/u.s)

#print(wave2doppler(waveclosetoHa, 656.489 * u.nm).to(u.km/u.s))
def dopp_shift_continuum_list(radial_velocity):
    dopp_cont_list = []
    for waves in continuum_list:
        shift_waves = get_doppler_shifted(waves, radial_velocity)
        dopp_cont_list.append(shift_waves)
    return dopp_cont_list


def chi_squared(observed, actual):
    return (observed - actual)**2/actual

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

def retrieve_target_spec(filename):
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
    file_spec = spt.trim_spec(file_spec, low_wave_cut, high_wave_cut)
    file_noise_spec = spt.trim_spec(file_noise_spec, low_wave_cut, high_wave_cut)
    file_spec = spt.poly_norm_spec(file_spec, continuum_list = target_continuum_list, poly_degree = poly_degree)
    file_noise_spec[1] = file_spec[1]*file_noise_spec[1]
    return file_spec, header, file_noise_spec


def make_rest_spectrum(filename, radial_velocity):
    i = fits.open(filename)
    header= fits.getheader(filename)
    file_waves= i[0].data
    file_flux = i[1].data
    file_sky = i[2].data
    file_noise = i[3].data
    rest_waves = get_doppler_shifted(file_waves, -1*radial_velocity)
    return
#===========================
##########################
#===========================


#model_spec = spt.poly_norm_spec(model_spec, poly_degree= poly_degree)


def fit_rv(target_file):
    return

def minimize_velocity(model_spec, target_spec, noise_spec, target_header, velocity_center, velocity_tests, plot_fit = False):
    """
    Test the whole grid and output the optimal radial velocity for the given target spectrum at the specified grid resolution
    """
    rv_dist_list=[]
    #print velocity_tests
    for radial_velocity in velocity_tests:
        test_model = np.copy(model_spec)
        test_model[0]=get_doppler_shifted(test_model[0], radial_velocity)
        test_model = convolve_model(test_model, target_spec, target_header)
        dopp_cont_list= dopp_shift_continuum_list(radial_velocity)
        #test_model = poly_norm_spec(test_model)
        test_model = spt.poly_norm_spec(test_model, continuum_list=dopp_cont_list, poly_degree= poly_degree)
        new_rv_dist= calc_sq_dist(target_spec, test_model, error_spec = noise_spec)
        rv_dist_list.append(new_rv_dist)
    rv_dist_array = np.array(rv_dist_list)
    
    min_rv_index= np.argmin(rv_dist_array)
    
    new_dist = np.copy(rv_dist_array[min_rv_index])
    new_rv = np.copy(velocity_tests[min_rv_index])
    if plot_fit:
        
        plt.plot(velocity_tests, rv_dist_array)
        plt.plot(new_rv, new_dist,marker = 'o', linestyle ='none', color = 'r')
        plt.show()
    return new_rv
    

def iterate_resolutions(model_spec, target_file ):
    target_spec, target_header, noise_spec = retrieve_target_spec(target_file) 
    best_rv = velocity_center
    for index in range(1, len(velocity_step_list)):
        #if index >0:
            #prev_velocity_step = velocity_step_list[index-1]
        #else:
            #prev_velocity_step = first_prev_velocity_step
        prev_velocity_step = velocity_step_list[index-1]
        velocity_step = velocity_step_list[index]
        velocity_tests= make_velocity_grid(best_rv, velocity_step, prev_velocity_step)
        if index == len(velocity_step_list)-1:
            best_rv = minimize_velocity(model_spec, target_spec, noise_spec, target_header, best_rv, velocity_tests, plot_fit = plot_fit)
        else:
            best_rv = minimize_velocity(model_spec, target_spec, noise_spec, target_header, best_rv, velocity_tests, plot_fit = False)
        print "==========="
        print target_file, "best_rv: ", best_rv
    if plot_fit:
        test_model = np.copy(model_spec)
        test_model[0]=get_doppler_shifted(test_model[0], best_rv)
        test_model = convolve_model(test_model, target_spec, target_header)
        dopp_cont_list= dopp_shift_continuum_list(best_rv)
        #test_model = poly_norm_spec(test_model)
        test_model = spt.poly_norm_spec(test_model, continuum_list=dopp_cont_list, poly_degree= poly_degree)
        test_model= spt.clean_spectrum(test_model, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
        #norm_target_spec = spt.poly_norm_spec(target_spec, continuum_list = target_continuum_list, poly_degree = poly_degree)
        plt.title("RV: " + str(best_rv) + str(" km/s ") + target_file)
        plt.plot(target_spec[0], target_spec[1], color = 'b', label = "Target")
        plt.plot(test_model[0], test_model[1], color = 'r', label = "Model")
        plt.legend()
        plt.show()
    return best_rv


rv_list = []
time_list= []
for target_file in target_list:
    best_rv = iterate_resolutions(model_spec, target_file)
    print "###############"
    print target_file, " best_rv:", best_rv
    print "###############"
    #i=fits.open(target_file)
    header = fits.getheader(target_file)
    rv_list.append(best_rv)
    time_list.append(header['BMJD_TDB'])
    #target_spec, target_header = retrieve_target_spec(target_file) #first spectrum in the list for testing.
    ##rv_dist_list=[]
    #velocity_tests = make_velocity_grid(velocity_center, velocity_step_list[0], velocity_grid_radius = velocity_grid_radius)
    #for radial_velocity in velocity_tests:
        #test_model = np.copy(model_spec)
        #test_model[0]=get_doppler_shifted(test_model[0], radial_velocity)
        #test_model = convolve_model(test_model, target_spec, target_header)
        #dopp_cont_list= dopp_shift_continuum_list(radial_velocity)
        ##test_model = poly_norm_spec(test_model)
        #test_model = spt.poly_norm_spec(test_model, continuum_list=dopp_cont_list, poly_degree= poly_degree)
        #new_rv_dist= calc_sq_dist(target_spec, test_model)
        #rv_dist_list.append(new_rv_dist)
    #rv_dist_array = np.array(rv_dist_list)
    #min_rv_index= np.argmin(rv_dist_array)
    #new_dist = np.copy(rv_dist_array[min_rv_index])
    #new_rv = np.copy(velocity_tests[min_rv_index])
    #new_rv = minimize_velocity(model_spec, target_spec, target_header, velocity_center, velocity_tests)
    #print new_rv
    #new_dist = calc_sq_dist(target_spec, model_spec)
    #rv_list.append(new_rv)
    #dist_list.append(new_dist)
    #dist_array = np.array(dist_list)
    #rv_array = np.array(rv_list)
    
rv_array = np.array(rv_list)
time_array = np.array(time_list)
print rv_array
print time_array

plt.scatter(time_array, rv_array)
plt.ylabel('RV (km/s)')
plt.xlabel("BMJD_TDB")
plt.show()

out_array = np.vstack([time_array,rv_array])
np.savetxt('rv_plot.txt', out_array.T, delimiter =',', header = 'Times(BMJD_TDB), RV (km/s)')
