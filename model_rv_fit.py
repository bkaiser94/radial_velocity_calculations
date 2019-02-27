"""
using a model-template already fitted to a median-combined group of observations, we now find the radial velocity of each individual observation and zero it out to the rest frame while appending a radial velocity value to the header in addition to the model that was used for the 
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
from astropy.table import Table, Column
import cosmics
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const
from astropy import convolution as conv
import scipy.interpolate as scinterp
import time
start = time.time()
#print start
import wdatmos
import spec_plot_tools as spt
import model_manipulation as mm
#plt.rc('font', size =18)
#plt.rc('lines', markersize=12)
#plt.rc('font', size = 12)
#plt.rc('lines', markersize = 12)


#teff = 7250
#logg = 5.25


#teff = 7250
#logg = 6.0

astropy_in= False #indicator of whether or not the logg and teff are customized for each observation

teff = 7500
logg = 5.5

#output_filename= 'rv_plot_second_iteration.txt'
#output_filename= 'rv_plot_teff7500_logg550.txt'
#output_filename= 'rv_plot_teff7500_logg550_20181011.txt'
#output_filename= '20190206_rv_rand_teff7500_logg550.txt'
output_filename= '20190213_rv_teff7500_logg550_full.csv'
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

plot_fit = True

mc_jump = 1 #number of layers of velocity grid to skip for the Monte Carlo evaluation. Probably want to be >0
#num_mc = 100 #number of randomized spectra to produce for each target spectrum
num_mc = 2 #number of randomized spectra to produce for each target spectrum
poly_degree = 5

first_conv_bin = 0.1 #width in angstroms of the first interpolation of the model to then be used in the convolution.
test_loc = 1200 #pixel location in the target spectrum to look to get a pixel to wavelength value to use for the seeing

#velocity_bound = 400 #km/s
velocity_step  = 100 #km
first_prev_velocity_step = 200
#velocity_step_list = [200., 100., 10., 1., 0.1] #km/s (the first one doesn't actually get used except to set the outer bounds of the grid)
velocity_step_list = [200., 100., 10.]#km/s (the first one doesn't actually get used except to set the outer bounds of the grid)
velocity_center = -100 #km/s
velocity_grid_radius = 8 #number of gridpoints away from the central one to include
overlap_radius = 4 #was 18
velocity_low_bound = -500 #km/s
velocity_high_bound = 300 #km/s
velocity_tests = np.arange(velocity_low_bound, velocity_high_bound+velocity_step, velocity_step)

sample_points=1000

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
##########################
col_names= ['filename',
            'bmjd_tdb', 
            'rv',
            'rv_error',
            'teff_used',
            'logg_used']

dtype_list=['S75',
            'f',
            'f',
            'f',
            'i',
            'f']
unit_list= ['none',
            'days',
            'km/s',
            'km/s',
            'K',
            'log in cgs']



######

def find_time_difference(begin, end):
    difference = end-begin
    hours = int(difference) / 3600
    minutes = int(difference%3600)/60
    seconds = difference%60
    print "Runtime: ", hours, 'h', minutes, 'm', seconds, 's'

def make_velocity_grid(velocity_center, velocity_step, prev_velocity_step, overlap_radius= overlap_radius):
#def make_velocity_grid(velocity_center, velocity_step, velocity_grid_radius= velocity_grid_radius):
    """
    Produce a np.arange() that covers the desired velocity range.
    """
    low_bound = velocity_center - prev_velocity_step*overlap_radius
    high_bound= velocity_center +prev_velocity_step *(overlap_radius +1)
    grid = np.arange(low_bound, high_bound, velocity_step)
    return grid

def make_rand_velocity_grid(velocity_center,  prev_velocity_step, overlap_radius= overlap_radius, sample_points=sample_points):
    """
    generate a random distribution of velocities to test in a certain range. Actually, this should only take bounds
    or something as arguments
    """
    low_bound = velocity_center - prev_velocity_step*overlap_radius
    high_bound= velocity_center +prev_velocity_step *(overlap_radius)
    print "high_bound", high_bound, "low_bound", low_bound
    scaling_factor= high_bound-low_bound
    grid = np.random.rand(sample_points)*scaling_factor+low_bound
    print "min in rand grid", np.min(grid)
    print "max in rand grid", np.max(grid)
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

def convolve_model(model_spec, target_spec, header):
    """
    receive the fits file input of the target because you need a number of things from the header.
    """
    wavelengths = np.arange(np.nanmin(model_spec[0]), np.nanmax(model_spec[0]), first_conv_bin)
    #fluxes = scinterp.interp1d(wavelengths)
    fluxes = np.interp(wavelengths, model_spec[0], model_spec[1])
    #interpolator = scinterp.CubicSpline(model_spec[0], model_spec[1])
    #fluxes= interpolator(wavelengths)
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
def remove_bad_noise(target_spec, noise_spec, scaled_noise):
    """
    Basically, this just takes the values from the spectrum where the noise was scaled to be negative for 
    whatever reason and reassigns them as a positive value that is really large so that they are ignored in chi-
    squared fitting. 
    """
    return



def minimize_velocity(model_spec, target_spec, noise_spec, target_header, velocity_tests, plot_fit = False, last_test= False):
    """
    Test the whole grid and output the optimal radial velocity for the given target spectrum at the specified grid resolution
    """
    rv_dist_list=[]
    red_rv_dist_list=[]
    #print velocity_tests
    for radial_velocity in velocity_tests:
        test_model = np.copy(model_spec)
        test_model[0]=get_doppler_shifted(test_model[0], radial_velocity)
        test_model = mm.convolve_model(test_model, target_spec, target_header)
        dopp_cont_list= dopp_shift_continuum_list(radial_velocity)
        #test_model = poly_norm_spec(test_model)
        test_model = spt.poly_norm_spec(test_model, continuum_list=dopp_cont_list, poly_degree= poly_degree)
        #new_rv_dist= calc_sq_dist(target_spec, test_model, error_spec = noise_spec)
        #new_red_rv_dist= mm.calc_sq_dist(target_spec, test_model, error_spec = noise_spec, free_parameters= 1, raw_chi=False)
        new_rv_dist= mm.calc_sq_dist(target_spec, test_model, error_spec = noise_spec, free_parameters= 1, raw_chi=True)
        rv_dist_list.append(new_rv_dist)
    rv_dist_array = np.array(rv_dist_list)
    
    min_rv_index= np.argmin(rv_dist_array)
    if last_test:
        test_model = np.copy(model_spec)
        test_model[0]=get_doppler_shifted(test_model[0], velocity_tests[min_rv_index])
        test_model = mm.convolve_model(test_model, target_spec, target_header)
        dopp_cont_list= dopp_shift_continuum_list(velocity_tests[min_rv_index])
        #test_model = poly_norm_spec(test_model)
        test_model = spt.poly_norm_spec(test_model, continuum_list=dopp_cont_list, poly_degree= poly_degree)
        #new_rv_dist= calc_sq_dist(target_spec, test_model, error_spec = noise_spec)
        #new_red_rv_dist= mm.calc_sq_dist(target_spec, test_model, error_spec = noise_spec, free_parameters= 1, raw_chi=False)
        chi_factor= mm.calc_sq_dist(target_spec, test_model, error_spec = noise_spec, free_parameters= 1, raw_chi=False) #this is the minimum reduced chi-square value that we'll use to scale the raw chi^2 values rescale chi-squared values
        print "Dividing raw chi2 values by", chi_factor
        rv_dist_array= rv_dist_array/chi_factor
        #plot_fit= True
    new_dist = np.copy(rv_dist_array[min_rv_index])
    new_rv = np.copy(velocity_tests[min_rv_index])
    print "new_rv:", new_rv, "new_dist", new_dist
    if plot_fit:
        #in_range= np.where(np.abs(velocity_tests-new_rv)<velocity_step_list[-1])
        #fit_params= np.polyfit(velocity_tests[in_range], rv_dist_array[in_range], 2)
        #fit_params= np.polyfit(velocity_tests, rv_dist_array, 2)
        ##xvals = np.linspace(np.min(velocity_tests[in_range]), np.max(velocity_tests[in_range]), 1000)
        #xvals = np.linspace(np.min(velocity_tests), np.max(velocity_tests), 1000)
        #yvals= np.polyval(fit_params, xvals)
        ##plt.plot(velocity_tests, rv_dist_array)
        #plt.scatter(velocity_tests, rv_dist_array)
        #plt.plot(xvals, yvals, color='r')
        #plt.plot(new_rv, new_dist,marker = 'o', linestyle ='none', color = 'r')
        plt.xlabel('Radial Velocity (km/s)')
        plt.ylabel(r'$\chi^2$')
        #mm.fit_fixed_parabola(velocity_tests, rv_dist_array)
        #plt.show()
    if last_test:
        sigma = mm.fit_fixed_parabola(velocity_tests, rv_dist_array, plot_fit= plot_fit, dof=1)
    
    if last_test:
        return new_rv, sigma
    else:
        return new_rv
    #return new_rv
    #return 

#def iterate_MC(model_spec, target_file, original_rv):
    #target_spec, target_header, noise_spec = retrieve_target_spec(target_file) 
    #wave_vals = target_spec[0]
    #scaled_noise = np.copy(noise_spec[1]*target_spec[1])
    #scaled_noise = np.abs(scaled_noise)
    #mc_rvs = []
    #print np.sort(noise_spec[1])[:5]
    #print np.sort(scaled_noise)[:5]
    #for jindex in range(0, num_mc):
        #best_rv = original_rv
        #random_flux = np.random.normal(target_spec[1], scaled_noise)
        #random_spec = np.vstack([wave_vals, random_flux])
        #for index in range(1+mc_jump, len(velocity_step_list)):
            #prev_velocity_step = velocity_step_list[index-1]
            #velocity_step = velocity_step_list[index]
            #velocity_tests= make_velocity_grid(best_rv, velocity_step, prev_velocity_step)
            #if index == len(velocity_step_list)-1:
                #best_rv = minimize_velocity(model_spec, random_spec, noise_spec, target_header, best_rv, velocity_tests, plot_fit = plot_fit)
            #else:
                #best_rv = minimize_velocity(model_spec, random_spec, noise_spec, target_header, best_rv, velocity_tests, plot_fit = False)
            ##print "==========="
            ##print target_file, "best_rv: ", best_rv
        #if plot_fit:
            #test_model = np.copy(model_spec)
            #test_model[0]=get_doppler_shifted(test_model[0], best_rv)
            #test_model = mm.convolve_model(test_model, random_spec, target_header)
            #dopp_cont_list= dopp_shift_continuum_list(best_rv)
            ##test_model = poly_norm_spec(test_model)
            #test_model = spt.poly_norm_spec(test_model, continuum_list=dopp_cont_list, poly_degree= poly_degree)
            #test_model= spt.clean_spectrum(test_model, np.min(random_spec[0]), np.max(random_spec[0]), mask_list)
            ##norm_target_spec = spt.poly_norm_spec(target_spec, continuum_list = target_continuum_list, poly_degree = poly_degree)
            #plt.title("RV: " + str(best_rv) + str(" km/s ") + target_file)
            #plt.plot(random_spec[0], random_spec[1], color = 'b', label = "Target")
            #plt.plot(test_model[0], test_model[1], color = 'r', label = "Model")
            #plt.legend()
            #plt.show()
        #mc_rvs.append(best_rv)
    #return mc_rvs

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
            #best_rv = minimize_velocity(model_spec, target_spec, noise_spec, target_header, velocity_tests, plot_fit = plot_fit)
            best_rv = minimize_velocity(model_spec, target_spec, noise_spec, target_header, velocity_tests, plot_fit = False)
        else:
            best_rv = minimize_velocity(model_spec, target_spec, noise_spec, target_header, velocity_tests, plot_fit = False)
        print "==========="
        print target_file, "best_rv: ", best_rv
    #now we do it for the random sampling...
    velocity_tests= make_rand_velocity_grid(best_rv, velocity_step_list[-1]) #it's just the very last step
    best_rv= minimize_velocity(model_spec, target_spec, noise_spec, target_header, velocity_tests, last_test= False)
    velocity_tests= make_rand_velocity_grid(best_rv, velocity_step_list[-1]) #it's just the very last step
    best_rv, sigma= minimize_velocity(model_spec, target_spec, noise_spec, target_header, velocity_tests, last_test= True, plot_fit= plot_fit)
    if plot_fit:
        test_model = np.copy(model_spec)
        test_model[0]=get_doppler_shifted(test_model[0], best_rv)
        test_model = mm.convolve_model(test_model, target_spec, target_header)
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
    return best_rv, sigma


rv_list = []
time_list= []
sigma_list = []
teff_list=[]
logg_list=[]
if not astropy_in:
    for target_file in target_list:
        very_begin = time.time()
        begin = time.time()
        best_rv, sigma = iterate_resolutions(model_spec, target_file)
        end= time.time()
        find_time_difference(begin, end)
        #i=fits.open(target_file)
        #begin= time.time()
        #mc_rvs = iterate_MC(model_spec, target_file, best_rv)
        #sigma = np.std(mc_rvs)
        #end = time.time()
        #find_time_difference(begin,end)
        print "###############"
        print target_file, " best_rv:", best_rv, '+/-', sigma
        #print target_file, " best_rv:", best_rv
        very_end = time.time()
        find_time_difference(very_begin, very_end)
        print "###############"
        header = fits.getheader(target_file)
        rv_list.append(best_rv)
        time_list.append(header['BMJD_TDB'])
        sigma_list.append(sigma)
        teff_list.append(teff)
        logg_list.append(logg)
    
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
sigma_array = np.array(sigma_list)
teff_array = np.array(teff_list)
logg_array = np.array(logg_list)
target_name_array= np.array(target_list)
array_list= [target_name_array,
             time_array,
             rv_array,
             sigma_array,
             teff_array,
             logg_array]
for thing in array_list:
    print thing


out_table = Table(array_list, names= col_names, dtype=dtype_list)
#out_table.meta['comments']=[','.join(unit_list)]
print rv_array
print time_array
stop = time.time()
out_table.pprint()
#out_array = np.vstack([time_array, rv_array])
#out_array = np.vstack([time_array,rv_array, sigma_array])

print "Saving the data... hopefully"
#np.savetxt('rv_plot.txt', out_array.T, delimiter =',', header = 'Times(BMJD_TDB), RV (km/s), Sigma (km/s)')
#np.savetxt(output_filename, out_array.T, delimiter =',', header = 'Times(BMJD_TDB), RV (km/s), Sigma (km/s)')
out_table.write(output_filename, format='ascii.csv')
#np.savetxt(output_filename, out_array.T, delimiter =',', header = 'Times(BMJD_TDB), RV (km/s)') 
print "Saved the data"


print "Start:", start
print "Stop:", stop
#difference = stop-start

#print "Difference: ", stop-start
find_time_difference(start,stop)
#except:
    #print 'astropy table failed.'
print "\n#################"
print "Scaled noise has absolute value used, so that needs to be fixed"
print "##################\n"
#plt.scatter(time_array, rv_array)
plt.errorbar(time_array, rv_array, yerr = sigma_array, color = 'b', marker= 'o', linestyle='none')
plt.ylabel('RV (km/s)')
plt.xlabel("BMJD_TDB")
plt.show()

#print "Difference: ", stop-start
#hours = int(difference) / 3600
#minutes = int(difference%3600)/60
#seconds = difference%60

#print "Runtime: ", hours, 'h', minutes, 'm', seconds, 's'

