"""
this script should open a model file (or all of them I suppose more accurately), and step through them.

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
import kernel_builder

plt.rc('font', size =18)
#plt.rc('lines', markersize=12)
#plt.rc('font', size = 11)
plt.rc('lines', markersize = 5)

target_list_name = 'listFWCTB'
target_list = np.genfromtxt(target_list_name, dtype = 'str')
combined_spec_file = 'combined_PSRJ1431m4715_new.fits'
scaling_range = [4600,4650]
slit_width = 1.0 #arcseconds
pixel_scale = 0.3 #arcseconds per pixel_scale
slit_width = slit_width/pixel_scale #slit width in pixels

spectrum_type= 'combined'
#options are 'combined', 'single run'

ca_mask = [3920,3946]
weird2_mask= [4485,4507]
weird_mask=[4563,4576]
red_metal_mask= [5162,5192]

free_parameters= 2

output_names = "teff, logg, rv, chi_square"
#output_filename= 'chi_square_values.csv'
#output_filename= 'chi_square_values_noca.csv'
output_filename= 'chi_square_values_modeldiv.csv'


#teff = 7250
#logg = 6.0

teff = 7500
logg = 5.5
#teff = 7250
#logg = 5.50
#teff = 14750
#logg = 3.75
#logg = 6.25

#teff = 6000
#logg = 3.75
plot_fit = False

poly_degree = 5

first_conv_bin = 0.1 #width in angstroms of the first interpolation of the model to then be used in the convolution.
#first_conv_bin = 0.01  #width in angstroms of the first interpolation of the model to then be used in the convolution.

test_loc = 1200 #pixel location in the target spectrum to look to get a pixel to wavelength value to use for the seeing

#velocity_bound = 400 #km/s
velocity_step  = 50 #km
#velocity_step  = 5 #km
#velocity_tests = np.arange(-1*velocity_bound, velocity_bound+velocity_step, velocity_step)
velocity_low_bound = -500 #km/s
velocity_high_bound = 300 #km/s
velocity_tests = np.arange(velocity_low_bound, velocity_high_bound+velocity_step, velocity_step)

low_wave_cut= 3800
#high_wave_cut= 5200
high_wave_cut= 5050


#low_wave_cut= 3670
#high_wave_cut= 5200

#low_wave_cut = 4000
#high_wave_cut = 5200


#####
                  
                  
#continuum_list = [[3809,3812],
                  #[3861,3864],
                  #[3907,3911],
                  #[4014,4017],
                  #[4036,4040],
                  #[4183, 4214],
                  #[4422,4427],
                  #[4427,4432],
                  #[4432,4437],
                  #[4589,4608],
                  #[4645, 4650],
                  #[4655, 4660],
                  #[4665, 4670],
                  #[4675,4680],
                  #[4720,4725],
                  #[4730,4735],
                  #[4740,4745],
                  #[4750,4755],
                  #[4760,4765],
                  #[4770,4775],
                  #[4970,4975],
                  #[5045, 5050],
                  #[5055,5060],
                  #[5065,5070],
                  #[5110,5130],
                  #[5175,5180],
                  #[5190,5195]]#Best one there is. Before 2018-08-10

#continuum_list = [[3809,3812],
                  #[3861,3864],
                  #[3907,3911],
                  #[4014,4017],
                  #[4036,4040],
                  #[4183, 4214],
                  #[4422,4427],
                  #[4427,4432],
                  #[4432,4437],
                  #[4589,4608],
                  #[4645, 4650],
                  #[4655, 4660],
                  #[4665, 4670],
                  #[4675,4680],
                  #[4720,4725],
                  #[4730,4735],
                  #[4740,4745],
                  #[4750,4755],
                  #[4760,4765],
                  #[4770,4775],
                  #[4970,4975],
                  #[5045, 5050],
                  #[5055,5060],
                  #[5065,5070],
                  #[5110,5130],
                  #[5190,5195]]#Best one there is. Removed a range that fell in red_metal_mask


continuum_list = [[3809,3812],
                  [3861,3864],
                  [3907,3911],
                  [4014,4017],
                  [4036,4040],
                  [4183, 4214],
                  [4422,4427],
                  [4427,4432],
                  [4432,4437],
                  [4450,4455],
                  [4455,4460],
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
                  [4930,4935],
                  [4940,4945],
                  [4945,4950],
                  [4970,4975],
                  [5035,5040],
                  [5040,5045],
                  [5045, 5050]]#Best one there is. shortened red side

#continuum_list = [[3678,3682],
                  #[3689,3692],
                  #[3698,3702],
                  #[3715,3718],
                  #[3725,3727],
                  #[3738,3741],
                  #[3754,3756],
                  #[3782,3785],
                  #[3809,3812],
                  #[3861,3864],
                  #[3907,3911],
                  #[4014,4017],
                  #[4036,4040],
                  #[4183, 4214],
                  #[4422,4427],
                  #[4427,4432],
                  #[4432,4437],
                  #[4589,4608],
                  #[4645, 4650],
                  #[4655, 4660],
                  #[4665, 4670],
                  #[4675,4680],
                  #[4720,4725],
                  #[4730,4735],
                  #[4740,4745],
                  #[4750,4755],
                  #[4760,4765],
                  #[4770,4775],
                  #[4970,4975],
                  #[5045, 5050],
                  #[5055,5060],
                  #[5065,5070],
                  #[5110,5130],
                  #[5175,5180],
                  #[5190,5195]] #extended wavelength range
   

#target_continuum_list = [[3861,3864],
                  #[3900,3915],
                  #[4014,4034],
                  #[4183, 4214],
                  #[4589,4608],
                  #[4645,4680],
                  #[4740,4760],
                  #[4930,4935],
                  #[5045,5070],
                  #[5110,5130],
                  #[5187,5192]]#old one before 2018-07-01
                  
                  
if spectrum_type== 'single run':
    mask_metals= False
        
    if mask_metals==False:
        mask_list = []
    elif mask_metals == True:
        mask_list = [ca_mask]+[weird2_mask]+[weird_mask]+[red_metal_mask]
        print "\n******************"
        print "WARNING: INDIVIDUAL RUN SPECTRA DON'T ACTUALLY HAVE METAL LINES MASKED!!"
        print "******************\n"
    target_continuum_list = [[3814,3820],
                    [3863,3870],
                    [3909,3924],
                    [4014,4034],
                    [4183, 4214],
                    [4589,4608],
                    [4645,4680],
                    [4740,4760],
                    [4930,4935],
                    [5045,5070],
                    [5110,5130],
                    [5187,5192]]
        
    flux_stack = []
    noise_stack = []
    #for index in range(25,31):
    #for index in range(3,9):
    for index in range(0,6):
    #for index in range(3,6):
        filename = target_list[index]
        print filename
        i=fits.open(filename)
        header = fits.getheader(filename)
        file_waves= i[0].data
        file_flux = i[1].data
        file_noise = i[3].data
        noise = file_noise #don't want to scale it yet since there will be the normalization later
        #noise = file_flux*file_noise
        flux_stack.append([file_flux])
        noise_stack.append([noise])
        
    target_waves = file_waves
    target_flux= np.nanmedian(flux_stack, axis=0)[0]
    target_noise = np.nanmedian(noise_stack, axis=0)[0]
    #print target_waves.shape
    #print target_flux.shape
    target_spec = np.vstack([target_waves, target_flux])
    noise_spec = np.vstack([target_waves, target_noise])
    target_spec = spt.trim_spec(target_spec, low_wave_cut, high_wave_cut)
    noise_spec = spt.trim_spec(noise_spec, low_wave_cut, high_wave_cut)
    target_file = target_list[0]
    print target_file

elif spectrum_type == 'combined':
    mask_metals = True
    if mask_metals==False:
        mask_list = []
    elif mask_metals == True:
        mask_list = [ca_mask]+[weird2_mask]+[weird_mask]+[red_metal_mask]
    target_continuum_list= continuum_list
    target_spec, header, noise_spec = spt.retrieve_spec(combined_spec_file)
    target_spec = spt.trim_spec(target_spec, low_wave_cut, high_wave_cut)
    noise_spec = spt.trim_spec(noise_spec, low_wave_cut, high_wave_cut)
    target_spec= spt.clean_spectrum(target_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
    noise_spec= spt.clean_spectrum(noise_spec, np.min(noise_spec[0]), np.max(noise_spec[0]), mask_list)
    print "\n*************"
    print "Combined spectrum, so only using RV=0; no RV stepping actually done."
    print "**************\n"
    velocity_tests= np.arange(0,velocity_step/2.,velocity_step) #this makes the only velocity that is used be 0 km/s, but you don't have to eliminate the for-loop.... I think/hope.

######


#flux_stack = []
#noise_stack = []
##for index in range(25,31):
##for index in range(3,9):
#for index in range(0,6):
##for index in range(3,6):
    #filename = target_list[index]
    #print filename
    #i=fits.open(filename)
    #header = fits.getheader(filename)
    #file_waves= i[0].data
    #file_flux = i[1].data
    #file_noise = i[3].data
    #noise = file_noise #don't want to scale it yet since there will be the normalization later
    ##noise = file_flux*file_noise
    #flux_stack.append([file_flux])
    #noise_stack.append([noise])
    
#target_waves = file_waves
#target_flux= np.nanmedian(flux_stack, axis=0)[0]
#target_noise = np.nanmedian(noise_stack, axis=0)[0]
##print target_waves.shape
##print target_flux.shape
#target_spec = np.vstack([target_waves, target_flux])
#noise_spec = np.vstack([target_waves, target_noise])
#target_spec = spt.trim_spec(target_spec, low_wave_cut, high_wave_cut)
#noise_spec = spt.trim_spec(noise_spec, low_wave_cut, high_wave_cut)
#target_file = target_list[0]
#print target_file
#i= fits.open(target_file)
#header = fits.getheader(target_file)
#target_waves= i[0].data
##target_flux = i[1].data
#target_spec, header, noise_spec = spt.retrieve_spec('combined_PSRJ1431m4715.fits')
#target_spec, header, noise_spec = spt.retrieve_spec('combined_PSRJ1431m4715_new.fits')
#target_spec = spt.trim_spec(target_spec, low_wave_cut, high_wave_cut)
#noise_spec = spt.trim_spec(noise_spec, low_wave_cut, high_wave_cut)
#######

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

def calc_sq_dist(target_spec, model_spec, error_spec = np.array([]), free_parameters= free_parameters):
    """
    Return the reduced chi-square value if provided using the error spectrum (already rescaled to the spectrum values) if provided; otherwise it will return the reduced chi-square values using the model values for the denominator.
    
    """
    #interp_model_flux = np.interp(target_spec[0], model_spec[0], model_spec[1])
    interpolator3= scinterp.CubicSpline(model_spec[0], model_spec[1])
    interp_model_flux= interpolator3(target_spec[0])
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
    #dif = np.sum(norm_difs)/norm_difs.shape[0]
    dif = np.sum(norm_difs)/(norm_difs.shape[0]-1-free_parameters) #based on Numerical Recipes in C page 621. (Section 14.3)

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
    #fluxes = np.interp(wavelengths, model_spec[0], model_spec[1])
    interpolator= scinterp.CubicSpline(model_spec[0], model_spec[1])
    fluxes = interpolator(wavelengths)
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

######

#David's instructions for loading the model
wd=wdatmos.wdmodel(filename='ELM.hdf5')

####3

#teff_array = np.arange(6000, 15000, 250)
#logg_array = np.arange(3.75, 6.5, 0.25)
teff_array=wd.Teffs
logg_array = wd.loggs
model = wd(Teff = teff, logg = logg)
#print wd.Teffs
#print wd.loggs
model_num =0

#####

model_waves = model['w']
model_flux = model['flux'] #since we'll be arbitrarily-ish scaling this it won't work.

model_spec  = np.vstack([model_waves, model_flux])
#target_spec = poly_norm_spec(target_spec, continuum_list=target_continuum_list)
#### Here's the target normalization step=========================
target_spec = spt.poly_norm_spec(target_spec, continuum_list=target_continuum_list, poly_degree = poly_degree, plot_all = True)

noise_spec[1]= noise_spec[1]*target_spec[1] #scale the noise spectrum with the flattened target spectrum.


rv_dist_list=[]
for radial_velocity in velocity_tests:
    test_model = np.copy(model_spec)
    test_model[0]=get_doppler_shifted(test_model[0], radial_velocity)
    test_model = convolve_model(test_model, target_spec, header)
    dopp_cont_list= dopp_shift_continuum_list(radial_velocity)
    #test_model = poly_norm_spec(test_model, continuum_list = dopp_cont_list)
    #### Here's the model normalization step==============================
    #test_model = spt.poly_norm_spec(test_model, continuum_list = dopp_cont_list, poly_degree = poly_degree, plot_all = plot_fit)
    test_model = spt.poly_norm_spec(test_model, continuum_list = dopp_cont_list, poly_degree = poly_degree, plot_all = False)
    #test_model = spt.rescale_spectrum(test_model, target_spec, scaling_range)
    new_rv_dist= calc_sq_dist(target_spec, test_model, error_spec = noise_spec)
    rv_dist_list.append(new_rv_dist)
rv_dist_array = np.array(rv_dist_list)
min_index = np.argmin(rv_dist_array)
#min_model = model_file_list[min_index]
#min_teff = teff_array[min_index]
#min_logg = logg_array[min_index]
#min_dist = dist_array[min_index]
min_rv = velocity_tests[min_index]
#mask_list=[]
#model_spec= get_model_fromfile(min_model)
min_model = wd(Teff= teff, logg = logg)

#model_waves= np.array(min_model['w'])
#print 'Model_waves'
#print model_waves
#wave_difs = model_waves-np.roll(model_waves, 1)
#print wave_difs
#plt.plot(model_waves, wave_difs)
#plt.show()
model_spec = np.vstack([min_model['w'], min_model['flux']])
model_spec[0] = get_doppler_shifted(model_spec[0], min_rv)
####model_spec = spt.trim_spec(model_spec, np.min(target_spec[0]), np.max(target_spec[0]))
model_spec = convolve_model(model_spec, target_spec, header)
####model_spec= poly_norm_spec(model_spec)
dopp_cont_list= dopp_shift_continuum_list(min_rv)
#######model normalization ==========================
model_spec= spt.poly_norm_spec(model_spec, continuum_list = dopp_cont_list, poly_degree= poly_degree)
#######model_spec = spt.rescale_spectrum(model_spec, target_spec, scaling_range)
model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
plot_overlays(target_spec, model_spec, model_string = 'Teff ' + str(teff) + ' logg ' +str(logg)+ ' RV '+ str(min_rv)+'km/s')
#plt.plot(model_waves, scale_model_flux, label = 'model'+str(teff) + ' ' + str(logg))
#plt.plot(target_waves, target_flux, label = 'Target')
#plt.legend()
#plt.xlabel('Angstroms')
#plt.ylabel('Flux in cgs 10**-16')
#plt.show()

def run_model_grid(target_spec):
    #mask_list = []
    #target_spec = spt.clean_spectrum(target_spec, min_wave, max_wave, mask_list)
    dist_list = []
    rv_list = []
    #for teff,logg in zip(teff_array, logg_array):
        #model = wd(Teff = teff , logg = logg)
        #model_spec = np.vstack([model['w'], model['flux']])
        #model_spec = convolve_model(model_spec, target_spec, header)
        #new_dist = calc_sq_dist(target_spec, model_spec)
        #dist_list.append(new_dist)
    for teff,logg in zip(teff_array, logg_array):
        print "Teff:", teff, "logg:", logg
        model = wd(Teff = teff , logg = logg)
        model_spec = np.vstack([model['w'], model['flux']])
        #insert the doppler shifting
        rv_dist_list=[]
        for radial_velocity in velocity_tests:
            test_model = np.copy(model_spec)
            test_model[0]=get_doppler_shifted(test_model[0], radial_velocity)
            test_model = convolve_model(test_model, target_spec, header)
            dopp_cont_list= dopp_shift_continuum_list(radial_velocity)
            #test_model = poly_norm_spec(test_model)
            #model normalization =================================
            test_model = spt.poly_norm_spec(test_model, continuum_list=dopp_cont_list, poly_degree= poly_degree)
            #test_model= spt.rescale_spectrum(test_model, target_spec, scaling_range)
            #new_rv_dist= calc_sq_dist(target_spec, test_model)
            #new_rv_dist = calc_sq_dist(target_spec, test_model, error_spec = noise_spec)
            new_rv_dist = calc_sq_dist(target_spec, test_model) #for error-free chi-square; uses model division

            rv_dist_list.append(new_rv_dist)
        rv_dist_array = np.array(rv_dist_list)
        min_rv_index= np.argmin(rv_dist_array)
        new_dist = np.copy(rv_dist_array[min_rv_index])
        new_rv = np.copy(velocity_tests[min_rv_index])
        if plot_fit:
            model_spec[0] = get_doppler_shifted(model_spec[0], new_rv)
            model_spec = convolve_model(model_spec, target_spec, header)
            dopp_cont_list= dopp_shift_continuum_list(new_rv)
            model_spec= spt.poly_norm_spec(model_spec, continuum_list = dopp_cont_list, poly_degree= poly_degree)
            model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
            plt.title(r'$\chi^2=$'+str(new_dist))
            plot_overlays(target_spec, model_spec, model_string = 'Teff ' + str(teff) + ' logg ' +str(logg)+ ' RV '+ str(new_rv)+'km/s')
        rv_list.append(new_rv)
        dist_list.append(new_dist)
    dist_array = np.array(dist_list)
    rv_array = np.array(rv_list)
    for teff,logg, dist_mod, rv in zip(teff_array, logg_array, dist_list, rv_list):
        print "Teff:", teff, "logg:", logg, "chi-squared:", dist_mod, "Radial_velocity:", rv
    output_array= np.vstack([teff_array, logg_array, rv_array, dist_array]).T
    print "Saving "+ output_filename
    np.savetxt(output_filename, output_array, header = output_names, delimiter= ',')
    print "Saved " + output_filename
    sorted_indices = np.argsort(dist_list)
    sorted_teff = teff_array[sorted_indices]
    sorted_logg= logg_array[sorted_indices]
    sorted_dist = dist_array[sorted_indices]
    sorted_rv = rv_array[sorted_indices]
    print "======== Sorted  by chi-squared =========="
    for teff,logg, dist_mod, rv in zip(sorted_teff, sorted_logg, sorted_dist, sorted_rv):
        print "Teff:", teff, "logg:", logg, "chi-squared:", dist_mod, "Radial_velocity:", rv
    min_index = np.argmin(dist_list)
    #min_model = model_file_list[min_index]
    min_teff = teff_array[min_index]
    min_logg = logg_array[min_index]
    min_dist = dist_array[min_index]
    min_rv = rv_array[min_index]
    print "best fit model:", "Teff", min_teff, "logg", min_logg, "chi-squared", min_dist, "Radial Velocity:", min_rv, 'km/s'
    #model_spec= get_model_fromfile(min_model)
    min_model = wd(Teff= min_teff, logg = min_logg)
    
    #model_waves= np.array(min_model['w'])
    #print 'Model_waves'
    #print model_waves
    #wave_difs = model_waves-np.roll(model_waves, 1)
    #print wave_difs
    #plt.plot(model_waves, wave_difs)
    #plt.show()
    dopp_cont_list= dopp_shift_continuum_list(min_rv)
    model_spec = np.vstack([min_model['w'], min_model['flux']])
    model_spec[0] = get_doppler_shifted(model_spec[0], min_rv)
    #model_spec = spt.trim_spec(model_spec, np.min(target_spec[0]), np.max(target_spec[0]))
    model_spec = convolve_model(model_spec, target_spec, header)
    #model_spec= poly_norm_spec(model_spec)
    #normalization of the spectrum================
    model_spec= spt.poly_norm_spec(model_spec, continuum_list=dopp_cont_list, poly_degree = poly_degree)
    #model_spec = spt.rescale_spectrum(model_spec, target_spec, scaling_range)
    model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
    plot_overlays(target_spec, model_spec, model_string = 'Teff ' + str(min_teff) + ' logg ' +str(min_logg)+ ' RV '+ str(min_rv)+'km/s')
    #interp_model_flux = np.interp(target_spec[0], model_spec[0], model_spec[1])
    interpolator2= scinterp.CubicSpline(model_spec[0], model_spec[1])
    interp_model_flux = interpolator2(target_spec[0])
    interp_model= np.vstack([np.copy(target_spec[0]),interp_model_flux])
    plot_overlays(target_spec,interp_model, model_string = 'interp Teff ' + str(min_teff) + ' logg ' +str(min_logg))
    plot_overlays(target_spec, noise_spec, model_string = 'noise')
    plot_overlays(model_spec,interp_model, model_string = 'interp Teff ' + str(min_teff) + ' logg ' +str(min_logg))
    plot_overlays_convolve(target_spec, model_spec, model_string = 'Teff ' + str(min_teff) + ' logg ' +str(min_logg))
    chi_square_countours(teff_array,logg_array, dist_array)
    plt.scatter(teff_array, dist_array)
    plt.xlabel("Teff (K)")
    plt.ylabel(r"red $\chi^2$")
    plt.show()
    
    teff_scale = -1*(teff_array-np.min(teff_array))
    plt.scatter(logg_array+np.random.rand(logg_array.shape[0]), dist_array, c = teff_scale)
    plt.xlabel("log(g)")
    plt.ylabel(r"red $\chi^2$")
    plt.show()
    #output_array = np.vstack([target_spec[0], target_spec[1], target_err[1]]).T
    #np.savetxt( 'output_spectrum.csv',output_array, header = 'Wavelength, Flux (cgs units), Error', delimiter = ',')

    #print model['w'][0]

print get_doppler_shifted(4000, -200)

#target_spec = poly_norm_spec(target_spec)
run_model_grid(target_spec)
