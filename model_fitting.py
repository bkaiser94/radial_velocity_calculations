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

#plt.rc('font', size =18)
#plt.rc('lines', markersize=10)
#plt.rc('lines', linewidth = 2)

target_list_name = 'listFWCTB'
target_list = np.genfromtxt(target_list_name, dtype = 'str')
scaling_range = [4600,4650]
slit_width = 1.0 #arcseconds
pixel_scale = 0.3 #arcseconds per pixel_scale
slit_width = slit_width/pixel_scale #slit width in pixels

first_conv_bin = 0.1 #width in angstroms of the first interpolation of the model to then be used in the convolution.
test_loc = 1200 #pixel location in the target spectrum to look to get a pixel to wavelength value to use for the seeing

#velocity_bound = 400 #km/s
velocity_step  = 50 #km
#velocity_tests = np.arange(-1*velocity_bound, velocity_bound+velocity_step, velocity_step)
velocity_low_bound = -500 #km/s
velocity_high_bound = 300 #km/s
velocity_tests = np.arange(velocity_low_bound, velocity_high_bound+velocity_step, velocity_step)

#low_wave_cut= 3800
#high_wave_cut= 5200

low_wave_cut = 4000
high_wave_cut = 5200


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
                  
                  
#continuum_list = [[3812,3815],
                  #[3861,3864],
                  #[3928,3929],
                  #[4014,4017],
                  #[4036,4040],
                  #[4052,4055],
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
                  #[4930,4935],
                  #[5045,5070],
                  #[5110,5130]]
                  
continuum_list = [[4014,4034],
                  [4183, 4214],
                  [4589,4608],
                  [4645,4680],
                  [4740,4760],
                  [4930,4935],
                  [5045,5070],
                  [5110,5130]]






######


flux_stack = []
for index in range(25,31):
#for index in range(3,9):
    filename = target_list[index]
    print filename
    i=fits.open(filename)
    header = fits.getheader(filename)
    file_waves= i[0].data
    file_flux = i[1].data
    flux_stack.append([file_flux])
    
target_waves = file_waves
target_flux= np.nanmedian(flux_stack, axis=0)[0]
print target_waves.shape
print target_flux.shape
target_spec = np.vstack([target_waves, target_flux])
target_spec = spt.trim_spec(target_spec, low_wave_cut, high_wave_cut)
target_file = target_list[0]
print target_file
#i= fits.open(target_file)
#header = fits.getheader(target_file)
#target_waves= i[0].data
#target_flux = i[1].data

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
    poly_coeffs= np.polyfit(continuum_spec[0], continuum_spec[1], 3)
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
    plt.axhline(y=1, label = 'y=1', color = 'cyan')
    #plt.plot(spec2[0], spec2[1], label = model_string, linestyle ='none', marker = 'o')
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
    #plt.plot(spec2[0], spec2conv, label = model_string, linestyle ='none', marker = 'o')

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

######

#David's instructions for loading the model
wd=wdatmos.wdmodel(filename='ELM.hdf5')
#teff = 9000
#logg = 5.25
teff = 14750
#logg = 3.75
logg = 6.25

#teff = 6000
#logg = 3.75
####3

#teff_array = np.arange(6000, 15000, 250)
#logg_array = np.arange(3.75, 6.5, 0.25)
teff_array=wd.Teffs
logg_array = wd.loggs
model = wd(Teff = teff, logg = logg)
print wd.Teffs
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
#target_spec = np.vstack([target_waves, target_flux])
target_spec = poly_norm_spec(target_spec)

#scale_factor= get_scale_factor(target_spec, model_spec)
#scale_model_flux = model_flux* scale_factor
#print scale_model_flux.mean()
#print target_flux.mean()
rv_dist_list=[]
for radial_velocity in velocity_tests:
    test_model = np.copy(model_spec)
    test_model[0]=get_doppler_shifted(test_model[0], radial_velocity)
    test_model = convolve_model(test_model, target_spec, header)
    dopp_cont_list= dopp_shift_continuum_list(radial_velocity)
    test_model = poly_norm_spec(test_model, continuum_list = dopp_cont_list)
    #scaling_coefficient= get_scale_factor(target_spec, test_model)
    #test_model[1]=test_model[1]*scaling_coefficient
    new_rv_dist= calc_sq_dist(target_spec, test_model)
    rv_dist_list.append(new_rv_dist)
rv_dist_array = np.array(rv_dist_list)
min_index = np.argmin(rv_dist_array)
#min_model = model_file_list[min_index]
#min_teff = teff_array[min_index]
#min_logg = logg_array[min_index]
#min_dist = dist_array[min_index]
min_rv = velocity_tests[min_index]
mask_list=[]
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
#model_spec = spt.trim_spec(model_spec, np.min(target_spec[0]), np.max(target_spec[0]))
model_spec = convolve_model(model_spec, target_spec, header)
model_spec= poly_norm_spec(model_spec)
#scaling_coefficient= get_scale_factor(target_spec, model_spec)
#model_spec[1]= model_spec[1]*scaling_coefficient
model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
#calc_rdist(scaling_coefficient)
plot_overlays(target_spec, model_spec, model_string = 'Teff ' + str(teff) + ' logg ' +str(logg)+ ' RV '+ str(min_rv)+'km/s')
#plt.plot(model_waves, scale_model_flux, label = 'model'+str(teff) + ' ' + str(logg))
#plt.plot(target_waves, target_flux, label = 'Target')
#plt.legend()
#plt.xlabel('Angstroms')
#plt.ylabel('Flux in cgs 10**-16')
#plt.show()

def run_model_grid(target_spec):
    mask_list = []
    #target_spec = spt.clean_spectrum(target_spec, min_wave, max_wave, mask_list)
    dist_list = []
    rv_list = []
    #for teff,logg in zip(teff_array, logg_array):
        #model = wd(Teff = teff , logg = logg)
        #model_spec = np.vstack([model['w'], model['flux']])
        #model_spec = convolve_model(model_spec, target_spec, header)
        #scaling_coefficient= get_scale_factor(target_spec, model_spec)
        #model_spec[1]=model_spec[1]*scaling_coefficient
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
            test_model = poly_norm_spec(test_model)
            #scaling_coefficient= get_scale_factor(target_spec, test_model)
            #test_model[1]=test_model[1]*scaling_coefficient
            new_rv_dist= calc_sq_dist(target_spec, test_model)
            rv_dist_list.append(new_rv_dist)
        rv_dist_array = np.array(rv_dist_list)
        min_rv_index= np.argmin(rv_dist_array)
        new_dist = np.copy(rv_dist_array[min_rv_index])
        new_rv = np.copy(velocity_tests[min_rv_index])
        #new_dist = calc_sq_dist(target_spec, model_spec)
        rv_list.append(new_rv)
        dist_list.append(new_dist)
    dist_array = np.array(dist_list)
    rv_array = np.array(rv_list)
    for teff,logg, dist_mod, rv in zip(teff_array, logg_array, dist_list, rv_list):
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
    model_spec = np.vstack([min_model['w'], min_model['flux']])
    model_spec[0] = get_doppler_shifted(model_spec[0], min_rv)
    #model_spec = spt.trim_spec(model_spec, np.min(target_spec[0]), np.max(target_spec[0]))
    model_spec = convolve_model(model_spec, target_spec, header)
    model_spec= poly_norm_spec(model_spec)
    scaling_coefficient= get_scale_factor(target_spec, model_spec)
    model_spec[1]= model_spec[1]*scaling_coefficient
    model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
    #calc_rdist(scaling_coefficient)
    plot_overlays(target_spec, model_spec, model_string = 'Teff ' + str(min_teff) + ' logg ' +str(min_logg)+ ' RV '+ str(min_rv)+'km/s')
    interp_model_flux = np.interp(target_spec[0], model_spec[0], model_spec[1])
    interp_model= np.vstack([np.copy(target_spec[0]),interp_model_flux])
    plot_overlays(target_spec,interp_model, model_string = 'interp Teff ' + str(min_teff) + ' logg ' +str(min_logg))
    plot_overlays(model_spec,interp_model, model_string = 'interp Teff ' + str(min_teff) + ' logg ' +str(min_logg))
    plot_overlays_convolve(target_spec, model_spec, model_string = 'Teff ' + str(min_teff) + ' logg ' +str(min_logg))
    chi_square_countours(teff_array,logg_array, dist_array)
    #output_array = np.vstack([target_spec[0], target_spec[1], target_err[1]]).T
    #np.savetxt( 'output_spectrum.csv',output_array, header = 'Wavelength, Flux (cgs units), Error', delimiter = ',')

    #print model['w'][0]

print get_doppler_shifted(4000, -200)

#target_spec = poly_norm_spec(target_spec)
run_model_grid(target_spec)
