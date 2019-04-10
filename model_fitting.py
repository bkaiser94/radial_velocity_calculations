"""
this script should open a model file (or all of them I suppose more accurately), and step through them.

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
from astropy.table import Table, Column
import scipy.interpolate as scinterp



import wdatmos
import spec_plot_tools as spt
#import kernel_builder
import model_manipulation as mm

plt.rc('font', size =18)
#plt.rc('lines', markersize=12)
#plt.rc('font', size = 11)
plt.rc('lines', markersize = 5)

target_list_name = 'listFWCTB'
target_list = np.genfromtxt(target_list_name, dtype = 'str')

combined_spec_file='fwctb.0220_J1431m4715_930_blue.fits'
#combined_spec_file='fwctb.0251_J1431m4715_930_blue_1arcsec.fits'

#astropy_input= '20190210_rv_teff7500_logg550_full.csv'
#astropy_input='20190301_rv_teff7500_logg550.csv'
#astropy_input='20190320_rv_teff7000_logg500.csv'
astropy_input ='20190408_rv_teff7000_logg500.csv'
#astropy_input='20190319_rv_teff7000_logg500.csv'
#scaling_range = [4600,4650]
rms_range= [4600,4650]
slit_width = 1.0 #arcseconds
pixel_scale = 0.3 #arcseconds per pixel_scale
slit_width = slit_width/pixel_scale #slit width in pixels

spectrum_type= 'individual'
#spectrum_type= 'single run'
#options are 'combined', 'single run', 'individual'

chi_norm= False
raw_chi= True
fixed_minimum =False

#ca_mask = [3920,3946]
ca_mask = [3920,4006]
#ca_mask=[1,1]
weird2_mask= [4485,4507]
weird_mask=[4563,4576]
red_metal_mask= [5162,5192]
#red_metal_mask= [1,1]

free_parameters= 2

#output_names = "teff, logg, rv, chi_square"
output_names = "teff, logg, rv, chi_square, revised_chi_square"

output_filename= '20190408_new_x2.csv'
output_astropy= '20190408_model_fits_nomask.csv'


cerro_pachon_location = coords.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)

#teff = 7250
#logg = 6.0

teff = 6000
#logg = 6.25
#teff = 7250
#logg = 5.50
#teff = 14750
logg = 3.75
#logg = 6.25

#teff = 6000
#logg = 3.75
plot_fit = False
output_grid= True
poly_degree = 5
#poly_degree = 7

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

#low_wave_cut= 3800
low_wave_cut=3600
#low_wave_cut= 2000
#high_wave_cut= 5500
high_wave_cut= 5200
#high_wave_cut= 5050


#low_wave_cut= 3670
#high_wave_cut= 5200

#low_wave_cut = 4000
#high_wave_cut = 5200


#####

#David's instructions for loading the model
wd_og=wdatmos.wdmodel(filename='ELM.hdf5')

def wd(Teff= teff, logg= logg):
    return np.copy(wd_og(Teff=Teff, logg=logg))

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
   

#################
#These are different from continuum_list because these are the large segments of the spectrum that should be ignored during the model minimization
#continuum_masks = [[3500,3815],
                   #[4008,4068],
                   #[4141,4293],
                   #[4363,4821],
                   #[4917,6000]] #from 8/23/18

continuum_masks=[[3500,3815],
                 [3845,3874],
                 [3897,3962],
                 [3980,4090],
                 [4110,4329],
                 [4346,4850],
                 [4884,6000]]



###################
col_names= ['filename',
            'bmjd_tdb', 
            'rv',
            'rv_error',
            'teff',
            'teff_error',
            'logg',
            'logg_error',
            'baryv_corr']

dtype_list=['S75',
            'f8',
            'f8',
            'f8',
            'f8',
            'f8',
            'f8',
            'f8',
            'f8']
unit_list= ['none',
            'days',
            'km/s',
            'km/s',
            'K',
            'K',
            'log in cgs',
            'log in cgs',
            'km/s']


##################
                  
                  
if spectrum_type== 'single run':
    balmer_only= False
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
    #for index in range(8,10):
    #for index in range(4,9):
    #for index in range(2,4):
    #for index in range(3,9):
    for index in range(0,6):
    #for index in range(3,6):
    #for index in range(10,21):
        filename = target_list[index]
        print filename
        i=fits.open(filename)
        header = fits.getheader(filename)
        file_waves= i[0].data
        file_flux = i[1].data
        file_noise = i[3].data
        noise = file_noise #don't want to scale it yet since there will be the normalization later
        flux_stack.append([file_flux])
        noise_stack.append([noise])
    noise_stack=np.array(noise_stack)
    target_waves = file_waves
    target_flux= np.nanmedian(flux_stack, axis=0)[0]
    target_noise = np.sum(noise_stack**2,axis=0)[0]
    print target_noise.shape
    target_noise= np.copy(np.sqrt(target_noise/np.float_(noise_stack.shape[0])))
    target_spec = np.vstack([target_waves, target_flux])
    noise_spec = np.vstack([target_waves, target_noise])
    noise_spec[1]=noise_spec[1]/target_spec[1]
    target_spec = spt.trim_spec(target_spec, low_wave_cut, high_wave_cut)
    noise_spec = spt.trim_spec(noise_spec, low_wave_cut, high_wave_cut)
    target_file = target_list[0]
    print target_file

elif spectrum_type == 'combined':
    #balmer_only = True 
    balmer_only = False

    mask_metals = True
    
    if mask_metals==False:
        mask_list = []
    elif mask_metals == True:
        mask_list = [ca_mask]+[weird2_mask]+[weird_mask]+[red_metal_mask]
    target_continuum_list= continuum_list
    target_spec, header, noise_spec = spt.retrieve_spec(combined_spec_file)
    target_spec = spt.trim_spec(target_spec, low_wave_cut, high_wave_cut)
    
    noise_spec = spt.trim_spec(noise_spec, low_wave_cut, high_wave_cut)
    #now undo the de-normalization of the noise spectrum that is done by retrieve_spec
    noise_spec[1]= noise_spec[1]/target_spec[1]
    target_spec= spt.clean_spectrum(target_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
    noise_spec= spt.clean_spectrum(noise_spec, np.min(noise_spec[0]), np.max(noise_spec[0]), mask_list)
    print "\n*************"
    print "Combined spectrum, so only using RV=0; no RV stepping actually done."
    print "**************\n"
    velocity_tests= np.arange(0,velocity_step/2.,velocity_step) #this makes the only velocity that is used be 0 km/s, but you don't have to eliminate the for-loop.... I think/hope.
    
elif spectrum_type == 'individual':
    #balmer_only = True 
    balmer_only = False
    output_grid= False
    mask_metals = False
    if mask_metals==False:
        mask_list = []
    elif mask_metals == True:
        mask_list = [ca_mask]+[weird2_mask]+[weird_mask]+[red_metal_mask]
    target_continuum_list= continuum_list
    target_spec, header, noise_spec = spt.retrieve_spec(combined_spec_file)
    target_spec = spt.trim_spec(target_spec, low_wave_cut, high_wave_cut)
    
    noise_spec = spt.trim_spec(noise_spec, low_wave_cut, high_wave_cut)
    #now undo the de-normalization of the noise spectrum that is done by retrieve_spec
    noise_spec[1]= noise_spec[1]/target_spec[1]
    target_spec= spt.clean_spectrum(target_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
    noise_spec= spt.clean_spectrum(noise_spec, np.min(noise_spec[0]), np.max(noise_spec[0]), mask_list)
    print "\n*************"
    print "individual spectra, so  using single RV provided in astropy file."
    print "**************\n"
    velocity_tests= np.arange(0,velocity_step/2.,velocity_step) #this makes the only velocity that is used be 0 km/s, but you don't have to eliminate the for-loop.... I think/hope.

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

def barycentric_vel(header):
    ra = header['RA']
    dec = header['DEC']
    radec = coords.SkyCoord(ra, dec, frame = 'icrs', unit= (u.hourangle, u.deg))
    bary_corr = radec.radial_velocity_correction(obstime= Time(header['DATE-OBS'], format = 'isot', scale= 'utc'), location = cerro_pachon_location)
    bary_corr = bary_corr.to(u.km/u.s)
    return bary_corr.value

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

def create_test_model(target_spec, header, logg, teff, radial_velocity):
    """
    
    """
    full_model= np.copy(wd(Teff=teff, logg=logg))
    test_waves= np.copy(full_model['w'])
    test_flux= np.copy(full_model['flux'])
    test_model=np.vstack([test_waves, test_flux])
    test_model[0]=get_doppler_shifted(test_model[0], radial_velocity)
    test_model = mm.convolve_model(test_model, target_spec, header)
    
    return

def extract_match_param(dist_array, min_index, param_array, other_param_array):
    """
    Make the row/column of the grid to be fit for the delta chi square parabola
    """
    try:
        param_array= param_array.value #make sure not a quantity
    except AttributeError:
        print "param_array didn't have value"
    try:
        other_param_array = other_param_array.value
    except:
        print "other_param_array didn't have value"
    min_param= param_array[min_index]
    matched_inds= np.where(param_array==min_param)
    other_param_row= np.copy(other_param_array[matched_inds])
    dist_row= np.copy(dist_array[matched_inds])
    if plot_fit:
        plt.plot(other_param_row, dist_row)
        plt.title(min_param)
        plt.show()
    else:
        pass
    return other_param_row, dist_row

def process_target_spec():
    
    return

######

##David's instructions for loading the model
#wd=wdatmos.wdmodel(filename='ELM.hdf5')

####3

teff_array=wd_og.Teffs
logg_array = wd_og.loggs
model = wd(Teff = teff, logg = logg)

model_num =0

#####

model_waves = model['w']
model_flux = model['flux'] #since we'll be arbitrarily-ish scaling this it won't work.

model_spec  = np.vstack([model_waves, model_flux])
#### Here's the target normalization step=========================
rescale_model= np.vstack([model_spec[0], model_spec[1]*(np.nanmax(target_spec[1])/np.nanmax(model_spec[1]))])
plot_overlays_convolve(target_spec, rescale_model)
target_spec = spt.poly_norm_spec(target_spec, continuum_list=target_continuum_list, poly_degree = poly_degree, plot_all = True)
noise_spec[1]= noise_spec[1]*target_spec[1] #scale the noise spectrum with the flattened target spectrum.

if balmer_only:
    print "\n*************\nOnly using the Balmer lines.\n****************\n"
    line_spec= spt.clean_spectrum(target_spec, np.min(target_spec[0]), np.max(target_spec[0]), continuum_masks) #the spectrum only including the balmer lines
    line_noise_spec= spt.clean_spectrum(noise_spec, np.min(noise_spec[0]), np.max(noise_spec[0]), continuum_masks)
else:
    pass

rv_dist_list=[]
for radial_velocity in velocity_tests:
    test_model = np.copy(model_spec)
    test_model[0]=get_doppler_shifted(test_model[0], radial_velocity)
    test_model = mm.convolve_model(test_model, target_spec, header)
    dopp_cont_list= dopp_shift_continuum_list(radial_velocity)
    #### Here's the model normalization step==============================
    test_model = spt.poly_norm_spec(test_model, continuum_list = dopp_cont_list, poly_degree = poly_degree, plot_all = False)
    if balmer_only:
        new_rv_dist= mm.calc_sq_dist(line_spec, test_model, error_spec = line_noise_spec, free_parameters= 2, norm=chi_norm, raw_chi= raw_chi)
    else:
        new_rv_dist= mm.calc_sq_dist(target_spec, test_model, error_spec = noise_spec, free_parameters=free_parameters, norm= chi_norm, raw_chi = raw_chi)
    rv_dist_list.append(new_rv_dist)
rv_dist_array = np.array(rv_dist_list)
min_index = np.argmin(rv_dist_array)

min_rv = velocity_tests[min_index]

min_model = wd(Teff= teff, logg = logg)

model_spec = np.vstack([min_model['w'], min_model['flux']])
model_spec[0] = get_doppler_shifted(model_spec[0], min_rv)
model_spec = mm.convolve_model(model_spec, target_spec, header)
dopp_cont_list= dopp_shift_continuum_list(min_rv)
#######model normalization ==========================
model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
model_spec= spt.poly_norm_spec(model_spec, continuum_list = dopp_cont_list, poly_degree= poly_degree, plot_all=plot_fit)
#model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
if balmer_only:
    model_spec = spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), continuum_masks)
    plot_overlays(line_spec, model_spec, model_string = 'Teff ' + str(teff) + ' logg ' +str(logg)+ ' RV '+ str(min_rv)+'km/s')
else:
    plot_overlays(target_spec, model_spec, model_string = 'Teff ' + str(teff) + ' logg ' +str(logg)+ ' RV '+ str(min_rv)+'km/s')

def run_model_grid(target_spec, velocity_tests= velocity_tests, noise_spec= np.array([]), plot_all = True, mask_list= mask_list):
    
    dist_list = []
    rv_list = []
    
    for teff,logg in zip(teff_array, logg_array):
        #print "Teff:", teff, "logg:", logg
        model = wd(Teff = teff , logg = logg)
        model_spec = np.vstack([model['w'], model['flux']])
        #insert the doppler shifting
        rv_dist_list=[]
        for radial_velocity in velocity_tests:
            test_model = np.copy(model_spec)
            test_model[0]=get_doppler_shifted(test_model[0], radial_velocity)
            test_model = mm.convolve_model(test_model, target_spec, header)
            dopp_cont_list= dopp_shift_continuum_list(radial_velocity)
            #model normalization =================================
            test_model= spt.clean_spectrum(test_model, np.min(target_spec[0]), np.max(target_spec[0]),mask_list)
            #test_model = spt.poly_norm_spec(test_model, continuum_list=dopp_cont_list, poly_degree= poly_degree)
            test_model = spt.poly_norm_spec(test_model, continuum_list=dopp_cont_list, poly_degree= poly_degree, radial_velocity= radial_velocity)
            if balmer_only:
                
                new_rv_dist = mm.calc_sq_dist(line_spec, test_model, error_spec = line_noise_spec, free_parameters=free_parameters, norm=chi_norm, raw_chi = raw_chi)

            else:
                #print "not balmer_only"
                #new_rv_dist = calc_sq_dist(target_spec, test_model, error_spec = noise_spec)
                new_rv_dist = mm.calc_sq_dist(target_spec, test_model, error_spec = noise_spec, free_parameters=free_parameters, norm= chi_norm, raw_chi= raw_chi)
            #new_rv_dist = calc_sq_dist(target_spec, test_model) #for error-free chi-square; uses model division

            rv_dist_list.append(new_rv_dist)
        rv_dist_array = np.array(rv_dist_list)
        min_rv_index= np.argmin(rv_dist_array)
        new_dist = np.copy(rv_dist_array[min_rv_index])
        new_rv = np.copy(velocity_tests[min_rv_index])
        if plot_fit:
            model_spec[0] = get_doppler_shifted(model_spec[0], new_rv)
            model_spec = mm.convolve_model(model_spec, target_spec, header)
            dopp_cont_list= dopp_shift_continuum_list(new_rv)
            model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
            model_spec= spt.poly_norm_spec(model_spec, continuum_list = dopp_cont_list, poly_degree= poly_degree, radial_velocity=new_rv)
            plt.title(r'$\chi^2=$'+str(new_dist))
            
            if balmer_only:
                model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), continuum_masks)
                plot_overlays(line_spec, model_spec, model_string = 'Teff ' + str(teff) + ' logg ' +str(logg)+ ' RV '+ str(new_rv)+'km/s')
            model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
            plot_overlays(target_spec, model_spec, model_string = 'Teff ' + str(teff) + ' logg ' +str(logg)+ ' RV '+ str(new_rv)+'km/s')
            #plot_overlays_convolve(target_spec, model_spec,  model_string = 'Teff ' + str(teff) + ' logg ' +str(logg)+ ' RV '+ str(new_rv)+'km/s')
        rv_list.append(new_rv)
        dist_list.append(new_dist)
    dist_array = np.array(dist_list)
    rv_array = np.array(rv_list)
    ##########################
    ###Rescale the chi-square values and make them actually be the version that's not divided by N
    #in both cases
    min_index = np.argmin(dist_list)
    #min_model = model_file_list[min_index]
    min_teff = np.copy(teff_array[min_index])
    min_logg = np.copy(logg_array[min_index])
    min_dist = np.copy(np.copy(dist_array[min_index]))
    min_rv = rv_array[min_index]
    model = wd(Teff= min_teff, logg= min_logg)
    model_waves= model['w']
    model_flux= model['flux']
    radial_velocity= min_rv
    model_waves = get_doppler_shifted(model_waves, radial_velocity)
    model_spec= np.vstack([model_waves, model_flux])
    model_spec = mm.convolve_model(model_spec, target_spec, header)
    dopp_cont_list= dopp_shift_continuum_list(radial_velocity)
    model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]),mask_list)
    model_spec = spt.poly_norm_spec(model_spec, continuum_list=dopp_cont_list, poly_degree= poly_degree, radial_velocity=radial_velocity)
    min_red_chi2= mm.calc_sq_dist(target_spec, model_spec, error_spec= noise_spec, free_parameters= free_parameters, raw_chi=False)
    #print 'rescaling shenanigans'
    #print "target_spec.shape[1]:",target_spec.shape[1]
    #print "min_dist:",min_dist
    if not raw_chi:
        dist_array= dist_array*(target_spec.shape[1]-1-free_parameters) #undoing the division by N
    else:
        pass
        #dist_array= dist_array*target_spec.shape[1] #undoing the division by N
    print "np.nanmin(dist_array):",np.nanmin(dist_array)
    #rescale_dist= np.copy(dist_array/min_dist) #dividing all of the values by the minimum pseudo-reduced chi-square value
    rescale_dist= np.copy(dist_array/min_red_chi2)
    min_rescale_dist= rescale_dist[min_index]
    #print "best fit model:", "Teff", min_teff, "logg", min_logg, "chi-squared", min_dist, "Radial Velocity:", min_rv, 'km/s, rescaled chi-square:', min_rescale_dist
    logg_teff, logg_dist = extract_match_param(rescale_dist, min_index, logg_array, teff_array) #row of Teffs around logg min
    logg_teff_lims= np.where(logg_teff < logg_teff[np.argmax(logg_dist)])
    logg_teff= logg_teff[logg_teff_lims]
    logg_dist = logg_dist[logg_teff_lims]
    teff_logg, teff_dist= extract_match_param(rescale_dist, min_index, teff_array, logg_array) #row of loggs around Teff min
    #teff_sigma= mm.fit_fixed_parabola(logg_teff, logg_dist, dof = free_parameters, plot_fit= plot_all)
    #logg_sigma= mm.fit_fixed_parabola(teff_logg, teff_dist, dof= free_parameters, plot_fit= plot_all)
    min_teff, teff_sigma= mm.fit_parabola(logg_teff, logg_dist, dof = free_parameters, plot_fit= plot_all, fixed_minimum=fixed_minimum)
    min_logg, logg_sigma= mm.fit_parabola(teff_logg, teff_dist, dof= free_parameters, plot_fit= plot_all, fixed_minimum=fixed_minimum)
    
    
    
    
    ###################
    #for teff,logg, dist_mod, rv, resc_chisq in zip(teff_array, logg_array, dist_list, rv_list, rescale_dist):
    #for teff,logg, dist_mod, rv, resc_chisq in zip(teff_array, logg_array, dist_array, rv_list, rescale_dist):
        #print "Teff:", teff, "logg:", logg, "chi-squared:", dist_mod, "Radial_velocity:", rv, "rescaled chi-square:", resc_chisq
    #output_array= np.vstack([teff_array, logg_array, rv_array, dist_array]).T
    if output_grid:
        output_array= np.vstack([teff_array, logg_array, rv_array, dist_array, rescale_dist]).T
        print "Saving "+ output_filename
        np.savetxt(output_filename, output_array, header = output_names, delimiter= ',')
        print "Saved " + output_filename
    else:
        pass
    #sorted_indices = np.argsort(dist_list)
    #sorted_teff = teff_array[sorted_indices]
    #sorted_logg= logg_array[sorted_indices]
    #sorted_dist = dist_array[sorted_indices]
    #sorted_rv = rv_array[sorted_indices]
    #sorted_resc_chisq = rescale_dist[sorted_indices]
    #print "======== Sorted  by chi-squared =========="
    #for teff,logg, dist_mod, rv, resc_chisq in zip(sorted_teff, sorted_logg, sorted_dist, sorted_rv, sorted_resc_chisq):
        #print "Teff:", teff, "logg:", logg, "chi-squared:", dist_mod, "Radial_velocity:", rv, "rescaled chi-square:", resc_chisq
    
    print "best fit model:", "Teff", min_teff,"+/-", teff_sigma, "logg", min_logg, "+/-", logg_sigma, "chi-squared", min_dist, "Radial Velocity:", min_rv, 'km/s, rescaled chi-square:', min_rescale_dist
    #model_spec= get_model_fromfile(min_model)
    #if plot_all or plot_fit:
    #need to find the nearest model point to our thing
    nearest_teff= np.copy(teff_array[np.argmin(np.abs(teff_array.value-min_teff))]) #should return nearest teff
    nearest_logg= np.copy(logg_array[np.argmin(np.abs(logg_array- min_logg))])#should return nearest logg
    print "nearest model:", "Teff", nearest_teff, "logg", nearest_logg
    #min_model = wd(Teff= min_teff, logg = min_logg)
    min_model= wd(Teff=nearest_teff, logg= nearest_logg)
    
    dopp_cont_list= dopp_shift_continuum_list(min_rv)
    model_spec = np.vstack([min_model['w'], min_model['flux']])
    model_spec[0] = get_doppler_shifted(model_spec[0], min_rv)
    
    model_spec = mm.convolve_model(model_spec, target_spec, header)
    #normalization of the spectrum================
    model_spec= spt.poly_norm_spec(model_spec, continuum_list=dopp_cont_list, poly_degree = poly_degree, radial_velocity=min_rv)
    model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), mask_list)
    if balmer_only:
        temp_model_spec= spt.clean_spectrum(model_spec, np.min(target_spec[0]), np.max(target_spec[0]), continuum_masks)
        plot_overlays(line_spec, temp_model_spec, model_string = 'Teff ' + str(min_teff) + ' logg ' +str(min_logg)+ ' RV '+ str(min_rv)+'km/s')
    elif plot_fit:
        plot_overlays(target_spec, model_spec, model_string = 'Teff ' + str(min_teff) + ' logg ' +str(min_logg)+ ' RV '+ str(min_rv)+'km/s')
        
    if plot_all:
        interp_model= mm.interpolate_model(target_spec, model_spec)
        plot_overlays(target_spec,interp_model, model_string = 'interp Teff ' + str(min_teff) + ' logg ' +str(min_logg))
        plot_overlays(target_spec, noise_spec, model_string = 'noise')
        
        
        
        residual_spec= mm.calc_residuals(target_spec, model_spec)
        
        rms_spec= spt.clean_spectrum(residual_spec, rms_range[0], rms_range[1],mask_list)
        rms_scatter= np.sqrt(np.nanmean(rms_spec[1]**2))
        rms_noise_spec= spt.clean_spectrum(noise_spec, rms_range[0], rms_range[1],mask_list)
        mean_noise= np.nanmean(rms_noise_spec[1])
        rms_target= spt.clean_spectrum(target_spec, rms_range[0],rms_range[1],mask_list)
        rms_about_mean= np.sqrt(np.nanmean((rms_target[1]-np.nanmean(rms_target[1]))**2))
        rms_about_median= np.sqrt(np.nanmean((rms_target[1]-np.nanmedian(rms_target[1]))**2))
        print "RMS in", rms_range, ":", rms_scatter, "compared to model"
        print "Mean sigma of", rms_range,":", mean_noise
        print "RMS in", rms_range,":", rms_about_mean, "compared to mean:", np.nanmean(rms_target[1])
        print "RMS in",rms_range,":", rms_about_median, "compared to median:", np.nanmedian(rms_target[1])
        print "max noise in", rms_range, ":", np.nanmax(rms_noise_spec[1])
        rms_model= spt.clean_spectrum(interp_model, rms_range[0], rms_range[1], mask_list)
        stat_std_dev = np.sqrt(np.sum(rms_target[1]**2-rms_model[1]**2)/rms_target[1].shape[0])
        print "statistical std dev (around model):", stat_std_dev
        mean_std_dev= np.std(rms_target[1])
        print "statistical std dev (around mean):", mean_std_dev
    
    
        chi_sq_singles= residual_spec[1]**2/noise_spec[1]**2
        chi_sq_spec= np.vstack([target_spec[0], chi_sq_singles])
        print "Sum of residual chi-square:", np.sum(chi_sq_spec[1])
        print "chi-square values:", dist_array[min_index]
        plt.plot(residual_spec[0], residual_spec[1])
        plt.xlabel('Wavelength')
        plt.ylabel('data-model')
        plt.title('Residuals')
        plt.show()
        
        #plt.plot(chi_sq_spec[0], chi_sq_spec[1])
        #plt.xlabel('Wavelength')
        #plt.ylabel('chi-square')
        #plt.title('chi-square values by wavelength bin')
        #plt.show()
        f, (ax1, ax2,ax3)= plt.subplots(3,1, sharex=True)
        ax1.plot(target_spec[0], target_spec[1])
        ax1.plot(interp_model[0], interp_model[1])
        ax1.set_ylabel('Spectrum Flux')
        ax3.plot(chi_sq_spec[0], chi_sq_spec[1])
        ax3.set_ylabel('Chi-square')
        ax3.set_xlabel('Wavelength')
        ax2.plot(residual_spec[0], residual_spec[1])
        ax2.set_ylabel('residuals')
        ax2.axhline(y=0,color='k')
        f.subplots_adjust(wspace=0)
        f.subplots_adjust(hspace = 0)  
        plt.show()
        
        med_resid= np.nanmedian(residual_spec[1])
        mean_resid= np.nanmean(residual_spec[1])
        print "\n**************\nmedian residual:"+str(med_resid) + "\n***********"
        print "\n**************\nmean residual:"+str(mean_resid) + "\n***********"

        plt.hist(residual_spec[1], bins=20)
        plt.axvline(x=med_resid, color= 'r')
        plt.axvline(x=med_resid, color= 'g')
        plt.xlabel('data-model')
        plt.ylabel('N')
        plt.title('Residual Distribution')
        plt.show()
        
        med_chisq= np.nanmedian(chi_sq_spec[1])
        mean_chisq= np.nanmean(chi_sq_spec[1])
        print "\n**************\nmedian chi-square:"+str(med_chisq) + "\n***********"
        print "\n**************\nmean chi-square:"+str(mean_chisq) + "\n***********"
    
        plt.hist(chi_sq_spec[1], bins=100, normed=1)
        plt.axvline(x=np.nanmedian(chi_sq_spec[1]), color= 'r', label='median')
        plt.axvline(x=np.nanmean(chi_sq_spec[1]), color= 'g', label='mean')
        #plt.plot(chi2_linspace, chi2_dist.pdf(chi2_linspace), label = 'chi-square')
        plt.xlabel('chi-square')
        plt.ylabel('N')
        plt.title('chi-square distribution')
        plt.legend()
        plt.show()
        
        plt.plot(noise_spec[0], 1/noise_spec[1], color='b', label = '1/noise_spec')
        plt.plot(noise_spec[0], target_spec[1]/noise_spec[1], color = 'r', label= 'target_spec/noise_spec')
        plt.ylabel('S/N')
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.legend()
        plt.show()
        
        plt.errorbar(target_spec[0], target_spec[1], noise_spec[1], color = 'b', label = 'target')
        plt.plot(interp_model[0],interp_model[1], color = 'r', label = 'interp Teff ' + str(min_teff) + ' logg ' +str(min_logg))
        plt.ylabel('Flux (arbitrary units of /angstrom)')
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.show()
        
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
    else:
        pass
    return min_teff, teff_sigma, min_logg, logg_sigma
print get_doppler_shifted(4000, -200)

#target_spec = poly_norm_spec(target_spec)
######################################
#I can probably just iterate starting all the way down here.
if spectrum_type=='individual':
    print "spectrum_type==",spectrum_type
    print "So, iterating through astropy table:", astropy_input
    input_table = Table.read(astropy_input, format='ascii.csv')
    target_name_list=[]
    time_list= []
    teff_list=[]
    logg_list=[]
    teff_err_list=[]
    logg_err_list=[]
    rv_list=[]
    rv_err_list=[]
    baryv_list=[]
    
    for row in input_table:
        target_continuum_list= np.copy(continuum_list)
        target_continuum_list = get_doppler_shifted(target_continuum_list, row['rv'])
        target_spec, header, noise_spec = spt.retrieve_spec(row['filename'])
        #target_spec, header, noise_spec = spt.retrieve_spec(combined_spec_file)
        target_spec = spt.trim_spec(target_spec, low_wave_cut, high_wave_cut)
        #plt.plot(target_spec[0], target_spec[1])
        #plt.show()
        shift_mask=[]
        print "File:", row['filename']
        for mask in mask_list:
            new_mask = get_doppler_shifted(mask, row['rv'])
            shift_mask.append(new_mask)
        shift_mask = mask_list
        #print "shift_mask", shift_mask
        noise_spec = spt.trim_spec(noise_spec, low_wave_cut, high_wave_cut)
        #now undo the de-normalization of the noise spectrum that is done by retrieve_spec
        noise_spec[1]= noise_spec[1]/target_spec[1]
        target_spec= spt.clean_spectrum(target_spec, np.min(target_spec[0]), np.max(target_spec[0]), shift_mask)
        #plt.plot(target_spec[0], target_spec[1])
        #plt.show()
        noise_spec= spt.clean_spectrum(noise_spec, np.min(noise_spec[0]), np.max(noise_spec[0]), shift_mask)
        target_spec = spt.poly_norm_spec(target_spec, continuum_list=target_continuum_list, poly_degree = poly_degree, plot_all = plot_fit, radial_velocity=row['rv'])
        #plt.plot(target_spec[0], target_spec[1])
        #plt.show()
        noise_spec[1]= noise_spec[1]*target_spec[1] #scale the noise spectrum with the flattened target spectrum.
        min_teff, teff_sigma, min_logg, logg_sigma= run_model_grid(target_spec, velocity_tests= np.array([row['rv']]), noise_spec= noise_spec, plot_all=plot_fit, mask_list= shift_mask)
        bary_vcorr= barycentric_vel(header)
        baryv_list.append(bary_vcorr)
        target_name_list.append(row['filename'])
        time_list.append(header['BMJD_TDB'])
        logg_list.append(min_logg)
        logg_err_list.append(logg_sigma)
        try:
            teff_list.append(min_teff.value)
        except AttributeError:
            teff_list.append(min_teff)
        teff_err_list.append(teff_sigma)
        rv_list.append(row['rv'])
        rv_err_list.append(row['rv_error'])
    print "Hopefully saving the data"
    print target_name_list
    print time_list
    print logg_list
    print logg_err_list
    print teff_list
    print teff_err_list
    print rv_list
    print rv_err_list
    #array_list= [np.array(target_name_list), np.array(time_list), np.array(rv_list), np.array(rv_err_list), np.array(teff_list), np.array(teff_err_list), np.array(logg_list), np.array(logg_err_list)]
    #array_list= [target_name_list, time_list, rv_list, rv_err_list, teff_list, teff_err_list, logg_list, logg_err_list]
    array_list= [target_name_list, time_list, rv_list, rv_err_list, teff_list, teff_err_list, logg_list, logg_err_list, baryv_list]
    out_table= Table(array_list, names=col_names, dtype= dtype_list)
    out_table.pprint()
    out_table.write(output_astropy, format='ascii.csv')
    print "Saved the data"
else:
    run_model_grid(target_spec, noise_spec= noise_spec, plot_all= plot_fit)
