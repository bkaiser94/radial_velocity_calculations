"""
Created by Ben Kaiser (UNC-Chapel Hill) (date not known of original creation.)
This should produce the sensitivity curve that is needed for calibrate_flux.py.

Uses the outputs of actual_reduction.py.

Have to manually change targets that are being used... probably could automate that to actually load all of the
necessary options after at least changing targets... doesn't seem really worth it at this point.

Step 5 of Reduction
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
import scipy.interpolate as scinterp


import spec_plot_tools as spt
import cal_params as cp

#poly_degree=3
#poly_degree= 5 #order of polynomials from before 20190506
#poly_degree=7
poly_degree=7
model_poly_degree= 5

ext_corr=True #correct extinction


#standard_directory = '~/Desktop/standards/'
#standard_directory = '/Users/BenKaiser/Desktop/standards/'
standard_directory= cp.standard_dir
#standard_file = "foke/fgd108.dat"
#standard_file = 'foke/ffeige67.dat'
#standard_file = 'fhamuy/fltt6248.dat'
#standard_file= 'fhamuy/feg274.dat'
#standard_file= 'ctio_standards/fltt3218.dat'


#output_filename = "GD108_sensitivity_curve.txt"
#output_filename = 'Feige67_sensitivity_curve.txt'
#output_filename= 'LTT6248_sensitivity_curve.txt'
#output_filename= 'EG274_sensitivity_curve.txt'
#output_filename = 'LTT3218_sensitivity_curve.txt'

#standard_file = standard_directory+standard_file

##standard_name = "GD108"
#standard_name = 'Feige67'
#standard_name = 'LTT6248'
standard_name='EG274'
#standard_name = 'GD153'
#standard_name= 'LTT3218'
#standard_name='Feige110'

##observed_file = "wcmtb.GD108930blue.fits"
##observed_file = 'wcmtb.feige67930blue.fits'
#observed_file  = 'wcmtb.ltt6248930blue.fits'
##observed_file = 'wcmtb.eg274930blue.fits'
##observed_file= 'wcmtb.ltt3218930blue.fits'


#observed_file = "wcmtb.GD108930blue.fits"
#observed_file = 'wctb.0232_feige67_930_blue.fits'
#observed_file  = 'wctb.0244_LTT6248_930_blue.fits'
#observed_file  = 'wctb.0269_LTT6248_930_blue_1arcsec.fits'
#observed_file = 'wctb.0272_eg274_930_blue.fits'
#observed_file= 'wcmtb.ltt3218930blue.fits'

#observed_file='avg_EG274_400m1.fits'
#observed_file='avg_wctb.EG274_400m2.fits'
#observed_file='avg_wctb.Feige110_400m2.fits'
#observed_file='avg_wctb.Feige110_400m1.fits'

#observed_file='avg_wctb.GD153_400m2.fits'
#observed_file='avg_wctb.EG274_400m1.fits'
#observed_file='avg_wctb.EG274_400m2.fits'
observed_file='avg_wctb.eg274_930_blue.fits'
#observed_file='avg_wctb.eg274am104_400m2.fits'
#observed_file='avg_wctb.eg274am126_400m2.fits'
#observed_file='avg_wctb.eg274am176_400m2.fits'
#observed_file='avg_wctb.GD153_400m1.fits'
#observed_file='avg_wctb.EG274_400m1_am13.fits'
#observed_file='avg_wctb.EG274_400m2_am13.fits'
##############################
##############################

def get_star_info(starname):
    standard_dict= cp.standard_dict[starname.lower()]
    standard_dict['filename']=standard_directory+standard_dict['filename']
    return standard_dict










############################
##########################


core_name= observed_file.split('.')[1] #get the part of the filename that follows the first period and exclude the extension
output_filename=core_name+'_sensitivity_curve.txt'


obs_fits = fits.open(observed_file)
header = fits.getheader(observed_file)
obs_waves1= obs_fits[0].data
obs_flux1 = obs_fits[1].data
airmass = header['AIRMASS']
obs_time = header['OPENTIME']
obs_date = header['OPENDATE']
obs_time = obs_date+'T'+obs_time
obs_time = Time(obs_time, format = 'isot', scale = 'utc').mjd
exptime = header['EXPTIME']

obs_spec= np.vstack([obs_waves1, obs_flux1])
try:
    dlambda= obs_fits[4].data
    obs_spec= np.copy(spt.counts_to_flambda(obs_spec, dlambda))
    print('Obseved spectrum in units of erg/s/cm^2/angstrom')
except IndexError:
    print('No dlambda extension in observed FITS file. You need to redo wave_cal.py to include that extension.')
    print('This file is most likely generated before 2019-07-16 when this change was implemented')

if ext_corr:
    obs_spec= spt.correct_extinction(obs_spec, header, plot_all=True)
else:
    pass
obs_waves1=obs_spec[0]
obs_flux1=obs_spec[1]
#obs_flux1= obs_flux1/np.float_(exptime) #converts to counts per second

#####
#wavelength_masks=[
    #[3792.92, 3811.62],
    #[3823.59, 3853.88],
    #[3867.34,3915.21],
    #[3939.52, 4006.45],
    #[4067.53, 4141.13],
    #[4315.3, 4378.2],
    #[4672.57,4706.4],
    #[4835.18, 4907.76]
    #] #for Feige67

#wavelength_masks=[
    #[3792.92, 3811.62],
    #[3823.59, 3853.88],
    #[3867.34,3915.21],
    #[3939.52, 4029.45],
    #[4046.53, 4189.13],
    #[4251.3, 4470.2],
    #[4672.57,4706.4],
    #[4781.18, 4994.76]
    #] #for EG274
#####

#standard_file = standard_directory+standard_file
standard_info = get_star_info(standard_name)
print type(standard_info['balmer_masks'])
print type(standard_info['other_masks'])
wavelength_masks=standard_info['balmer_masks']+standard_info['other_masks']
print "wavelength_masks:", wavelength_masks
#stand_array = np.genfromtxt(glob(standard_file)[0]).T
#print(glob(standard_info['filename']))
stand_array = np.genfromtxt(glob(standard_info['filename'])[0]).T
#output_filename= standard_dict['sens_filename']

stand_waves1 = stand_array[0]
#stand_flux1 = stand_array[1] *1e16 #ergs/cm/cm/s/A 10**16 (That's exactly how it's written in the README, and it isn't -16, as one would assume...)
#stand_flux1 = stand_array[1]  #ergs/cm/cm/s/A 10**16 (That's exactly how it's written in the README, and it isn't -16, as one would assume...)

stand_flux1 = stand_array[1]  #ergs/cm/cm/s/A (That's exactly how it's written in the README for X-shooter)


def rescale_flux(stand_flux1):
    """
    Make sure the flux is in the correct units for your purposes, which should be 1e-16
    
    """
    standard_type= standard_info['filename'].split('/')[-2]
    print('standard_type:', standard_type)
    if standard_type == 'xshooter_standards':
        print('should be resetting fluxes')
        stand_flux1= stand_flux1*1e16 #converts to 10**-16 flux vals hopefully
    else:
        stand_flux1= stand_flux1
    return stand_flux1

#stand_flux1= rescale_flux(stand_flux1)

#stand_bins = stand_array[3]

#print stand_bins[0]

#plt.plot(stand_waves, stand_array[1]/np.mean(stand_array[1]), label='per angstrom flux')
#plt.plot(stand_waves, stand_array[2]/np.mean(stand_array[2]), label='jansky')
#plt.ylabel('flux divided by mean flux value')
#plt.legend()
#plt.show()

#min_wave = np.nanmin(obs_waves1)
#max_wave = np.nanmax(obs_waves1)

#upper_cut = np.where(stand_waves < max_wave)
#stand_waves = stand_waves[upper_cut]
#stand_flux = stand_flux[upper_cut]
#lower_cut = np.where(stand_waves > min_wave)
#stand_waves = stand_waves[lower_cut]
#stand_flux=stand_flux[lower_cut]




plt.title('model versus observed')
plt.plot(stand_waves1, stand_flux1/np.nanmedian(stand_flux1),label = 'model')
plt.plot(obs_waves1, obs_flux1/np.nanmedian(obs_flux1), label = 'observed')
#plt.scatter(stand_waves1, stand_flux1,label = 'model')
#plt.scatter(obs_waves1, obs_flux1, label = 'observed', color='r')
#plt.legend()
#plt.show()
plt.xlabel('Wavelengths (Angstroms)')
plt.ylabel('Flux(erg/s/cm^2/A)')
spt.show_plot()

do_offset= bool(raw_input("Do you need to do a wavelength offset?(True/False)>>>"))
if do_offset:
    print "Enter the approximate wavelength for the same feature in the model and observed spectra for offset"
    model_wavelength = float(raw_input("Model spec wavelength>>>"))
    obs_wavelength= float(raw_input("Observed spec wavelength>>>"))
    #dotted_pixel=0
    #emission_pixel=0
    offset = model_wavelength-obs_wavelength
    obs_waves1=obs_waves1+offset
    plt.title('model versus observed')
    plt.plot(stand_waves1, stand_flux1/np.nanmedian(stand_flux1),label = 'model')
    plt.plot(obs_waves1, obs_flux1/np.nanmedian(obs_flux1), label = 'observed')
    plt.xlabel('Wavelengths (Angstroms)')
    plt.ylabel('Flux(erg/s/cm^2/A)')
    #plt.legend()
    #plt.show()
    spt.show_plot()
else:
    print "Skipping offsetting"
    

min_wave = np.nanmin(obs_waves1)
max_wave = np.nanmax(obs_waves1)


plt.title('interpolated model versus standard model')
#interp_model_flux = np.interp(obs_waves1, stand_waves1, stand_flux1) #
interpolator = scinterp.CubicSpline(stand_waves1, stand_flux1)
interp_model_flux= interpolator(obs_waves1)
plt.plot(obs_waves1, interp_model_flux, label = 'interpolated')
plt.plot(stand_waves1, stand_flux1, label = 'model')
plt.legend()
plt.show()

obs_spec = np.vstack([obs_waves1, obs_flux1])
#stand_spec= np.vstack([stand_waves1, stand_flux1])
stand_spec= np.vstack([obs_waves1, interp_model_flux ])


obs_spec = spt.clean_spectrum(obs_spec, min_wave, max_wave, wavelength_masks)
unmasked_stand_spec=spt.clean_spectrum(stand_spec, min_wave, max_wave,standard_info['balmer_masks'])
stand_spec =spt.clean_spectrum(stand_spec, min_wave, max_wave, wavelength_masks)
stand_waves = stand_spec[0]
stand_flux= stand_spec[1]
obs_waves = obs_spec[0]
obs_flux = obs_spec[1]



plt.title('model versus observed')
plt.plot(stand_waves, stand_flux,label = 'model')
plt.plot(obs_waves, obs_flux, label = 'observed')
#plt.legend()
#plt.show()
spt.show_plot()

#plt.title('wavelength values')
#plt.plot(obs_waves,np.ones(obs_waves.shape), label = 'observed', marker='o')
#plt.plot(stand_waves, np.zeros(stand_waves.shape), label = 'model', marker= 'o')
#plt.legend()
#plt.show()


#interp_obs_flux = np.interp(stand_waves, obs_waves, obs_flux)

#sens_curve_points = obs_flux/interp_model_flux
#sens_curve_fit = np.polyfit(obs_waves, sens_curve_points, 5)
#sens_curve_points= interp_obs_flux/stand_flux
#sens_curve_fit= np.polyfit(stand_waves, sens_curve_points,5)

obs_curve= np.polyfit(obs_waves, obs_flux, poly_degree)
#model_curve = np.polyfit(stand_waves, stand_flux, poly_degree)
#model_curve = np.polyfit(stand_waves, stand_flux, model_poly_degree)
model_curve = np.polyfit(unmasked_stand_spec[0], unmasked_stand_spec[1], model_poly_degree)

calc_waves=np.linspace(min_wave, max_wave,1000)
#sens_curve_points = np.polyval(obs_curve, stand_waves)/np.polyval(model_curve, stand_waves)
#sens_curve_points = np.polyval(obs_curve, obs_waves)/np.polyval(model_curve, obs_waves)
sens_curve_points= np.polyval(obs_curve,calc_waves)/np.polyval(model_curve, calc_waves)
#sens_curve_fit= np.polyfit(stand_waves, sens_curve_points,5)
#sens_curve_fit= np.polyfit(obs_waves, sens_curve_points,poly_degree)
sens_curve_fit= np.polyfit(calc_waves, sens_curve_points,poly_degree+2)

#sens_curve_fit=np.polyfit(obs_waves, obs_flux/stand_flux,poly_degree) #20190618


plt.plot(obs_waves1, obs_flux1, label= 'observed', marker= 'o', linestyle='none')
plt.plot(obs_waves, obs_flux, label= 'observed used', marker= 'o', linestyle='none')
plt.plot(obs_waves1, np.polyval(obs_curve, obs_waves1), label = 'curve')
plt.legend()
plt.show()

plt.plot(obs_waves1, interp_model_flux, label= 'model', marker= 'o', linestyle = 'none')
plt.plot(stand_waves, stand_flux, label= 'model used', marker= 'o', linestyle = 'none')
plt.plot(obs_waves1, np.polyval(model_curve, obs_waves1), label = 'curve')
plt.legend()
plt.show()



poly_curve = np.polyval(sens_curve_fit,obs_waves1)

#plt.plot(stand_waves, sens_curve_points, label = 'data points', marker = 'o', linestyle = 'none')
plt.plot(obs_waves1, obs_flux1/interp_model_flux, label='all obs/model', marker='o', linestyle='none')
#plt.plot(obs_waves, sens_curve_points, label = 'used obs/model', marker = 'o', linestyle = 'none')
plt.plot(calc_waves, sens_curve_points, label = 'used obs curve/model curve', marker = 'o', linestyle = 'none')
#plt.plot(calc_waves, sens_curve_points, label = 'data points', marker = 'o', linestyle = 'none')
#plt.plot(obs_waves1, poly_curve, label= 'polynomial fit')
plt.plot(obs_waves1, poly_curve, label= 'polynomial fit')
plt.legend()
plt.show()

#fcal_obs = obs_flux1/poly_curve
fcal_obs = obs_flux1/poly_curve


#plt.plot(obs_waves1, fcal_obs, label ='flux calibrated observation', marker= 'o', linestyle = 'none')
#plt.plot(stand_waves1, stand_flux1, label = 'model', marker = 'o', linestyle = 'none')
#plt.plot(obs_waves1, fcal_obs, label ='flux calibrated observation')
#plt.plot(stand_waves1, stand_flux1, label = 'model')
plt.plot(obs_waves1, fcal_obs, label ='flux calibrated observation')
#plt.plot(stand_waves1, stand_flux1, label = 'model')
plt.plot(obs_waves1, interp_model_flux, label='model')
plt.xlabel('wavelength ($\AA$)')
plt.ylabel('Flux (ergs/cm/cm/s/A 1e-16)')
#plt.legend()
#plt.show()
spt.show_plot()

def get_residuals(plot_all = False):
    residuals= fcal_obs/interp_model_flux
    if plot_all:
        plt.plot(obs_waves1, residuals)
        plt.xlabel('wavelength ($\AA$)')
        #plt.ylabel('Flux (ergs/cm/cm/s/A 1e-16)')
        plt.ylabel('F_obscal/ F_model')
        plt.title('residuals')
        plt.axhline(y=1, color='k')
        plt.xlim(np.nanmin(obs_waves1), np.nanmax(obs_waves1))
        #plt.show()
        spt.show_plot()
    return residuals

def limit_to_telluric(obs_waves1,residuals):
    """
    Take the calculated residuals and limit them to the region of telluric absorption and set the value to be 1 everywhere else, so the remnants of other regions aren't messed up too.
    
    """
    #io_telluric_lines= 
    
    return residuals

residuals= get_residuals(plot_all=True)

plt.title('residual corrected spectrum')
plt.plot(obs_waves1, fcal_obs/residuals, label='flux-calibrated observation/residuals')
plt.plot(obs_waves1, interp_model_flux, label='model')
plt.xlabel('wavelength ($\AA$)')
plt.ylabel('Flux (ergs/cm/cm/s/A 1e-16)')
#plt.legend()
#plt.show()
spt.show_plot()

#np.savetxt(output_filename, sens_curve_fit, header = 'Airmass: ' +str(airmass) + '\tMJD: ' +str(obs_time))
#np.savetxt(standard_info['sens_filename'], sens_curve_fit, header = 'Airmass: ' +str(airmass) + '\tMJD: ' +str(obs_time))

#np.savetxt('residuals_' + standard_info['sens_filename'], residuals, header='Airmass: ' +str(airmass) + '\tMJD: ' +str(obs_time))

np.savetxt(output_filename, sens_curve_fit, header = 'Airmass: ' +str(airmass) + '\tMJD: ' +str(obs_time))

np.savetxt('residuals_' +output_filename, residuals, header='Airmass: ' +str(airmass) + '\tMJD: ' +str(obs_time))




