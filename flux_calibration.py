"""
This should produce the sensitivity curve that is needed for 
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


import spec_plot_tools as spt


#standard_directory = '~/Desktop/standards/'
standard_directory = '/Users/BenKaiser/Desktop/standards/'
#standard_file = "foke/fgd108.dat"
#standard_file = 'foke/ffeige67.dat'
standard_file = 'fhamuy/fltt6248.dat'
#standard_file= 'fhamuy/feg274.dat'


#output_filename = "GD108_sensitivity_curve.txt"
#output_filename = 'Feige67_sensitivity_curve.txt'
output_filename= 'LTT6248_sensitivity_curve.txt'
#output_filename= 'EG274_sensitivity_curve.txt'

standard_file = standard_directory+standard_file

#standard_name = "GD108"
#standard_name = 'Feige67'
standard_name = 'LTT6248'
#standard_name = 'EG274'

#observed_file = "wcmtb.GD108930blue.fits"
#observed_file = 'wcmtb.feige67930blue.fits'
observed_file  = 'wcmtb.ltt6248930blue.fits'
#observed_file = 'wcmtb.eg274930blue.fits'

obs_fits = fits.open(observed_file)
header = fits.getheader(observed_file)
obs_waves1= obs_fits[0].data
obs_flux1 = obs_fits[1].data
airmass = header['AIRMASS']
obs_time = header['OPENTIME']
obs_date = header['OPENDATE']
obs_time = obs_date+'T'+obs_time
obs_time = Time(obs_time, format = 'isot', scale = 'utc').mjd

#####
wavelength_masks=[
    [3792.92, 3811.62],
    [3823.59, 3853.88],
    [3867.34,3915.21],
    [3939.52, 4006.45],
    [4067.53, 4141.13],
    [4315.3, 4378.2],
    [4672.57,4706.4],
    [4835.18, 4907.76]
    ] #for Feige67

#####


stand_array = np.genfromtxt(glob(standard_file)[0]).T

stand_waves1 = stand_array[0]
#stand_flux1 = stand_array[1] *1e16 #ergs/cm/cm/s/A 10**16 (That's exactly how it's written in the README, and it isn't -16, as one would assume...)
stand_flux1 = stand_array[1]  #ergs/cm/cm/s/A 10**16 (That's exactly how it's written in the README, and it isn't -16, as one would assume...)

stand_bins = stand_array[3]

print stand_bins[0]

#plt.plot(stand_waves, stand_array[1]/np.mean(stand_array[1]), label='per angstrom flux')
#plt.plot(stand_waves, stand_array[2]/np.mean(stand_array[2]), label='jansky')
#plt.ylabel('flux divided by mean flux value')
#plt.legend()
#plt.show()

min_wave = np.nanmin(obs_waves1)
max_wave = np.nanmax(obs_waves1)

#upper_cut = np.where(stand_waves < max_wave)
#stand_waves = stand_waves[upper_cut]
#stand_flux = stand_flux[upper_cut]
#lower_cut = np.where(stand_waves > min_wave)
#stand_waves = stand_waves[lower_cut]
#stand_flux=stand_flux[lower_cut]




plt.title('model versus observed')
plt.plot(stand_waves1, stand_flux1,label = 'model')
plt.plot(obs_waves1, obs_flux1, label = 'observed')
plt.legend()
plt.show()

plt.title('interpolated model versus standard model')
interp_model_flux = np.interp(obs_waves1, stand_waves1, stand_flux1) #
plt.plot(obs_waves1, interp_model_flux, label = 'interpolated')
plt.plot(stand_waves1, stand_flux1, label = 'model')
plt.legend()
plt.show()

obs_spec = np.vstack([obs_waves1, obs_flux1])
stand_spec= np.vstack([stand_waves1, stand_flux1])

obs_spec = spt.clean_spectrum(obs_spec, min_wave, max_wave, wavelength_masks)
stand_spec =spt.clean_spectrum(stand_spec, min_wave, max_wave, wavelength_masks)
stand_waves = stand_spec[0]
stand_flux= stand_spec[1]
obs_waves = obs_spec[0]
obs_flux = obs_spec[1]

plt.title('model versus observed')
plt.plot(stand_waves, stand_flux,label = 'model')
plt.plot(obs_waves, obs_flux, label = 'observed')
plt.legend()
plt.show()


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

obs_curve= np.polyfit(obs_waves, obs_flux, 5)
model_curve = np.polyfit(stand_waves, stand_flux, 5)

sens_curve_points = np.polyval(obs_curve, stand_waves)/np.polyval(model_curve, stand_waves)
sens_curve_fit= np.polyfit(stand_waves, sens_curve_points,5)

plt.plot(obs_waves, obs_flux, label= 'observed', marker= 'o', linestyle='none')
plt.plot(obs_waves, np.polyval(obs_curve, obs_waves), label = 'curve')
plt.legend()
plt.show()

plt.plot(stand_waves, stand_flux, label= 'model', marker= 'o', linestyle = 'none')
plt.plot(stand_waves, np.polyval(model_curve, stand_waves), label = 'curve')
plt.legend()
plt.show()




np.savetxt(output_filename, sens_curve_fit, header = 'Airmass: ' +str(airmass) + '\tMJD: ' +str(obs_time))
#poly_curve=  sens_curve_fit[-1]+sens_curve_fit[-2]*obs_waves + sens_curve_fit[-3]*(obs_waves**2)+sens_curve_fit[-4]*(obs_waves**3)+sens_curve_fit[-5]*(obs_waves**4)+sens_curve_fit[-6]*(obs_waves**5)
poly_curve = np.polyval(sens_curve_fit,obs_waves1)

#plt.plot(obs_waves, sens_curve_points, label = 'data points', marker = 'o', linestyle = 'none')
#plt.plot(obs_waves, poly_curve, label= 'polynomial fit')
#plt.legend()
#plt.show()

plt.plot(stand_waves, sens_curve_points, label = 'data points', marker = 'o', linestyle = 'none')
plt.plot(obs_waves1, poly_curve, label= 'polynomial fit')
plt.legend()
plt.show()

fcal_obs = obs_flux1/poly_curve

#plt.plot(obs_waves1, fcal_obs, label ='flux calibrated observation', marker= 'o', linestyle = 'none')
#plt.plot(stand_waves1, stand_flux1, label = 'model', marker = 'o', linestyle = 'none')
plt.plot(obs_waves1, fcal_obs, label ='flux calibrated observation')
plt.plot(stand_waves1, stand_flux1, label = 'model')
plt.xlabel('wavelength ($\AA$)')
plt.ylabel('Flux (ergs/cm/cm/s/A 1e-16)')
plt.legend()
plt.show()
