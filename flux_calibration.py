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


#standard_directory = '~/Desktop/standards/'
standard_directory = '/Users/BenKaiser/Desktop/standards/'
standard_file = "foke/fgd108.dat"

standard_file = standard_directory+standard_file

standard_name = "GD108"

observed_file = "wcmtb.GD108930blue.fits"
obs_fits = fits.open(observed_file)
header = fits.getheader(observed_file)
obs_waves = obs_fits[0].data
obs_flux = obs_fits[1].data


stand_array = np.genfromtxt(glob(standard_file)[0]).T

stand_waves = stand_array[0]
stand_flux = stand_array[1] *1e-16 #ergs/cm/cm/s/A 10**16 (That's exactly how it's written in the README, and it isn't -16, as one would assume...)

min_wave = np.nanmin(obs_waves)
max_wave = np.nanmax(obs_waves)

upper_cut = np.where(stand_waves < max_wave)
stand_waves = stand_waves[upper_cut]
stand_flux = stand_flux[upper_cut]
lower_cut = np.where(stand_waves > min_wave)
stand_waves = stand_waves[lower_cut]
stand_flux=stand_flux[lower_cut]

plt.plot(stand_waves, stand_flux,label = 'model')
plt.plot(obs_waves, obs_flux, label = 'observed')
plt.legend()
plt.show()

interp_model_flux = np.interp(obs_waves, stand_waves, stand_flux) #
plt.plot(obs_waves, interp_model_flux, label = 'interpolated')
plt.plot(stand_waves, stand_flux, label = 'model')
plt.legend()
plt.show()


plt.plot(obs_waves, interp_model_flux, label = 'interpolated')
plt.plot(stand_waves, stand_flux, label = 'model')
plt.legend()
plt.show()

interp_obs_flux = np.interp(stand_waves, obs_waves, obs_flux)

#sens_curve_points = obs_flux/interp_model_flux
#sens_curve_fit = np.polyfit(obs_waves, sens_curve_points, 5)
sens_curve_points= interp_obs_flux/stand_flux
sens_curve_fit= np.polyfit(stand_waves, sens_curve_points,5)


#poly_curve=  sens_curve_fit[-1]+sens_curve_fit[-2]*obs_waves + sens_curve_fit[-3]*(obs_waves**2)+sens_curve_fit[-4]*(obs_waves**3)+sens_curve_fit[-5]*(obs_waves**4)+sens_curve_fit[-6]*(obs_waves**5)
poly_curve = np.polyval(sens_curve_fit,obs_waves)

#plt.plot(obs_waves, sens_curve_points, label = 'data points', marker = 'o', linestyle = 'none')
#plt.plot(obs_waves, poly_curve, label= 'polynomial fit')
#plt.legend()
#plt.show()

plt.plot(stand_waves, sens_curve_points, label = 'data points', marker = 'o', linestyle = 'none')
plt.plot(obs_waves, poly_curve, label= 'polynomial fit')
plt.legend()
plt.show()

fcal_obs = obs_flux/poly_curve

plt.plot(obs_waves, fcal_obs, label ='flux calibrated observation')
plt.plot(stand_waves, stand_flux, label = 'model')
plt.legend()
plt.show()
