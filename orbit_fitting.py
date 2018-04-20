import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
from astropy.units import cds
import scipy.optimize as sciop
cds.enable()


plotfile = glob(sys.argv[1])
print plotfile
all_array = np.genfromtxt(plotfile[0], names = True, delimiter = ',')

parkes_location = coord.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
cerro_pachon_location = coord.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)

#e = 0.000023# median
#e = 0.000031# max
omega = 97 * u.degree
period = 0.4497391377  #median
#period = (0.449739137 * u.day).to(u.second) #low
#period = (0.4497391384 * u.day).to(u.second) #high

t0 = Time(55756.23, format = 'mjd', scale= 'utc', location = parkes_location)
tasc = Time(55756.1047771, format = 'mjd', scale= 'utc', location = parkes_location) #median
#tasc = Time(55756.1047767, format = 'mjd', scale= 'utc', location = parkes_location) #min
tconj= Time(55756.21712, format = 'mjd', scale ='utc', location = parkes_location)

#epoch_difference = obs_times[0].mjd-tasc.mjd


bmjd_array = all_array['BMJD_TDB']
H_delta = all_array['H_delta']
H_gamma = all_array["H_gamma"]
H_beta = all_array['H_beta']
#print Time(mjd_array, format = 'mjd').utc.isot
#print all_array

p0_list = [-100, period, 300,bmjd_array[0]-0.1]



#bmjd_obs = to_barycenter(obs_times) #corrected to barycenter to use against the rv curve
#bmjd_t0 = to_barycenter(t0) #corrected initial epoch
#bmjd_tasc= to_barycenter(tasc)
#bmjd_tconj = to_barycenter(tconj)

#zero_point = bmjd_tconj

def sine_function(times, systemic_vel, period, amplitude, zero_point):
    return amplitude*np.sin((2*np.pi)/period*(times-zero_point))+systemic_vel
fitted_curve_delta, fitted_curve_cov = sciop.curve_fit(sine_function, bmjd_array, H_delta, p0_list)
fitted_curve_gamma, fitted_curve_cov = sciop.curve_fit(sine_function, bmjd_array, H_gamma, p0_list)
fitted_curve_beta, fitted_curve_cov = sciop.curve_fit(sine_function, bmjd_array, H_beta, p0_list)

print "sys vel.   period    amplitude    zero_point"
print fitted_curve_delta
print fitted_curve_gamma
print fitted_curve_beta

def zero_rvs(rv_array):
    print "systemic velocity (includes Earth's motion):", np.mean([rv_array.max(),rv_array.min()])
    return rv_array-np.mean([rv_array.max(),rv_array.min()])

#H_delta = zero_rvs(H_delta)
#H_gamma = zero_rvs(H_gamma)
#H_beta = zero_rvs(H_beta)
#remean_rv = np.mean([H_delta, H_gamma, H_beta], axis = 0)
#remean_std = np.std([H_delta,H_gamma, H_beta], axis=0)
#mean_rv = zero_rvs(mean_rv)
plt.plot(bmjd_array, H_delta, label = r"H-$\delta$", linestyle = 'none', marker = '*')
plt.plot(bmjd_array, H_gamma, label = r"H-$\gamma$", linestyle = 'none', marker = '*')
plt.plot(bmjd_array, H_beta, label = r"H-$\beta$", linestyle = 'none', marker = '*')

plot_times = np.linspace(bmjd_array[0]-0.2, bmjd_array[-1]+0.2, 1000)
plt.plot(plot_times, sine_function(plot_times, fitted_curve_delta[0],fitted_curve_delta[1],fitted_curve_delta[2], fitted_curve_delta[3]), label = 'fitted curve delta')
plt.plot(plot_times, sine_function(plot_times, fitted_curve_gamma[0],fitted_curve_gamma[1],fitted_curve_gamma[2], fitted_curve_gamma[3]), label = 'fitted curve gamma')
plt.plot(plot_times, sine_function(plot_times, fitted_curve_beta[0],fitted_curve_beta[1],fitted_curve_beta[2], fitted_curve_beta[3]), label = 'fitted curve beta')
plt.legend()
plt.title("Period: " + str(fitted_curve_beta[1] )+ " Sys Vel.: " + str(fitted_curve_beta[0]) + " K: " + str(fitted_curve_beta[2]))
#plt.errorbar(mjd_array, mean_rv, std_dev, label = 'Mean RV', linestyle = 'none', marker = 'o')
#plt.errorbar(mjd_array, remean_rv, remean_std, label = r"Mean of zeroed RV's", linestyle = 'none', marker = 'o')
plt.show()
