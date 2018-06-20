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
import scipy.stats as scistats
import scipy.optimize as sciop
cds.enable()

#p0_list = [-100, 300, 0]
p0_list = [-100, 300]

precision = 3

parkes_location = coord.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
cerro_pachon_location = coord.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)
c= 2.998E8 *u.m/u.s

#open one of the original files and get the RA and Dec from it
listnames = np.genfromtxt('listFWCTB', dtype ='str')
filename = listnames[0]
print filename
header = fits.getheader(filename)
ra = header['RA']
dec = header['DEC']
target_coord = coord.SkyCoord(ra, dec, frame = 'icrs', unit= (u.hourangle, u.deg))




def to_barycenter(input_times):
    bary_corr =input_times.tdb.light_travel_time(target_coord)
    return (input_times.tdb+ bary_corr.tdb).mjd

#def sine_function(times, systemic_vel,  amplitude, phase):
#def sine_function(times, systemic_vel,  amplitude, phase, period):
def sine_function(times, systemic_vel,  amplitude):
    #return amplitude*np.sin((2*np.pi)*(times)+phase)+systemic_vel
    return amplitude*np.sin((2*np.pi)*(times))+systemic_vel



#e = 0.000023# median
#e = 0.000031# max
omega = 97 * u.degree
period = 0.4497391377  #median
#period = (0.449739137 * u.day).to(u.second) #low
#period = (0.4497391384 * u.day).to(u.second) #high

x_psr = 0.550061*u.second #the semi-major axis of the pulsar's orbit (projected)


t0 = Time(55756.23, format = 'mjd', scale= 'utc', location = parkes_location)
tasc = Time(55756.1047771, format = 'mjd', scale= 'utc', location = parkes_location) #median
#tasc = Time(55756.1047767, format = 'mjd', scale= 'utc', location = parkes_location) #min
tconj= Time(55756.21712, format = 'mjd', scale ='utc', location = parkes_location)


all_array = np.genfromtxt('rv_plot.txt', names = True, delimiter= ',')
bmjd_array = all_array['TimesBMJD_TDB']
rv_array = all_array['RV_kms']


#all_arrayB = np.genfromtxt('rv_plotB.txt', names = True, delimiter= ',')
#bmjd_arrayB = all_arrayB['TimesBMJD_TDB']
#rv_arrayB = all_arrayB['RV_kms']


#bmjd_obs = to_barycenter(obs_times) #corrected to barycenter to use against the rv curve
bmjd_t0 = to_barycenter(t0) #corrected initial epoch
bmjd_tasc= to_barycenter(tasc)
bmjd_tconj = to_barycenter(tconj)

zero_point = bmjd_tconj
#zero_point = bmjd_tasc



zero_times = bmjd_array - zero_point
folded_times = np.mod(zero_times , period)/period
#folded_timesB = np.mod(bmjd_arrayB-zero_point, period)/period


fitted_curve_all, fitted_cov_all = sciop.curve_fit(sine_function, folded_times, rv_array, p0= p0_list)
residuals = rv_array- sine_function(folded_times, fitted_curve_all[0], fitted_curve_all[1])
print fitted_curve_all


def get_mass_ratio(K_c, P_B, x_psr):
    """
    Mass ratio M_psr/M_companion
    """
    print (K_c.to(u.m/u.s)*P_B.to(u.s)/(2*np.pi*x_psr.to(u.s)))
    return K_c.to(u.m/u.s)*P_B.to(u.s)/(2*np.pi*x_psr.to(u.s)*c)

period_wunits = period*u.day

all_mratio= get_mass_ratio(fitted_curve_all[1]*(u.km/u.second),period_wunits, x_psr)

x_vals = np.linspace(0,1,100)
#sine_vals = sine_function(x_vals, fitted_curve_all[0], fitted_curve_all[1], fitted_curve_all[2])
sine_vals = sine_function(x_vals, fitted_curve_all[0], fitted_curve_all[1])



fig = plt.figure()
ax = fig.add_subplot(111)

ax.scatter(folded_times, rv_array, color = 'b', label = '3/18/18')
#ax.scatter(folded_timesB, rv_arrayB, color = 'r', label = '2/13/18')
ax.plot(x_vals, sine_vals, color = 'k', linestyle = '--', label = 'Model')

ax.set_xlim([0,1])
ax.set_xlabel('Phase')
ax.set_ylabel('RV (km/s)')
ax.set_title(str(np.round(1.4/all_mratio,precision)) + ' M_sun companion assuming 1.4Msun NS')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
text_string =  r'$K_c =$ ' +str(np.round(fitted_curve_all[1], precision)) + r'$\pm$ ' + str(np.round(np.sqrt(fitted_cov_all[1,1]),precision)) + 'km/s' +'\n'+ r'$v_{sys}=$' +str(np.round(fitted_curve_all[0],precision)) + r'$\pm$' + str(np.round(np.sqrt(fitted_cov_all[0,0]),precision))+ ' km/s'
ax.text(0.05, 0.05,text_string, transform=ax.transAxes, fontsize=14, verticalalignment='bottom', bbox=props)
#ax.legend()

plt.show()


#######
#fig = plt.figure()
#ax1 = fig.add_subplot(121)

#ax1.scatter(folded_times, rv_array, color = 'b', label = '3/18/18')
##ax.scatter(folded_timesB, rv_arrayB, color = 'r', label = '2/13/18')
#ax1.plot(x_vals, sine_vals, color = 'k', linestyle = '--', label = 'Model')

#ax1.set_xlim([0,1])
#ax1.set_xlabel('Phase')
#ax1.set_ylabel('RV (km/s)')
#ax1.set_title(str(np.round(1.4/all_mratio,precision)) + ' M_sun companion assuming 1.4Msun NS')
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#text_string =  r'$K_c =$ ' +str(np.round(fitted_curve_all[1], precision)) + r'$\pm$ ' + str(np.round(np.sqrt(fitted_cov_all[1,1]),precision)) + 'km/s' +'\n'+ r'$v_{sys}=$' +str(np.round(fitted_curve_all[0],precision)) + r'$\pm$' + str(np.round(np.sqrt(fitted_cov_all[0,0]),precision))+ ' km/s'
#ax1.text(0.05, 0.05,text_string, transform=ax.transAxes, fontsize=14, verticalalignment='bottom', bbox=props)
##ax.legend()

#plt.show()
#fig = plt.figure(figsize=(8, 6)) 
#gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
#ax0 = plt.subplot(gs[0])
#ax0.plot(x, y)
#ax1 = plt.subplot(gs[1])
#ax1.plot(y, x)

#plt.tight_layout()
#plt.savefig('grid_figure.pdf')
