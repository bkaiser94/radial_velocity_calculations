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
plt.rc('font', size =18)
#plt.rc('lines', markersize=12)
#plt.rc('font', size = 11)
plt.rc('lines', markersize = 5)
plotting_offset = 0.0005
#p0_list = [-100, 300, 0]
p0_list = [-100, 300]
photo_bounds = ([-np.inf, 0, -np.inf],[np.inf, np.inf, np.inf])
precision = 1
precision2= 4
parkes_location = coord.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
cerro_pachon_location = coord.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)
c= 2.998E8 *u.m/u.s

output_filename = "model_rvs.txt"

#open one of the original files and get the RA and Dec from it
listnames = np.genfromtxt('listFWCTB', dtype ='str')
filename = listnames[0]
print filename
header = fits.getheader(filename)
ra = header['RA']
dec = header['DEC']
target_coord = coord.SkyCoord(ra, dec, frame = 'icrs', unit= (u.hourangle, u.deg))
#rv_file = 'rv_plot.txt'
rv_file = 'rv_plot_filled_in.txt'
photometry_t0 = 2458231.5950237 #BJD_TDB, so need to convert to MJD
photometry_t0= Time(photometry_t0, format = 'jd', scale = 'tdb')
photometry_file = 'psrj1431m4715_lightcurve.dat'
photometry_all = np.genfromtxt(photometry_file, skip_header = 2).T
#print photometry_all

photometry_times = photometry_all[0]*u.s
photometry_flux = photometry_all[1]
photometry_error = photometry_all[2]

photometry_times = (photometry_times + photometry_t0).tdb.mjd
###################################


def to_barycenter(input_times):
    bary_corr =input_times.tdb.light_travel_time(target_coord)
    return (input_times.tdb+ bary_corr.tdb).mjd

#def sine_function(times, systemic_vel,  amplitude, phase):
#def sine_function(times, systemic_vel,  amplitude, phase, period):
def sine_function(times, systemic_vel,  amplitude):
    #return amplitude*np.sin((2*np.pi)*(times)+phase)+systemic_vel
    return amplitude*np.sin((2*np.pi)*(times))+systemic_vel

def photo_sine_function(times, systemic_vel,  amplitude, phase):
#def photo_sine_function(times, systemic_vel,  amplitude):
    return amplitude*np.sin(2*(2*np.pi*times+phase))+systemic_vel
    #return amplitude*np.sin(2*(2*np.pi) *(times+1.45267158))+systemic_vel
    #return amplitude*np.cos(2*(2*np.pi)*(times))+systemic_vel
    
    
#def photo_harmonic_function(times, zeropoint, A, B, C, D, phiA, phiB, phiC, phiD):
#def photo_harmonic_function(times, zeropoint, A, B, C, D):
def photo_harmonic_function(times, zeropoint, B, D):

    #return A*np.sin(2*np.pi*times+phiA) + B*np.cos(2*np.pi*times+phiB) + C*np.sin(4*np.pi*times+phiC) + D*np.cos(4*np.pi*times+phiD) + zeropoint
    #return A*np.sin(2*np.pi*times) + B*np.cos(2*np.pi*times) + C*np.sin(4*np.pi*times) + D*np.cos(4*np.pi*times) + zeropoint
    return  B*np.cos(2*np.pi*times)+ D*np.cos(4*np.pi*times) + zeropoint



def make_value_strings(fitted_coeffs, fitted_cov):
    diagonals = []
    for i in range(fitted_cov.shape[0]):
        diagonals.append(fitted_cov[i,i])
    diagonals = np.sqrt(np.array(diagonals))
    string_list = []
    for coeff, cov in zip(fitted_coeffs, diagonals):
        new_string = r'$($'+str(np.round(coeff, precision2)) + r'$\pm$' + str(np.round(cov, precision2))+r'$)$'
        string_list.append(new_string)
    print string_list
    return string_list


#############################################
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


all_array = np.genfromtxt(rv_file, names = True, delimiter= ',')
bmjd_array = all_array['TimesBMJD_TDB']
rv_array = all_array['RV_kms']
sigma_array = all_array['Sigma_kms']+0.5

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

photo_zero_times = photometry_times - zero_point
#photo_period = period/2.
photo_folded_times = np.mod(photo_zero_times, period)/period
#fitted_photo_curve, fitted_photo_cov = sciop.curve_fit(photo_sine_function, photo_folded_times, photometry_flux, sigma= photometry_error, p0= [0, 0.08])
#fitted_photo_curve, fitted_photo_cov = sciop.curve_fit(photo_sine_function, photo_folded_times, photometry_flux, sigma= photometry_error, bounds = photo_bounds)
#photo_residuals= photometry_flux - photo_sine_function(photo_folded_times, fitted_photo_curve[0], fitted_photo_curve[1])
fitted_photo_curve, fitted_photo_cov = sciop.curve_fit(photo_harmonic_function, photo_folded_times, photometry_flux, sigma=photometry_error)
print fitted_photo_curve
#photo_residuals= photometry_flux - photo_sine_function(photo_folded_times, fitted_photo_curve[0], fitted_photo_curve[1], fitted_photo_curve[2])
#photo_residuals = photometry_flux - photo_harmonic_function(photo_folded_times, fitted_photo_curve[0], fitted_photo_curve[1], fitted_photo_curve[2], fitted_photo_curve[3], fitted_photo_curve[4], fitted_photo_curve[5], fitted_photo_curve[6], fitted_photo_curve[7], fitted_photo_curve[8])
#photo_residuals = photometry_flux - photo_harmonic_function(photo_folded_times, fitted_photo_curve[0], fitted_photo_curve[1], fitted_photo_curve[2], fitted_photo_curve[3], fitted_photo_curve[4])
photo_residuals = photometry_flux - photo_harmonic_function(photo_folded_times, fitted_photo_curve[0], fitted_photo_curve[1], fitted_photo_curve[2])


coeffs = fitted_photo_curve

#fitted_curve_all, fitted_cov_all = sciop.curve_fit(sine_function, folded_times, rv_array, p0= p0_list)
fitted_curve_all, fitted_cov_all = sciop.curve_fit(sine_function, folded_times, rv_array, sigma = sigma_array, p0= p0_list)
residuals = rv_array- sine_function(folded_times, fitted_curve_all[0], fitted_curve_all[1])
print fitted_curve_all

####
model_rvs = sine_function(folded_times, fitted_curve_all[0], fitted_curve_all[1])
output_array = np.vstack([bmjd_array, model_rvs]).T
print "Saving model RV's"
np.savetxt(output_filename, output_array, delimiter= '\t', header= 'TimesBMJD_TDB\tRV_kms')

#####

def get_mass_ratio(K_c, P_B, x_psr):
    """
    Mass ratio M_psr/M_companion
    """
    print (K_c.to(u.m/u.s)*P_B.to(u.s)/(2*np.pi*x_psr.to(u.s)))
    return K_c.to(u.m/u.s)*P_B.to(u.s)/(2*np.pi*x_psr.to(u.s)*c)

period_wunits = period*u.day

all_mratio= get_mass_ratio(fitted_curve_all[1]*(u.km/u.second),period_wunits, x_psr)

x_vals = np.linspace(0,1,1000)
#sine_vals = sine_function(x_vals, fitted_curve_all[0], fitted_curve_all[1], fitted_curve_all[2])
sine_vals = sine_function(x_vals, fitted_curve_all[0], fitted_curve_all[1])
#photo_sine_vals = photo_sine_function(x_vals, fitted_photo_curve[0], fitted_photo_curve[1])
#photo_sine_vals = photo_sine_function(x_vals, fitted_photo_curve[0], fitted_photo_curve[1], fitted_photo_curve[2])
#photo_harmonic_vals = photo_harmonic_function(x_vals, fitted_photo_curve[0], fitted_photo_curve[1], fitted_photo_curve[2], fitted_photo_curve[3], fitted_photo_curve[4], fitted_photo_curve[5], fitted_photo_curve[6], fitted_photo_curve[7], fitted_photo_curve[8])
#photo_harmonic_vals = photo_harmonic_function(x_vals, fitted_photo_curve[0], fitted_photo_curve[1], fitted_photo_curve[2], fitted_photo_curve[3], fitted_photo_curve[4])
photo_harmonic_vals = photo_harmonic_function(x_vals, fitted_photo_curve[0], fitted_photo_curve[1], fitted_photo_curve[2])



harmonic_strings = make_value_strings(fitted_photo_curve, fitted_photo_cov)


#### RVs 
#fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True)

fig = plt.figure()
#ax = fig.add_subplot(311)
ax = plt.subplot2grid((6, 1), (0,0), rowspan = 2)
#ax.axhline(y= fitted_curve_all[0], color = 'r', linestyle = ':', alpha = 0.4, label= str(np.round(fitted_curve_all[0],precision))+' km/s')
ax.axhline(y= fitted_curve_all[0], color = 'r', linestyle = ':', alpha = 1, label= str(np.round(fitted_curve_all[0],precision))+' km/s')

#ax.scatter(folded_times, rv_array, color = 'b', label = '3/18/18')
#ax.plot(folded_times, rv_array, color = 'b', marker = 'o', linestyle = 'none', label = '3/18/18')
ax.errorbar(folded_times, rv_array, yerr= sigma_array, color = 'b', marker = 'o', linestyle = 'none', label = '3/18/18')
#ax.scatter(folded_timesB, rv_arrayB, color = 'r', label = '2/13/18')
ax.plot(x_vals, sine_vals, color = 'k', linestyle = '--', label = 'Model')
ax.set_xticklabels([])
ax.set_xlim([0,1])
#ax.set_xlabel('Phase')
ax.set_ylabel('RV (km/s)')
ax.set_title(str(np.round(1.4/all_mratio,precision2)) + ' M_sun companion assuming 1.4Msun NS')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
text_string =  r'$K_c =$ ' +str(np.round(fitted_curve_all[1], precision)) + r'$\pm$ ' + str(np.round(np.sqrt(fitted_cov_all[1,1]),precision)) + 'km/s' +'\n'+ r'$v_{sys}=$' +str(np.round(fitted_curve_all[0],precision)) + r'$\pm$' + str(np.round(np.sqrt(fitted_cov_all[0,0]),precision))+ ' km/s'
#ax.text(0.2, 0.05,text_string, transform=ax.transAxes, fontsize=14, verticalalignment='bottom', bbox=props)
ax.text(0.2, 0.05,text_string, transform=ax.transAxes, verticalalignment='bottom', bbox=props)

ax.legend()

#plt.show()

#### Residuals

#fig = plt.figure()
#ax2 = fig.add_subplot(613)
ax2 =plt.subplot2grid((6,1), (2, 0), rowspan= 1)

#ax2.scatter(folded_times, residuals, color = 'b', label = '3/18/18')
#ax2.plot(folded_times, residuals, color = 'b', label = '3/18/18', marker = 'o', linestyle = 'none')
ax2.errorbar(folded_times, residuals, yerr= sigma_array, color = 'b', label = '3/18/18', marker = 'o', linestyle = 'none')

#ax.scatter(folded_timesB, rv_arrayB, color = 'r', label = '2/13/18')
#ax.plot(x_vals, sine_vals, color = 'k', linestyle = '--', label = 'Model')
#ax2.plot(x_vals, np.zeros(x_vals.shape), color = 'k', linestyle = '--', label = 'Model')
ax2.axhline(y=0, color = 'k', linestyle = '--')
ax2.set_xticklabels([])
ax2.set_xlim([0,1])
#ax2.set_xlabel('Phase')
#ax2.set_ylabel('Residuals (km/s)')
ax2.set_ylabel('Residuals')

#ax2.set_title(str(np.round(1.4/all_mratio,precision)) + ' M_sun companion assuming 1.4Msun NS')
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#text_string =  r'$K_c =$ ' +str(np.round(fitted_curve_all[1], precision)) + r'$\pm$ ' + str(np.round(np.sqrt(fitted_cov_all[1,1]),precision)) + 'km/s' +'\n'+ r'$v_{sys}=$' +str(np.round(fitted_curve_all[0],precision)) + r'$\pm$' + str(np.round(np.sqrt(fitted_cov_all[0,0]),precision))+ ' km/s'
#ax.text(0.05, 0.05,text_string, transform=ax.transAxes, fontsize=14, verticalalignment='bottom', bbox=props)
#ax.legend()

#plt.show()

###3 Photometry

#fig = plt.figure()
#ax3 = fig.add_subplot(313)
ax3 = plt.subplot2grid((6,1), (3, 0), rowspan=2)
#ax3.axhline(y= fitted_photo_curve[0], color = 'r', linestyle = ':', alpha = 0.4, label= str(np.round(fitted_photo_curve[0],precision2)))
ax3.axhline(y= fitted_photo_curve[0], color = 'r', linestyle = ':', alpha = 1, label= str(np.round(fitted_photo_curve[0],precision2)))

ax3.errorbar(photo_folded_times, photometry_flux, photometry_error, color = 'b', label = '4/22/18', linestyle= 'None', marker = 'o')
#ax.scatter(folded_timesB, rv_arrayB, color = 'r', label = '2/13/18')
#ax3.plot(x_vals, photo_sine_vals, color = 'k', linestyle = '--', label = 'Model')
ax3.plot(x_vals, photo_harmonic_vals, color = 'k', linestyle = '--', label = 'Model')
ax3.set_xticklabels([])
ax3.set_xlim([0,1])
#ax3.set_xlabel('Phase')
ax3.set_ylabel('Normalized Flux')
#ax3.set_title(str(np.round(1.4/all_mratio,precision)) + ' M_sun companion assuming 1.4Msun NS')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#text_string =  r'$A =$ ' +str(np.round(fitted_photo_curve[1], precision)) + r'$\pm$ ' + str(np.round(np.sqrt(fitted_photo_cov[1,1]),precision)) + '' +'\n'+ r'$b=$' +str(np.round(fitted_photo_curve[0],precision)) + r'$\pm$' + str(np.round(np.sqrt(fitted_photo_cov[0,0]),precision))+ ''
#text_string =  r'$A =$ ' +str(np.round(fitted_photo_curve[1], precision)) + r'$\pm$ ' + str(np.round(np.sqrt(fitted_photo_cov[1,1]),precision)) + '' +'\n'+ r'$b=$' +str(np.round(fitted_photo_curve[0],precision)) + r'$\pm$' + str(np.round(np.sqrt(fitted_photo_cov[0,0]),precision))+ '' + '\n' + 'phase = ' + str(np.round(fitted_photo_curve[2], precision)) + r'$\pm$' + str(np.round(np.sqrt(fitted_photo_cov[2,2]), precision))
#text_string =  harmonic_strings[1] + r'$*\sin(2*\pi*t + $' +harmonic_strings[5]+ r'$) +$' + harmonic_strings[2] + r'$*\cos(2*\pi*t + $' +harmonic_strings[6] +  r'$) +$' + harmonic_strings[3]+ r'$*\sin(4*\pi*t + $' +harmonic_strings[7]+  r'$) +$ ' + harmonic_strings[4] + r'$*\cos(4*\pi*t + $'+ harmonic_strings[8]+ r'$) + $' + harmonic_strings[0]
#text_string =  harmonic_strings[1] + r'$*\sin(2*\pi*t/P) +$' + harmonic_strings[2] + r'$*\cos(2*\pi*t/P) +$' + harmonic_strings[3]+ r'$*\sin(4*\pi*t/P) +$ ' + harmonic_strings[4] + r'$*\cos(4*\pi*t/P) + $' + harmonic_strings[0]
#text_string =   harmonic_strings[1] + r'$*\cos(2*\pi*t/P) +$' + harmonic_strings[2] + r'$*\cos(4*\pi*t/P) + $' + harmonic_strings[0]
text_string =   "F= " + harmonic_strings[1] + r'$*\cos(2*\pi*t/P) +$' + harmonic_strings[2] + r'$*\cos(4*\pi*t/P) + $' + harmonic_strings[0]

#ax3.text(0.2, 0.05,text_string, transform=ax3.transAxes, fontsize=14, verticalalignment='bottom', bbox=props)
#ax3.text(0.05, 0.05,text_string, transform=ax3.transAxes, fontsize=14, verticalalignment='bottom', bbox=props)
ax3.text(0.05, 0.05,text_string, transform=ax3.transAxes, verticalalignment='bottom', bbox=props)

ax3.legend()

#plt.subplots_adjust(wspace = 0, hspace = 0, top = 0.95, bottom = 0.05, left = 0.05, right = 0.95)
#plt.show()

#### Photometric Residuals

#fig = plt.figure()
#ax2 = fig.add_subplot(613)
ax4 =plt.subplot2grid((6,1), (5, 0), rowspan= 1)

ax4.errorbar(photo_folded_times, photo_residuals, photometry_error, color = 'b', label = '4/22/18', linestyle= 'None', marker = 'o')
#ax4.errorbar(photo_folded_times, photo_residuals, photometry_error, color = 'b', label = '4/22/18')
#ax.scatter(folded_timesB, rv_arrayB, color = 'r', label = '2/13/18')
#ax.plot(x_vals, sine_vals, color = 'k', linestyle = '--', label = 'Model')
#ax2.plot(x_vals, np.zeros(x_vals.shape), color = 'k', linestyle = '--', label = 'Model')
ax4.axhline(y=0, color = 'k', linestyle = '--')
ax4.set_xlabel('Phase')
#ax4.set_xticklabels([])
ax4.set_xlim([0,1])
#ax2.set_xlabel('Phase')
ax4.set_ylabel('Residuals')
#ax2.set_title(str(np.round(1.4/all_mratio,precision)) + ' M_sun companion assuming 1.4Msun NS')
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#text_string =  r'$K_c =$ ' +str(np.round(fitted_curve_all[1], precision)) + r'$\pm$ ' + str(np.round(np.sqrt(fitted_cov_all[1,1]),precision)) + 'km/s' +'\n'+ r'$v_{sys}=$' +str(np.round(fitted_curve_all[0],precision)) + r'$\pm$' + str(np.round(np.sqrt(fitted_cov_all[0,0]),precision))+ ' km/s'
#ax.text(0.05, 0.05,text_string, transform=ax.transAxes, fontsize=14, verticalalignment='bottom', bbox=props)
#ax.legend()
plt.subplots_adjust(wspace = 0, hspace = 0, top = 0.93, bottom = 0.07, left = 0.07, right = 0.93)

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
