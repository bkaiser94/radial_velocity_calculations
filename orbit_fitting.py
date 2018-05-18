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
#plt.rc('font', size =18)
#plt.rc('lines', markersize=12)
plt.rc('font', size = 11)
plt.rc('lines', markersize = 5)
plotting_offset = 0.0005



plotfile = glob(sys.argv[3])
print plotfile
c= 2.998E8 *u.m/u.s
ra = float(sys.argv[1]) #values in decimal degrees
dec = float(sys.argv[2])

target_coord = coord.SkyCoord(ra, dec, unit= (u.deg, u.deg), frame= 'icrs')

def to_barycenter(input_times):
    bary_corr =input_times.tdb.light_travel_time(target_coord)
    return (input_times.tdb+ bary_corr.tdb).mjd


all_array = np.genfromtxt(plotfile[0], names = True, delimiter = ',')

parkes_location = coord.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
cerro_pachon_location = coord.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)

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

#epoch_difference = obs_times[0].mjd-tasc.mjd

bmjd_array = all_array['BMJD_TDB']
H_delta = all_array['H_delta']
H_gamma = all_array["H_gamma"]
H_beta = all_array['H_beta']
H_delta_s = all_array["H_delta_s"]
H_gamma_s = all_array['H_gamma_s']
H_beta_s = all_array['H_beta_s']
#print Time(mjd_array, format = 'mjd').utc.isot
#print all_array

#p0_list = [-100, period, 300,bmjd_array[0]-0.1]
#p0_list = [-100, 300,bmjd_array[0]-0.1]
p0_list = [-100, 300]


merged_lines = np.hstack([H_beta, H_gamma, H_delta])
merged_bmjd = np.hstack([bmjd_array, bmjd_array, bmjd_array])
merged_errors = np.hstack([H_beta_s, H_gamma_s, H_delta_s])

#bmjd_obs = to_barycenter(obs_times) #corrected to barycenter to use against the rv curve
bmjd_t0 = to_barycenter(t0) #corrected initial epoch
bmjd_tasc= to_barycenter(tasc)
bmjd_tconj = to_barycenter(tconj)

zero_point = bmjd_tconj

#def sine_function(times, systemic_vel,  amplitude, zero_point, period):
def sine_function(times, systemic_vel,  amplitude, zero_point):
#def sine_function(times, systemic_vel,  amplitude):
    return amplitude*np.sin((2*np.pi)/period*(times-zero_point))+systemic_vel


def calc_chi_square(rvs, sigmas, sine_function_vals):
    return np.sum((rvs-sine_function_vals)**2/(sigmas**2))

def calc_red_chi_square(rvs, sigmas, sine_function_vals, dof = 0):
    print "degrees of freedom", rvs.shape[0]-dof
    return calc_chi_square(rvs, sigmas, sine_function_vals)/np.float_(rvs.shape[0]- dof)



#fitted_curve_delta, fitted_cov_delta = sciop.curve_fit(sine_function, bmjd_array, H_delta, sigma=H_delta_s, p0=p0_list)
#fitted_curve_gamma, fitted_cov_gamma = sciop.curve_fit(sine_function, bmjd_array, H_gamma, sigma= H_gamma_s, p0=p0_list)
#fitted_curve_beta, fitted_cov_beta = sciop.curve_fit(sine_function, bmjd_array, H_beta, sigma = H_beta_s, p0=p0_list)
#fitted_curve_all, fitted_cov_all = sciop.curve_fit(sine_function, merged_bmjd, merged_lines, p0= p0_list)
#fitted_curve_all, fitted_cov_all = sciop.curve_fit(sine_function, merged_bmjd, merged_lines, p0= p0_list+[bmjd_array[0]])
fitted_curve_all, fitted_cov_all = sciop.curve_fit(sine_function, merged_bmjd, merged_lines, p0= p0_list+[bmjd_array[0]])
#fitted_curve_all, fitted_cov_all = sciop.curve_fit(sine_function, merged_bmjd, merged_lines, p0= p0_list+[bmjd_array[0],period])

print fitted_curve_all
print "difference in ephemeris (tconj-fit) by phase"
print  ((zero_point - fitted_curve_all[2])%period)/period

#delta_fit_errs = np.sqrt(np.diag(fitted_cov_delta))
#gamma_fit_errs = np.sqrt(np.diag(fitted_cov_gamma))
#beta_fit_errs = np.sqrt(np.diag(fitted_cov_beta))
all_fit_errs= np.sqrt(np.diag(fitted_cov_all))
#print "degrees of freedom: ", H_beta.shape[0]-2

all_lines_bmjd = np.hstack([bmjd_array, bmjd_array, bmjd_array]).T
all_lines_bmjd = all_lines_bmjd.ravel()


#beta_chi_square = scistats.chisquare(H_beta, sine_function(bmjd_array,fitted_curve_beta[0],fitted_curve_beta[1]),ddof=-2)
#gamma_chi_square = scistats.chisquare(H_gamma, sine_function(bmjd_array,fitted_curve_gamma[0],fitted_curve_gamma[1]),ddof=-2)
#delta_chi_square = scistats.chisquare(H_delta, sine_function(bmjd_array,fitted_curve_delta[0],fitted_curve_delta[1]),ddof=-2)
#all_chi_square = scistats.chisquare(merged_lines.ravel(), sine_function(all_lines_bmjd, fitted_curve_all[0],fitted_curve_all[1]) , ddof = -2)


#print "chi-square values"
#print "beta", "gamma", "delta"
#print beta_chi_square[0], gamma_chi_square[0], delta_chi_square[0]
#print all_chi_square

#all_chi_square= calc_chi_square(merged_lines.ravel(), merged_errors.ravel(), sine_function(all_lines_bmjd, fitted_curve_all[0], fitted_curve_all[1]))
all_chi_square= calc_chi_square(merged_lines.ravel(), merged_errors.ravel(), sine_function(all_lines_bmjd, fitted_curve_all[0], fitted_curve_all[1], fitted_curve_all[2]))
#all_chi_square= calc_chi_square(merged_lines.ravel(), merged_errors.ravel(), sine_function(all_lines_bmjd, fitted_curve_all[0], fitted_curve_all[1], fitted_curve_all[2],fitted_curve_all[3]))


print all_chi_square

#red_chi =  calc_red_chi_square(merged_lines.ravel(), merged_errors.ravel(), sine_function(all_lines_bmjd, fitted_curve_all[0], fitted_curve_all[1]), dof=2)
red_chi =  calc_red_chi_square(merged_lines.ravel(), merged_errors.ravel(), sine_function(all_lines_bmjd, fitted_curve_all[0], fitted_curve_all[1], fitted_curve_all[2]), dof=3)
#red_chi =  calc_red_chi_square(merged_lines.ravel(), merged_errors.ravel(), sine_function(all_lines_bmjd, fitted_curve_all[0], fitted_curve_all[1], fitted_curve_all[2], fitted_curve_all[3]), dof=4)


print red_chi
#print calc_red_chi_square(H_beta,H_beta_s, sine_function(bmjd_array, fitted_curve_beta[0], fitted_curve_beta[1]), dof=2)
print '-------'


#print "sys vel.   period    amplitude    zero_point"
#print fitted_curve_delta
#print fitted_curve_gamma
#print fitted_curve_beta

def zero_rvs(rv_array):
    print "systemic velocity (includes Earth's motion):", np.mean([rv_array.max(),rv_array.min()])
    return rv_array-np.mean([rv_array.max(),rv_array.min()])

def get_mass_ratio(K_c, P_B, x_psr):
    """
    Mass ratio M_psr/M_companion
    """
    print (K_c.to(u.m/u.s)*P_B.to(u.s)/(2*np.pi*x_psr.to(u.s)))
    return K_c.to(u.m/u.s)*P_B.to(u.s)/(2*np.pi*x_psr.to(u.s)*c)

period_wunits = period*u.day

#H_delta_mratio= get_mass_ratio(fitted_curve_delta[1]*(u.km/u.second),period_wunits, x_psr)
#H_gamma_mratio= get_mass_ratio(fitted_curve_gamma[1]*(u.km/u.second),period_wunits, x_psr)
#H_beta_mratio= get_mass_ratio(fitted_curve_beta[1]*(u.km/u.second),period_wunits, x_psr)
all_mratio= get_mass_ratio(fitted_curve_all[1]*(u.km/u.second),period_wunits, x_psr)

#H_delta_mratio_high= get_mass_ratio((fitted_curve_delta[1]+delta_fit_errs[1])*(u.km/u.second),period_wunits, x_psr)
#H_gamma_mratio_high= get_mass_ratio((fitted_curve_gamma[1]+gamma_fit_errs[1])*(u.km/u.second),period_wunits, x_psr)
#H_beta_mratio_high= get_mass_ratio((fitted_curve_beta[1]+beta_fit_errs[1])*(u.km/u.second),period_wunits, x_psr)
all_mratio_high =get_mass_ratio((fitted_curve_all[1]+ all_fit_errs[1])*(u.km/u.second),period_wunits, x_psr)



#H_delta_mratio_low= get_mass_ratio((fitted_curve_delta[1]-delta_fit_errs[1])*(u.km/u.second),period_wunits, x_psr)
#H_gamma_mratio_low= get_mass_ratio((fitted_curve_gamma[1]-gamma_fit_errs[1])*(u.km/u.second),period_wunits, x_psr)
#H_beta_mratio_low= get_mass_ratio((fitted_curve_beta[1]-beta_fit_errs[1])*(u.km/u.second),period_wunits, x_psr)
all_mratio_low =get_mass_ratio((fitted_curve_all[1]-all_fit_errs[1])*(u.km/u.second),period_wunits, x_psr)
print "Companion masses assuming 1.4 Msun pulsar"
#print "delta", 1.4/H_delta_mratio
#print "gamma", 1.4/H_gamma_mratio
#print "beta", 1.4/H_beta_mratio
print "all", 1.4/all_mratio
print "Companion masses uncertainty assuming 1.4 Msun pulsar"
#print "delta", 1.4/H_delta_mratio_high-1.4/H_delta_mratio
#print "gamma", 1.4/H_gamma_mratio_high-1.4/H_gamma_mratio
#print "beta", 1.4/H_beta_mratio_high- 1.4/H_beta_mratio
print "all", 1.4/all_mratio_high-1.4/all_mratio
print "Companion masses uncertainty assuming 1.4 Msun pulsar"
#print "delta", 1.4/H_delta_mratio_low-1.4/H_delta_mratio
#print "gamma", 1.4/H_gamma_mratio_low-1.4/H_gamma_mratio
#print "beta", 1.4/H_beta_mratio_low-1.4/H_beta_mratio
print "all", 1.4/all_mratio_low-1.4/all_mratio
#H_delta = zero_rvs(H_delta)
#H_gamma = zero_rvs(H_gamma)
#H_beta = zero_rvs(H_beta)
#remean_rv = np.mean([H_delta, H_gamma, H_beta], axis = 0)
#remean_std = np.std([H_delta,H_gamma, H_beta], axis=0)
#mean_rv = zero_rvs(mean_rv)
#plt.plot(bmjd_array, H_delta, label = r"H-$\delta$", linestyle = 'none', marker = '*')
#plt.plot(bmjd_array, H_gamma, label = r"H-$\gamma$", linestyle = 'none', marker = '*')
#plt.plot(bmjd_array, H_beta, label = r"H-$\beta$", linestyle = 'none', marker = '*')
plt.errorbar(bmjd_array, H_delta, H_delta_s, label = r"H-$\delta$", linestyle = 'none', marker = '*', color = 'b')
plt.errorbar(bmjd_array+plotting_offset, H_gamma, H_gamma_s, label = r"H-$\gamma$", linestyle = 'none', marker = '*', color = 'g')
plt.errorbar(bmjd_array+plotting_offset*2, H_beta, H_beta_s, label = r"H-$\beta$", linestyle = 'none', marker = '*', color = 'r')

plot_times = np.linspace(bmjd_array[0]-0.2, bmjd_array[-1]+0.2, 1000)
#plt.plot(plot_times, sine_function(plot_times, fitted_curve_delta[0],fitted_curve_delta[1],fitted_curve_delta[2], fitted_curve_delta[3]), label = 'fitted curve delta')
#plt.plot(plot_times, sine_function(plot_times, fitted_curve_gamma[0],fitted_curve_gamma[1],fitted_curve_gamma[2], fitted_curve_gamma[3]), label = 'fitted curve gamma')
#plt.plot(plot_times, sine_function(plot_times, fitted_curve_beta[0],fitted_curve_beta[1],fitted_curve_beta[2], fitted_curve_beta[3]), label = 'fitted curve beta')
#plt.plot(plot_times, sine_function(plot_times, fitted_curve_delta[0],fitted_curve_delta[1],fitted_curve_delta[2]), label = 'fitted curve delta')
#plt.plot(plot_times, sine_function(plot_times, fitted_curve_gamma[0],fitted_curve_gamma[1],fitted_curve_gamma[2]), label = 'fitted curve gamma')
#plt.plot(plot_times, sine_function(plot_times, fitted_curve_beta[0],fitted_curve_beta[1],fitted_curve_beta[2]), label = 'fitted curve beta')

#delta_label = 'fit delta K:'+str(np.round(fitted_curve_delta[1],2))+'(' + str(np.round(delta_fit_errs[1],2))+')'+'(km/s) sys:' +str(np.round(fitted_curve_delta[0],2)) +'(' + str(np.round(delta_fit_errs[0],2))+')'+'(km/s)'
#gamma_label='fit gamma K:'+str(np.round(fitted_curve_gamma[1],2))+'(' + str(np.round(gamma_fit_errs[1],2))+')'+'(km/s) sys:' +str(np.round(fitted_curve_gamma[0],2)) +'(' + str(np.round(gamma_fit_errs[0],2))+')'+'(km/s)'
#beta_label= 'fit beta K:'+str(np.round(fitted_curve_beta[1],2))+'(' + str(np.round(beta_fit_errs[1],2))+')'+'(km/s) sys:' +str(np.round(fitted_curve_beta[0],2))+'(' + str(np.round(beta_fit_errs[0],2))+')'+'(km/s)'
all_label = 'fit all K:'+str(np.round(fitted_curve_all[1],2))+'(' + str(np.round(all_fit_errs[1],2))+')'+'(km/s) sys:' +str(np.round(fitted_curve_all[0],2))+'(' + str(np.round(all_fit_errs[0],2))+')'+'(km/s)'

#plt.plot(plot_times, sine_function(plot_times, fitted_curve_delta[0],fitted_curve_delta[1]), label = delta_label, color = 'b')
#plt.plot(plot_times, sine_function(plot_times, fitted_curve_gamma[0],fitted_curve_gamma[1]), label = gamma_label, color = 'g')
#plt.plot(plot_times, sine_function(plot_times, fitted_curve_beta[0],fitted_curve_beta[1]), label = beta_label, color = 'r')
#plt.plot(plot_times, sine_function(plot_times, fitted_curve_all[0],fitted_curve_all[1]), label = all_label , color = 'k')
plt.plot(plot_times, sine_function(plot_times, fitted_curve_all[0],fitted_curve_all[1], fitted_curve_all[2]), label = all_label , color = 'k')
#plt.plot(plot_times, sine_function(plot_times, fitted_curve_all[0],fitted_curve_all[1], fitted_curve_all[2], fitted_curve_all[3]), label = all_label , color = 'k')


#plt.legend()
#plt.title("Period: " + str(fitted_curve_beta[1] )+ " Sys Vel.: " + str(fitted_curve_beta[0]) + " K: " + str(fitted_curve_beta[2]))
plt.ylabel('RV (km/s)')
plt.xlabel('BMJD_TDB')
#plt.title("Period: " + str(period)+ " Sys Vel.: " + str(fitted_curve_beta[0]) + " K: " + str(fitted_curve_beta[1]))
#plt.title("Period: " +str(period) + "days, Phase = 0 at Conjunction: " + str(bmjd_tconj))
plt.title('PSR J1431-4715')
#plt.errorbar(mjd_array, mean_rv, std_dev, label = 'Mean RV', linestyle = 'none', marker = 'o')
#plt.errorbar(mjd_array, remean_rv, remean_std, label = r"Mean of zeroed RV's", linestyle = 'none', marker = 'o')
plt.show()


f, (ax1, ax2) = plt.subplots(2,1, sharex=True, sharey = False)
ax1.errorbar(bmjd_array, H_delta, H_delta_s, label = r"H-$\delta$", linestyle = 'none', marker = '*', color = 'b')
ax1.errorbar(bmjd_array+plotting_offset, H_gamma, H_gamma_s, label = r"H-$\gamma$", linestyle = 'none', marker = 'o', color = 'g')
ax1.errorbar(bmjd_array+plotting_offset*2, H_beta, H_beta_s, label = r"H-$\beta$", linestyle = 'none', marker = 's', color = 'r')

plot_times = np.linspace(bmjd_array[0]-0.2, bmjd_array[-1]+0.2, 1000)
#delta_label = 'fit delta K:'+str(np.round(fitted_curve_delta[1],2))+'(' + str(np.round(delta_fit_errs[1],2))+')'+'(km/s) sys:' +str(np.round(fitted_curve_delta[0],2)) +'(' + str(np.round(delta_fit_errs[0],2))+')'+'(km/s)'
#gamma_label='fit gamma K:'+str(np.round(fitted_curve_gamma[1],2))+'(' + str(np.round(gamma_fit_errs[1],2))+')'+'(km/s) sys:' +str(np.round(fitted_curve_gamma[0],2)) +'(' + str(np.round(gamma_fit_errs[0],2))+')'+'(km/s)'
#beta_label= 'fit beta K:'+str(np.round(fitted_curve_beta[1],2))+'(' + str(np.round(beta_fit_errs[1],2))+')'+'(km/s) sys:' +str(np.round(fitted_curve_beta[0],2))+'(' + str(np.round(beta_fit_errs[0],2))+')'+'(km/s)'
all_label = 'fit all K:'+str(np.round(fitted_curve_all[1],2))+'(' + str(np.round(all_fit_errs[1],2))+')'+'(km/s) sys:' +str(np.round(fitted_curve_all[0],2))+'(' + str(np.round(all_fit_errs[0],2))+')'+'(km/s)'

#ax1.plot(plot_times, sine_function(plot_times, fitted_curve_delta[0],fitted_curve_delta[1]), label = delta_label, color = 'b')
#ax1.plot(plot_times, sine_function(plot_times, fitted_curve_gamma[0],fitted_curve_gamma[1]), label = gamma_label, color = 'g')

#ax1.plot(plot_times, sine_function(plot_times, fitted_curve_beta[0],fitted_curve_beta[1]), label = beta_label, color = 'r')

#ax1.plot(plot_times, sine_function(plot_times, fitted_curve_all[0],fitted_curve_all[1]), label = all_label , color = 'k')
#ax1.plot(plot_times, sine_function(plot_times, fitted_curve_all[0],fitted_curve_all[1], fitted_curve_all[2]), label = all_label , color = 'k')
ax1.plot(plot_times, sine_function(plot_times, fitted_curve_all[0],fitted_curve_all[1], fitted_curve_all[2], fitted_curve_all[3]), label = all_label , color = 'k')


ax1.set_ylabel('RV (km/s)')
#ax1.set_xlabel('BMJD_TDB')


ax2.axhline(y=0, color = 'k')
#fitted_points = sine_function(bmjd_array, fitted_curve_all[0], fitted_curve_all[1])
#fitted_points = sine_function(bmjd_array, fitted_curve_all[0], fitted_curve_all[1], fitted_curve_all[2])
fitted_points = sine_function(bmjd_array, fitted_curve_all[0], fitted_curve_all[1], fitted_curve_all[2], fitted_curve_all[3])


residuals_delta = H_delta-fitted_points
residuals_gamma = H_gamma-fitted_points
residuals_beta = H_beta- fitted_points

ax2.errorbar(bmjd_array, residuals_delta, H_delta_s, marker = '*', color= 'b', linestyle = 'none')
ax2.errorbar(bmjd_array, residuals_gamma, H_gamma_s, marker = 'o', color = 'g', linestyle= 'none')
ax2.errorbar(bmjd_array, residuals_beta, H_beta_s, marker = 's', color= 'r', linestyle = 'none')
ax2.set_xlabel('BMJD_TDB')
ax2.set_ylabel('Residuals (km/s)')
f.subplots_adjust(wspace=0)
f.subplots_adjust(hspace = 0)
plt.show()

########### attempt with the different subplot sizes
#f, (ax1, ax2) = plt.subplots(2,1, sharex=True, sharey = False)
range_x = (np.min(bmjd_array)-0.05, np.max(bmjd_array)+0.05)
plt.subplot(413)
plt.axhline(y=0, color = 'k', linestyle = '--')
#fitted_points = sine_function(bmjd_array, fitted_curve_all[0], fitted_curve_all[1])
fitted_points = sine_function(bmjd_array, fitted_curve_all[0], fitted_curve_all[1] , fitted_curve_all[2])

residuals_delta = H_delta-fitted_points
residuals_gamma = H_gamma-fitted_points
residuals_beta = H_beta- fitted_points

plt.errorbar(bmjd_array, residuals_delta, H_delta_s, marker = '*', color= 'b', linestyle = 'none')
plt.errorbar(bmjd_array, residuals_gamma, H_gamma_s, marker = 'o', color = 'g', linestyle= 'none')
plt.errorbar(bmjd_array, residuals_beta, H_beta_s, marker = 's', color= 'r', linestyle = 'none')
plt.xlim(range_x)
plt.xlabel('BMJD_TDB')
plt.ylabel('Residuals')


plt.subplot(211)
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
plt.errorbar(bmjd_array, H_delta, H_delta_s, label = r"H-$\delta$", linestyle = 'none', marker = '*', color = 'b')
plt.errorbar(bmjd_array+plotting_offset, H_gamma, H_gamma_s, label = r"H-$\gamma$", linestyle = 'none', marker = 'o', color = 'g')
plt.errorbar(bmjd_array+plotting_offset*2, H_beta, H_beta_s, label = r"H-$\beta$", linestyle = 'none', marker = 's', color = 'r')

plot_times = np.linspace(bmjd_array[0]-0.2, bmjd_array[-1]+0.2, 1000)
delta_label = 'fit delta K:'+str(np.round(fitted_curve_delta[1],2))+'(' + str(np.round(delta_fit_errs[1],2))+')'+'(km/s) sys:' +str(np.round(fitted_curve_delta[0],2)) +'(' + str(np.round(delta_fit_errs[0],2))+')'+'(km/s)'
gamma_label='fit gamma K:'+str(np.round(fitted_curve_gamma[1],2))+'(' + str(np.round(gamma_fit_errs[1],2))+')'+'(km/s) sys:' +str(np.round(fitted_curve_gamma[0],2)) +'(' + str(np.round(gamma_fit_errs[0],2))+')'+'(km/s)'
beta_label= 'fit beta K:'+str(np.round(fitted_curve_beta[1],2))+'(' + str(np.round(beta_fit_errs[1],2))+')'+'(km/s) sys:' +str(np.round(fitted_curve_beta[0],2))+'(' + str(np.round(beta_fit_errs[0],2))+')'+'(km/s)'
all_label = 'fit all K:'+str(np.round(fitted_curve_all[1],2))+'(' + str(np.round(all_fit_errs[1],2))+')'+'(km/s) sys:' +str(np.round(fitted_curve_all[0],2))+'(' + str(np.round(all_fit_errs[0],2))+')'+'(km/s)'

#ax1.plot(plot_times, sine_function(plot_times, fitted_curve_delta[0],fitted_curve_delta[1]), label = delta_label, color = 'b')
#ax1.plot(plot_times, sine_function(plot_times, fitted_curve_gamma[0],fitted_curve_gamma[1]), label = gamma_label, color = 'g')

#ax1.plot(plot_times, sine_function(plot_times, fitted_curve_beta[0],fitted_curve_beta[1]), label = beta_label, color = 'r')

plt.plot(plot_times, sine_function(plot_times, fitted_curve_all[0],fitted_curve_all[1]), label = all_label , color = 'k', linestyle = '--')
plt.ylabel('RV (km/s)')
plt.xlim(range_x)
#ax1.set_xlabel('BMJD_TDB')

plt.subplots_adjust(wspace=0)
plt.subplots_adjust(hspace = 0)
#plt.tight_layout()
#plt.title('PSR J1431-4715')
plt.show()


#plt.subplot(211)
#plt.plot([1,2,3])

#plt.subplot(413)
#plt.plot([1,2,3])
#plt.subplots_adjust(wspace=0)
#plt.show()



