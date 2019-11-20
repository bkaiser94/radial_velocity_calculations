"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-10-04

Plot the WD cooling models with Teff and logg values hopefully in an attempt to figure out which one best fits the 
target.



"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from astropy.table import Table, Column
import scipy.interpolate as scinterp


import spec_plot_tools as spt
import cal_params as cp


cooling_model_file='COModel_ThinH.csv'


cooling_model_file=cp.ref_dir+'WD_cooling_models/'+cooling_model_file
#target_logg= 7.85
#target_teff= 4030. #K

#target_logg=7.98
#target_teff= 4040. #K

n=100000
target_logg=7.77
target_logg_err= 0.23
target_teff= 3830.
target_teff_err= 230.

simon_mass= 0.45
simon_mass_err= 0.12

target_logg_dist= np.copy(np.random.normal(loc=target_logg, scale=target_logg_err, size=n))
target_teff_dist= np.copy(np.random.normal(loc=target_teff, scale=target_teff_err, size=n))
simon_mass_dist= np.copy(np.random.normal(loc=simon_mass, scale=simon_mass_err, size=n))


#target_logg=7.77
#target_teff= 4000.
#given_target_mass= 0.6
given_target_mass= 0.5


interp_kind='cubic'
#interp_kind='linear'

#target_logg=8.26
#target_teff= 4310. #K
##############################
def get_ms_lifetime(mass_wd):
    return 10*(8 *np.log(mass_wd/0.4))**(-2.5)


def operate_on_dist(dist1, dist2, function):
    """
    dist1 is the distribution of the first input to 'function'
    
    dist2 is the distribution of the second input to 'function'
    
    function is the most likely scinterp.interp2d() function output or object or whatever that you've defined that you'd 
    like to feed distributions through and be able to track which inputs yielded them.
    
    For whatever reason, scinterp.interp2d() -generated objects seem to completely randomize the indices of the 
    outputs relative to the inputs. I have no idea why. It also wants to take the 2 input (10,) shape arrays and make 
    a (10,10) output, whose indices have no correlation to the indices of either input (10,) array. I know, that 
    seems ridiculous, but it's what I've been experiencing, so I just made a damn for-loop.
    
    """
    output_dist= []
    for el1, el2 in zip(dist1, dist2):
        out_el= function(el1, el2)
        output_dist.append(out_el)
    output_dist=np.array(output_dist).T[0]
    return output_dist
##########################3

cooling_table= Table.read(cooling_model_file)

cooling_table.pprint()

loggteff_to_m = scinterp.interp2d(cooling_table['Teff'], cooling_table['logg'], cooling_table['Mass'], kind=interp_kind)
loggteff_to_logTc = scinterp.interp2d([cooling_table['Teff']], [cooling_table['logg']], cooling_table['Log(Tc)'], kind=interp_kind)

target_mass=  loggteff_to_m(target_teff, target_logg)
print('Target mass:', target_mass)

#target_mass_dist=loggteff_to_m(target_teff_dist, target_logg_dist)
target_mass_dist=operate_on_dist(target_teff_dist, target_logg_dist,loggteff_to_m)

print(target_teff_dist.shape)
print(np.where(target_mass_dist== loggteff_to_m(target_teff_dist[5], target_logg_dist[5])))
print('comparison', target_mass_dist[5], loggteff_to_m(target_teff_dist[5], target_logg_dist[5]))
print('target_mass_dist.shape', target_mass_dist.shape)
mean_mass= np.nanmean(target_mass_dist)
std_mass= np.std(target_mass_dist)
median_mass= np.nanmedian(target_mass_dist)
print('Mass:', mean_mass, '+/-', std_mass, 'or', median_mass)
plt.hist(target_mass_dist, bins=50)
plt.hist(simon_mass_dist, bins=50, alpha=0.5)
plt.axvline(x=mean_mass, color='r', linestyle='--')
plt.axvline(x=median_mass, color='g', linestyle= '--')
plt.xlabel('Mass')
plt.show()

teffm_to_age= scinterp.interp2d(cooling_table['Teff'], cooling_table['Mass'], cooling_table['Age'], kind=interp_kind)



target_age= teffm_to_age(target_teff, target_mass)

#target_age_dist=operate_on_dist(target_teff_dist, target_mass_dist, teffm_to_age)*1e-9 #Gyr units
target_age_dist=operate_on_dist(target_teff_dist, simon_mass_dist, teffm_to_age)*1e-9 #Gyr units
ms_age_dist=get_ms_lifetime(target_mass_dist)
total_age_dist= target_age_dist+ms_age_dist
clean_total_age_dist= np.copy(total_age_dist)
clean_total_age_dist[np.isnan(clean_total_age_dist)]=20.
clean_total_age_dist[np.where(clean_total_age_dist> 20.)] = 20. #setting a max
mean_age= np.nanmean(target_age_dist)
std_age= np.std(target_age_dist)
mean_total_age=np.nanmean(total_age_dist)
std_total_age=np.nanstd(total_age_dist)

med_total_age=np.nanmedian(clean_total_age_dist)
upper_total_age=np.nanpercentile(clean_total_age_dist, 84)
lower_total_age= np.nanpercentile(clean_total_age_dist, 16)

print('med total age', med_total_age, 'up to', upper_total_age, 'or down to', lower_total_age)
print('\n99\%\ chance that total age > ', np.nanpercentile(clean_total_age_dist, 1),'Gyr\n')
print('mean cooling age', mean_age, '+/-', std_age)
print('mean total age', mean_total_age, '+/-', std_total_age)
plt.hist(target_age_dist, bins=np.arange(0,21, 0.25),label='cooling ages', normed=True)
#plt.hist(ms_age_dist[~np.isnan(ms_age_dist)], bins=50, alpha=0.4, label='ms ages')
plt.hist(clean_total_age_dist, bins=np.arange(0,21, 0.25), alpha=0.5,label= 'Total Ages', normed=True)
plt.xlabel('Age (Gyr)')
plt.axvline(x=med_total_age, linestyle='--', color='k')
plt.axvline(x=upper_total_age, linestyle='--', color='k')
plt.axvline(x=lower_total_age, linestyle='--', color='k')
plt.legend()
#plt.yscale('log')
plt.show()

log_total= np.log10(total_age_dist)
clean_log_dist= np.copy(log_total)
clean_log_dist[np.isnan(clean_log_dist)]=20.
#clean_log_dist[np.where(clean_total_age_dist> 20.)] = 20.
plt.hist(clean_log_dist, bins=200, normed=True, color='g', alpha=0.5, label='Total Ages')
plt.xlabel('log10(age(Gyr))')
#plt.yscale('log')
plt.show()

target_age_gmass= teffm_to_age(target_teff, given_target_mass)
print('Target age:', target_age)
print('Target age assuming mass=', given_target_mass, ':', target_age_gmass)
target_logTc=  loggteff_to_logTc(target_teff, target_logg)
target_Tc= 10.** target_logTc

print('Target core temperature:', target_Tc)


loggteff_to_age= scinterp.interp2d(cooling_table['Teff'], cooling_table['logg'], cooling_table['Age'], kind=interp_kind)

target_age2= loggteff_to_age(target_teff, target_logg)
ms_lifetime= get_ms_lifetime(target_mass)

print("Target age from logg and teff:", target_age2)
print("MS lifetime:", ms_lifetime, "Gyr")
print("Total age from logg and teff:", ms_lifetime+(target_age2*1e-9))
print("Total age from given mass:", get_ms_lifetime(given_target_mass)+(target_age_gmass*1e-9), 'Gyr')

approx_inds= np.where((cooling_table['Teff']< 4000) & (cooling_table['Teff']> 3500))
approx_masses= cooling_table['Mass'][approx_inds]
approx_ages=cooling_table['Age'][approx_inds]*1e-9

wd_mass_vals= np.linspace(0.3, 1.3, 100)
cooling_ages=(teffm_to_age(target_teff, wd_mass_vals)*1e-9).T[0]
ms_ages= get_ms_lifetime(wd_mass_vals)
total_ages= cooling_ages+ms_ages
print('total_ages.shape', total_ages.shape)
print((teffm_to_age(target_teff, wd_mass_vals)*1e-9).shape)
print(get_ms_lifetime(wd_mass_vals).shape)
#plt.plot(wd_mass_vals, get_ms_lifetime(wd_mass_vals))
plt.axvline(x=0.5, linestyle='--', color='k')
plt.axhline(y=10, linestyle='--', color='k')
plt.plot(wd_mass_vals, total_ages, label='Total Age')
plt.plot(wd_mass_vals, cooling_ages, label='WD Cooling Age')
plt.plot(wd_mass_vals, ms_ages, label='MS lifetime')
plt.scatter(approx_masses, approx_ages, color='r', label='Grid vals with Teff ~3800K')
plt.xlabel('M_wd in solar masses')
#plt.ylabel('MS lifetime (Gyr)')
plt.ylabel('Age (Gyr)')
plt.legend()
#plt.yscale('log')
plt.ylim(0,15)
plt.show()

plt.scatter(cooling_table['Teff'], cooling_table['logg'], label='cooling models')
plt.plot(target_teff, target_logg, marker='*', color='r', markersize=12)
plt.xlabel('Teff (K)')
plt.ylabel('log(g)')
plt.title(cooling_model_file.split('/')[-1])
plt.show()


plt.scatter(cooling_table['Teff'], cooling_table['Mass'], label='cooling models')
plt.plot(target_teff, target_mass, marker='*', color='r', markersize=12)
plt.xlabel('Teff (K)')
plt.ylabel('Mass')
plt.title(cooling_model_file.split('/')[-1])

plt.show()


given_mass= given_target_mass
given_inds= np.where(cooling_table["Mass"]==given_mass)
plt.scatter(cooling_table['Teff'][given_inds], cooling_table['Age'][given_inds], label='table vals for '+ str(given_mass))
plt.scatter(cooling_table['Teff'][given_inds], teffm_to_age(cooling_table['Teff'][given_inds], given_mass), label='interpolated vals', color='r')
plt.plot(target_teff, target_age, marker='*', color='r', markersize=12, label='Target')
plt.plot(target_teff, target_age2, marker='*', color='g', markersize=12, label='Target (from teff+logg)')
plt.legend()
plt.xlabel('Teff')
plt.xlim(0,70000)
plt.ylabel('Age')
plt.show()


mass_vals= np.linspace(0.2,1.2, 100)
plt.plot(mass_vals, teffm_to_age(target_teff, mass_vals), label='interpolated vals for teff='+str(target_teff))
plt.plot(target_mass, target_age, marker='*', color='r', markersize=12, label='Target')
plt.legend()
plt.xlabel('mass')
plt.ylabel('age')
plt.show()

teff_vals= np.linspace(3000., 30000., 100)
allowed_inds= np.where(cooling_table['Mass']==given_target_mass)
plt.plot(teff_vals, loggteff_to_logTc(teff_vals, target_logg), label='log(Tc)')
plt.scatter(cooling_table['Teff'][allowed_inds], cooling_table['Log(Tc)'][allowed_inds], label='cooling table direct values for '+str(given_target_mass)+ 'M')
plt.plot(target_teff, target_logTc, marker='*', label='target')
plt.xlabel('Teff')

plt.ylabel('log(Tc)')
plt.legend()
plt.show()


logg_vals= np.linspace(7., 8.5, 100)
allowed_inds= np.where(cooling_table['Mass']==given_target_mass)
plt.plot(logg_vals, loggteff_to_logTc(target_teff, logg_vals), label='log(Tc)')
plt.scatter(cooling_table['logg'][allowed_inds], cooling_table['Log(Tc)'][allowed_inds], label='cooling table direct values for '+str(given_target_mass)+ 'M')
plt.plot(target_logg, target_logTc, marker='*', label='target')
plt.xlabel('logg')

plt.ylabel('log(Tc)')
plt.legend()
plt.show()








