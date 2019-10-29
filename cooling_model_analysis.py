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

target_logg=7.77
target_teff= 3830.

#target_logg=7.77
#target_teff= 4000.
given_target_mass= 0.6


#target_logg=8.26
#target_teff= 4310. #K


cooling_table= Table.read(cooling_model_file)

cooling_table.pprint()

loggteff_to_m = scinterp.interp2d(cooling_table['Teff'], cooling_table['logg'], cooling_table['Mass'], kind='cubic')
loggteff_to_logTc = scinterp.interp2d([cooling_table['Teff']], [cooling_table['logg']], cooling_table['Log(Tc)'], kind='cubic')

target_mass=  loggteff_to_m(target_teff, target_logg)
print('Target mass:', target_mass)

teffm_to_age= scinterp.interp2d(cooling_table['Teff'], cooling_table['Mass'], cooling_table['Age'], kind='cubic')

target_age= teffm_to_age(target_teff, target_mass)
target_age_gmass= teffm_to_age(target_teff, given_target_mass)
print('Target age:', target_age)
print('Target age assuming mass=', given_target_mass, ':', target_age_gmass)
target_logTc=  loggteff_to_logTc(target_teff, target_logg)
target_Tc= 10.** target_logTc

print('Target core temperature:', target_Tc)


loggteff_to_age= scinterp.interp2d(cooling_table['Teff'], cooling_table['logg'], cooling_table['Age'], kind='cubic')

target_age2= loggteff_to_age(target_teff, target_logg)
print("Target age from logg and teff:", target_age2)

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


given_mass= 0.5
given_inds= np.where(cooling_table["Mass"]==given_mass)
plt.scatter(cooling_table['Teff'][given_inds], cooling_table['Age'][given_inds], label='table vals for '+ str(given_mass))
plt.scatter(cooling_table['Teff'][given_inds], teffm_to_age(cooling_table['Teff'][given_inds], given_mass), label='interpolated vals')
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
allowed_inds= np.where(cooling_table['Mass']==0.5)
plt.plot(teff_vals, loggteff_to_logTc(teff_vals, target_logg), label='log(Tc)')
plt.scatter(cooling_table['Teff'][allowed_inds], cooling_table['Log(Tc)'][allowed_inds], label='cooling table direct values for 0.5M')
plt.plot(target_teff, target_logTc, marker='*', label='target')
plt.xlabel('Teff')

plt.ylabel('log(Tc)')
plt.legend()
plt.show()


logg_vals= np.linspace(7., 8.5, 100)
allowed_inds= np.where(cooling_table['Mass']==0.5)
plt.plot(logg_vals, loggteff_to_logTc(target_teff, logg_vals), label='log(Tc)')
plt.scatter(cooling_table['logg'][allowed_inds], cooling_table['Log(Tc)'][allowed_inds], label='cooling table direct values for 0.5M')
plt.plot(target_logg, target_logTc, marker='*', label='target')
plt.xlabel('logg')

plt.ylabel('log(Tc)')
plt.legend()
plt.show()








