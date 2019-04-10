"""
Created by Ben Kaiser (UNC-Chapel Hill) 201902-26

This should be able to read in the astropy tables output by model_fitting.py that contain the best-fit model logg 
and Teff and their associated errors, and then it should plot them all in the model space together with the mean
value and its uncertainties from the addition of the other errors in quadrature.

"""
from __future__ import print_function
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



#import wdatmos
import spec_plot_tools as spt
#import kernel_builder
import model_manipulation as mm

######################################3

min_teff= 6000
max_teff= 14750
min_logg= 3.75
max_logg= 6.25

def get_mean(input_table):
    mean_logg= np.average(input_table['logg'], weights= 1/(input_table['logg_error']**2)) #weighted-average of the values
    mean_teff= np.average(input_table['teff'], weights=1/(input_table['teff_error']**2))
    mean_logg_error= np.sqrt(np.sum(input_table['logg_error']**2)/input_table['logg_error'].shape[0]**2)
    mean_teff_error= np.sqrt(np.sum(input_table['teff_error']**2)/input_table['teff_error'].shape[0]**2)
    print('Mean Values:')
    print('Teff:', mean_teff, '+/-', mean_teff_error, 'logg:', mean_logg, '+/-', mean_logg_error)
    print('std Teff:', np.std(input_table['teff']), 'std logg', np.std(input_table['logg']))
    print('unweighted mean values', np.mean(input_table['teff']), np.mean(input_table['logg']))
    return [ mean_teff, mean_teff_error,mean_logg, mean_logg_error]

########################################3

#input_file='20190212_new_model_fits.csv'
#input_file='20190303_model_fits_nomask.csv'
input_file= '20190408_model_fits_nomask.csv'

print("input_file:", input_file)
input_table = Table.read(input_file, format='ascii.csv')

mean_list= get_mean(input_table)
plt.axvline(x=min_teff, color='k', linestyle=':')
plt.axvline(x=max_teff, color='k', linestyle=':')
plt.axhline(y=min_logg, color='k', linestyle=':')
plt.axhline(y=max_logg, color='k', linestyle=':')
plt.errorbar(input_table['teff'], input_table['logg'], xerr= input_table['teff_error'], yerr=input_table['logg_error'], linestyle='none', marker='o')
plt.errorbar(mean_list[0], mean_list[2], xerr=mean_list[1], yerr=mean_list[3], linestyle='none', color='r', marker='o')
plt.xlabel(r'$T_{eff}$')
plt.ylabel('log(g)')
plt.show()


plt.errorbar(input_table['bmjd_tdb'], input_table['rv'], yerr=input_table['rv_error'], linestyle='none', marker='o')
plt.xlabel('BMJD_TDB')
plt.ylabel('Radial Velocity (km/s)')
plt.show()

plt.hist(input_table['teff'])
plt.xlabel('Teff')
plt.show()

plt.hist(input_table['logg'])
plt.xlabel('logg')
plt.show()


plt.errorbar(input_table['bmjd_tdb'], input_table['teff'], yerr=input_table['teff_error'], linestyle='none', marker='o')
plt.xlabel('BMJD_TDB')
plt.ylabel(r'$T_{eff} (K)$')
plt.show()

plt.errorbar(input_table['bmjd_tdb'], input_table['logg'], yerr=input_table['logg_error'], linestyle='none', marker='o')
plt.xlabel('BMJD_TDB')
plt.ylabel('log(g)')
plt.show()

plt.errorbar(input_table['rv'], input_table['teff'], yerr=input_table['teff_error'], xerr=input_table['rv_error'], linestyle='none', marker='o')
plt.xlabel('RV (km/s)')
plt.ylabel(r'$T_{eff} (K)$')
plt.show()
