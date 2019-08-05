"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-06-06, D-Day

This should be a test-bed for the more complicated atmospheric correction processes that I'll be undertaking
I'm not sure if this is going to evolve into a separate script for atmospheric corrections or if it will be copied and 
pasted into flux_calibration.py and calibrate_flux.py; I just don't want to clutter those even further with
unnecessary plotting as I'm about to do here.

"""


from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const
from astropy import convolution as conv
import scipy.interpolate as scinterp
import csv
import time
start = time.time()

import spec_plot_tools as spt
import cal_params as cp
import get_cal_params as gcp

sens_names = glob('E*sensitivity_curve.txt')+glob('G*sensitivity_curve.txt')
#sens_names = glob('*sensitivity*.txt')

resid_names= glob('resid*eg274*sensitivity*.txt')

do_fnu=True

#wavelengths = np.linspace(4940,8980, 8080) #400M2 approximately
wavelengths = np.linspace(3800,7200, 8080) #400M1 approximately


def extract_AM_MJD(sens_curve_file):
    """
    INPUT: filename string for one of the sensitivity curve files
    
    OUTPUT: tuple of the airmass and MJD value for the sensitivity curve file from its 'header'
    
    """
    with open(sens_curve_file, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        #print('reader[0]',reader[0])
        index=0
        for row in reader:
            #print(row)
            if index==0:
                #print(row[0])
                #airmass_string=row[0][-5:]
                airmass_string= row[0].split(':')[1]
                #print('airmass:',airmass_string)
                #airmass_string=airmass_string.replace(' ', '')
                airmass= float(airmass_string)
                mjd_string= row[1].split(':')[1]
                mjd=float(mjd_string)
                #print('mjd:', mjd_string)
                #print('mjd:', mjd)
            else:
                pass
            index+=1
    
    return airmass, mjd

for sens_name in sens_names:
    airmass, mjd= extract_AM_MJD(sens_name)
    sens_curve_coeffs = np.genfromtxt(sens_name)
    sens_curve = np.polyval(sens_curve_coeffs,wavelengths)
    max_index= np.argmax(sens_curve)
    label= ','.join([sens_name, str(airmass), str(mjd)])
    plt.plot(wavelengths, sens_curve, label=label)
    plt.scatter(wavelengths[max_index], sens_curve[max_index], marker='*')
    print("\n=============")
    print(sens_name)
    print('airmass:', airmass, 'mjd:', mjd)
    print('Peak flux at', wavelengths[max_index], 'angstroms')
plt.xlabel(r'wavelength ($\AA$)')
plt.xlim(np.nanmin(wavelengths), np.nanmax(wavelengths))
plt.ylim(0,1)
spt.show_plot()

for resid_name in resid_names:
    airmass, mjd= extract_AM_MJD(resid_name)
    resid_array = np.genfromtxt(resid_name)
    max_index= np.argmax(sens_curve)
    label= ','.join([resid_name, str(airmass), str(mjd)])
    plt.plot(resid_array, label=label)
    print("\n=============")
    print(sens_name)
    print('airmass:', airmass, 'mjd:', mjd)
plt.xlabel(r'Pixel')
spt.show_plot()


def get_star_info(starname):
    standard_dict= cp.standard_dict[starname.lower()]
    standard_dict['filename']=cp.standard_dir+standard_dict['filename']
    return standard_dict

standard_name= 'GD153'
#standard_fits= 'avg_fwctb.GD153_400m2.fits'
standard_fits= 'avg_fwctb.GD153_400m1.fits'

standard_info = get_star_info(standard_name)
stand_array = np.genfromtxt(glob(standard_info['filename'])[0]).T
stand_waves1=stand_array[0]
stand_flux1=stand_array[1]
model_spec= np.vstack([stand_waves1, stand_flux1])

obs_spec, header, obs_noise= spt.retrieve_spec(standard_fits)
model_spec=spt.clean_spectrum(model_spec, np.nanmin(obs_spec[0]), np.nanmax(obs_spec[0]), [])
plt.plot(obs_spec[0], obs_spec[1]/np.nanmean(obs_spec[1]), label='obs')
if do_fnu:
    hdu= fits.open(standard_fits)
    dlambda= hdu[4].data
    obs_spec=spt.flambda_to_fnu(obs_spec, dlambda)
    model_spec[1]=model_spec[1]*1e16
    plt.plot(model_spec[0], model_spec[1]/np.mean(model_spec[1]), label='model flambda')

    model_spec= spt.flambda_to_fnu(model_spec)
    plt.ylabel(r'$f_{\nu}$')
else:
    obs_spec[1]=obs_spec[1]*10**-16




#model_spec=spt.clean_spectrum(model_spec, np.nanmin(obs_spec[0]), np.nanmax(obs_spec[0]), [])

plt.plot(obs_spec[0], obs_spec[1]/np.nanmean(obs_spec[1]), label='obs')
plt.plot(model_spec[0], model_spec[1]/np.nanmean(model_spec[1]), label='model')
plt.legend()
plt.title(standard_name)
plt.show()


