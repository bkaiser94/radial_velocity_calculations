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


sens_names = glob('*sensitivity*.txt')

wavelengths = np.linspace(4940,8980, 8080) #400M2 approximately


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
plt.ylim(0,3)
spt.show_plot()
