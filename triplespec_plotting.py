"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-12-16



Plot the output triplespec spectra from Spextool that I had to run in IDL on Maytag.





"""
from __future__ import print_function


#import matplotlib
#matplotlib.use('pdf')


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
from astropy.table import Table, Column
import scipy.interpolate as scinterp
import time
start = time.time()

plt.rc('lines',linewidth=0.5)

count=1
input_files=glob('NLTT*merge*')
for input_file in input_files:
    hdu=fits.open(input_file)
    print(hdu)
    print(len(hdu))
    print(hdu[0].data)
    spec=hdu[0].data
    plt.plot(spec[0],spec[1],label='run '+str(count))
    plt.xlabel('Wavelength (microns)')
    plt.ylabel('Flux (cgs units)')
    count+=1
plt.title('NLTT2478 TripleSpec Orders 4 through 7 (Order 3 had to be truncated in the merged version)')
plt.legend()
plt.show()
