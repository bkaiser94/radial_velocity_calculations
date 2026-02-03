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

import plot_spec as ps
import spec_plot_tools as spt

plt.rc('lines',linewidth=0.5)
norm_range=[8000,10000]

count=1
input_files=glob('*orders*fits')
for input_file in input_files:
    hdu=fits.open(input_file)
    print(hdu)
    print(len(hdu))
    print(hdu[0].data)
    spec=hdu[0].data
    #if spec[0][0]<100:
        #spec[0]=spec[0]*10**4.
        #spec[1]=spec[1]*10.**16
    spec=ps.norm_spectrum(spec,[1.04,1.104],show_norm_range=True)
    #spec=ps.norm_spectrum(spec,norm_range,show_norm_range=True)
    plt.plot(spec[0],spec[1],label=input_file,alpha=1)
    #plt.xlabel('Wavelength (microns)')
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux (cgs units)')
    #plt.ylabel('Flux_lambda (normalized)')
    count+=1
plt.title('TripleSpec Orders 4 through 7 (Order 3 had to be truncated in the merged version)')
plt.legend()
plt.show()
m2_spec_file=glob('*400m2*.fits')[0]
m2_spec,m2_header,m2_noise=spt.retrieve_spec(m2_spec_file)
m2_spec[1]=m2_spec[1]*10.**(-16.)
m2_spec[0]=m2_spec[0]*10.**(-4.)
ps.plot_spectrum(m2_spec,m2_spec_file,m2_header,norm=False,norm_range=norm_range)
plt.legend()
plt.show()
#spt.show_plot()







count=1
print('\n\nStarting orders plotting.\n\n')
order_files=glob('*J0212*xtc*fits')
print(input_files)
for input_file in order_files:
    hdu=fits.open(input_file)
    print(hdu)
    print(len(hdu))
    print(hdu[0].data)
    for thing in hdu[0].data:
        order_spec=thing
        #if spec[0][0]<100:
            #spec[0]=spec[0]*10**4.
            #spec[1]=spec[1]*10.**16
        #spec=ps.norm_spectrum(spec,[1.04,1.104],show_norm_range=True)
        #spec=ps.norm_spectrum(spec,norm_range,show_norm_range=True)
        plt.plot(order_spec[0],order_spec[1],alpha=1)
        count+=1
    plt.xlabel('Wavelength (microns)')
    #plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux (cgs units)')
    #plt.ylabel('Flux_lambda (normalized)')
    plt.plot(spec[0],spec[1],label=input_files[0],alpha=1)
    ps.plot_spectrum(m2_spec,m2_spec_file,m2_header,norm=False,norm_range=norm_range)

plt.title('J0212-5522A All orders hopefully')
plt.legend()
plt.show()






