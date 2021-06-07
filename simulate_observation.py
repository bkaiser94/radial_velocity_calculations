"""
Created by Ben Kaiser (UNC-Chapel Hill) 2021-06-06


Take a model spectrum (probably provided by Simon Blouin) and convolve it with the seeing, then truncate by
the slit width. Then apply noise pixel-by-pixel.

Output a simulated observation with the specified signal-to-noise.

First pass at this will not recalibrate the SNR for different parts of the spectrum. It's just going to apply the
same level of noise at all wavelengths regardless of flux level.

Later if I care enough, I'll mess with it to do neat stuff like flux-dependent noise levels.

Maybe it will eventually become an actual exposure time calculator for Goodman... but probably not.




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
from astropy.table import Table, Column
import scipy.interpolate as scinterp
import time
start = time.time()



import spec_plot_tools as spt
import cal_params as cp
import model_manipulation as mm
import plot_spec as ps






#wd_dir='WDJ2317p1830/'
#model_file='J2317_KCa-0.5.txt'
#model_file='J2317_KCa0.5.txt'
#model_file='J2317_KCa0.0.txt'
#model_file='J2317_KCa-1.0.txt'


wd_dir='WDJ1824p1213/'
model_file='J1824_KCa1.7.txt'
#model_file='J1824_KCa2.2.txt'
#model_file='J1824_KCa2.7.txt'



reference_observation='ravg_fwctb.GaiaJ1644m0449_20190825_tellcorr_400m2.fits'

obs_spec, header, obs_noise= spt.retrieve_spec(reference_observation)

wd_dir=cp.model_spectra_dir+wd_dir


model_file=wd_dir+model_file



snr_per_pixel=20.
#pixel_width=1.5 #angstroms of pixel in dispersion
#slit_width=3.0 #slit width in arc seconds



#def degrade_model(model_vals, obs_vals, header):
    #"""
    #Convolve model spectrum with the seeing and rebin it (using flux-conservative method) to the pixel-scale of
    #the observation.
    
    #INPUTS:
        #model_vals - [wavelengths, fluxes, wavelength bin widths]
        #obs_vals - [wavelengths, fluxes, wavelength bin widths]
        #header - observation header, which should include the pixel width, slit width, seeing, etc.
        
    #OUTPUTS:
        #model_spec 
    
    #"""
    #rebinned_model_spec= spt.rebin_generic_spec(model_vals, model_vals[1][1]-model_vals[1][0], obs_vals[0], obs_vals[2])
    #slit_width = spt.get_slit_width(header)
    #rebinned_model_spec= mm.convolve_model_new(rebinned_model_spec, header, slit_width=slit_width)
    #output_spec= rebinned_model_spec
    #return output_spec



model_array=np.genfromtxt(model_file).T
print(model_array[0])

#model_array[1]=model_array[1]*const.c/model_array[0]**2
#models from Simon are in f_nu it seems, but I don't have the conversion to go from f_nu to f_lambda readily available, so... I'll just test it in f_nu I guess

deg_model_spec=mm.degrade_model(model_array,obs_spec,header)

noisy_flux=np.random.normal(loc=deg_model_spec[1],scale=deg_model_spec[1]/snr_per_pixel)
noisy_spec=np.vstack([deg_model_spec[0],noisy_flux])
smoothed_spec=ps.convolve_spectrum(noisy_spec,header,kernel_type='gaussian')

plt.title(model_file.split('/')[-1])
plt.plot(noisy_spec[0],noisy_spec[1],label='model with resoution degraded and pixel SNR '+ str(snr_per_pixel))
plt.plot(smoothed_spec[0],smoothed_spec[1],label='smoothed version of noisy spec')
plt.plot(model_array[0], model_array[1], label='model')
plt.plot(deg_model_spec[0][200:-200],deg_model_spec[1][200:-200], label='model with resolution degraded')



plt.legend()
plt.xlabel('Wavelength (Angstrom)')
plt.xlim([7400,8000])
plt.show()
















