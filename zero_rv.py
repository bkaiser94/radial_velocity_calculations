"""
Created by Ben Kaiser (UNC-Chapel Hill) 06-30-2018 

This is supposed to take all of the reduced target spectra and in combination with the model radial velocities from
phase_folding.py, it should then take all of the spectra to a radial velocity of zero. That is the hope anyway... I 
am not certain that I am going about this the correct way.

"""




import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
import astropy.constants as const
from astropy.units import cds
import scipy.stats as scistats
import scipy.optimize as sciop
cds.enable()

import spec_plot_tools as spt


listfile = 'listFWCTB'
target_list = np.genfromtxt(listfile, dtype ='str')
#output_list = []
#output_list_name = 'listZFWCTB' #zero rv is the new letter tacked onto the front. I'm pretty sure the only thing we should be changing is the wavelength values in the zero index anyway
model_rv_input = np.genfromtxt('model_rvs.txt', names = True)


model_rvs = model_rv_input['RV_kms']


def rv_correct(wavelengths, radial_velocity):
    return (wavelengths*const.c.si/((radial_velocity*(u.km/u.s)).si+const.c.si)).value


count = 0
for rv, target_file in zip(model_rvs, target_list):
    new_name = 'z'+target_file
    #output_list.append(new_name)
    i=fits.open(target_file)
    header = fits.getheader(target_file)
    file_waves= i[0].data
    file_flux = i[1].data
    file_bkg_flux = i[2].data
    file_noise = i[3].data
    file_waves = rv_correct(file_waves, rv)
    header.append(card = ('RVused', rv, '(km/s) (used to change wavelengths)'))
    hdu=fits.PrimaryHDU(file_waves, header = header)
    hdu1= fits.ImageHDU(file_flux)
    hdu2 = fits.ImageHDU(file_bkg_flux)
    hdu3 = fits.ImageHDU(file_noise)
    hdulist = fits.HDUList([hdu, hdu1, hdu2, hdu3])
    hdulist.writeto(new_name, overwrite = True)
    plt.plot(file_waves, file_flux+count, label = target_file)
    #plt.plot(file_waves, file_flux+count, label = target_file, linestyle = 'none', marker = '.')
    #plt.plot(file_waves, linestyle = 'none', marker = 'o')
    count+=1
    
#np.savetxt(output_list_name, output_list, fmt= '%.|S36')

#plt.legend()
plt.xlabel(r'Wavelength $(\AA)$')
plt.ylabel(r'Flux +N (ergs/cm/cm/s/A 10**-16)')
plt.show()



    
    
    
