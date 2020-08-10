#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 11:55:47 2019

@author: edennihy and Ben Kaiser

Starting from the script Erik sent, this will take his format of reduced Gemini files and convert them to match the format of the Goodman spectral fits files that Ben Kaiser has made all of his code run with.

Major things that will need to be addressed are as follows:

- filenames need to match Ben's convention

- header needs to be retained and some fields probably renamed or generated (i.e. see_sig)

- extension storage needs to match up with Ben's convention

- probably have to add in the 
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import glob


def clean_filename(filename):
    file_parts=filename.split('.') #0 is object, 1 is ms, 2 is fits
    new_filename='.'.join([ 'ravg_wctb',file_parts[0]+'_gemini_600B', file_parts[-1]])
    return new_filename


def convert_file(filename):
    
    hdu=fits.open(filename)
    data_all=hdu[0].data
    header=hdu[0].header
    flux_opt=data_all[0,0,:]
    wavelength=data_all[5,0,:]
    new_filename=clean_filename(filename)
    output_hdu=fits.PrimaryHDU(wavelength, header=header)
    hdu1=fits.ImageHDU(flux_opt)
    standin_ones=np.ones(flux_opt.shape)
    hdu2=fits.ImageHDU(standin_ones)
    hdu3=fits.ImageHDU(standin_ones)
    hdu4=fits.ImageHDU(standin_ones)
    hdulist=fits.HDUList([output_hdu,hdu1,hdu2,hdu3,hdu4])
    hdulist.writeto(new_filename, overwrite=True)
    return








specfiles=glob.glob('*.ms.fits')
for fname in specfiles:
    print(fname)
    print(clean_filename(fname))
    convert_file(fname)
    hdu=fits.open(fname)
    data_all=hdu[0].data
    header=hdu[0].header
    flux_opt=data_all[0,0,:]
    wavelength=data_all[5,0,:]

    plt.figure()
    plt.plot(wavelength,flux_opt)
    plt.title(fname)
    plt.show()
