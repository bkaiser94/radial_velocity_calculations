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
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const

cerro_pachon_location = coords.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)

def clean_filename(filename):
    file_parts=filename.split('.') #0 is object, 1 is ms, 2 is fits
    new_filename='.'.join([ 'ravg_wctb',file_parts[0]+'_gemini_600B', file_parts[-1]])
    return new_filename


def convert_file(filename):
    
    hdu=fits.open(filename)
    data_all=hdu[0].data
    header=hdu[0].header
    print('opening ',filename)
    try:
        print('inside try statement')
        if header['convertd']==True:
            print('\n\n',filename, 'is already converted to Goodman format. No conversion needed. Moving on.\n\n')
        else:
            print(filename, "'convertd' header does exist, but somehow the conversion hasn't been done... This really should never print. \n\n\n****")
    except KeyError:
        print("Key error caught")
        print('converting', filename)
        header.append(card=('convertd', True, 'Gemini to Goodman conversion performed'))
        try:
            flux_opt=data_all[0,0,:]
            wavelength=data_all[5,0,:]
            background=data_all[2,0,:]
            noise= data_all[3,0,:]
            new_filename=clean_filename(filename)
            #need to delete headers that are inaccurate in the new file, namely the things pointing to the labels of the spectrum
            for index in range(1,7):
                del header['BANDID'+str(index)]
            
            ### headers that need to be in the file for functions to work properly for flux calibration
            #AIRMASS, which is already present in the way we would want it.
            
            #OPENTIME
            header.append(card=('OPENTIME',header['UTSTART'],'time for shutter open in UTC'))
            #OPENDATE
            header.append(card=('OPENDATE',header['DATE-OBS'],'date for shutter open in UTC'))
            #PIX_SCAL pixel scale of CCD with the binning applied, i.e. if default scale is 0.15"/pix but we binned 2x2, this value should be 0.3
            header.append(card=('PIX_SCAL', header['PIXSCALE'],' "/pixel'))
            #SEE_SIG sigma value of the gaussian in units of binned pixels
            pixel_scale=float(header['PIXSCALE'])
            see_FWHM=header['SPECFWHM'] #I realize that without the ADC and since we positioned at the parallactic, this value is most likely technically slightly too large for the actual seeing, but given that it's going to be used on a DC for the flux calibration anyway, it's not like it really matters that it's slightly too large.
            seeing_sig=(see_FWHM*0.5)/np.sqrt(2*np.log(2.))
            print("Seing FWHM:", see_FWHM)
            print("Seeing Sigma:", seeing_sig)
            print('Seeing in "', see_FWHM*pixel_scale)
            header.append(card=("SEE_FWHM",see_FWHM,'FWHM of spectrum, approx. seeing in pixels'))
            header.append(card=('SEE_SIG', seeing_sig, 'Sigma value of gaussian fit to spectrum trace in pixels'))
            #SLIT the name of the slit that was used
            header.append(card=('SLIT', header['MASKNAME'],'slit name'))
            #EXPTIME this is already handled in the default headers of Gemini
            
            #I realized there were a number of headers that actually get called at the end of flux_calibration.py and are listed in "in_headers[]" in cal_params.py:
            
            #WIDTH there isn't a match to this since Erik does optimal extraction. I'll just fill in 7 as the default so the code runs
            header.append(card=('WIDTH', 7, "made-up width of extraction to allow code to run"))
            
            #SEE_FWHM already handled
            
            #ENVHUM
            header.append(card=('ENVHUM', header['HUMIDITY'], "Relative Humidity at start of exposure  "))
            
            #ENVPRE
            header.append(card=('ENVPRE', header['PRESSUR2']/100., "Atmospheric Pressure [hPS] at start of exposur "))
            #ENVTEM
            header.append(card=('ENVTEM', header['TAMBIENT'], "Outside Temperature [C] at start of exposure "))
            #ENVWIN
            header.append(card=('ENVWIN', header['WINDSPEE']*3.6, "Wind Speed [km/hr] at start of exposure "))
            #ENVDIR
            header.append(card=('ENVDIR', header['WINDDIRE'], "Wind Direction at start of exposure"))
            
            
            #I need-ish the BMJD_TDB version of the time in advance of calibrate_flux.py
            input_year = header['OPENDATE'] #gps-synched date
            input_hours = header['OPENTIME'] #gps-synched time
            exp_time= header['EXPTIME']*u.s
            input_times = input_year+'T'+input_hours #formatting correctly
            obs_time = Time(input_times, format = 'isot', scale = 'utc',location = cerro_pachon_location)
            obs_time= obs_time+exp_time/2.
            ra = header['RA']
            dec = header['DEC']
            target_coord = coords.SkyCoord(ra, dec, frame = 'icrs', unit= (u.hourangle, u.deg), )
            bary_corr =obs_time.tdb.light_travel_time(target_coord)
            bmjd_tdb_val = (obs_time.tdb+ bary_corr.tdb).mjd
            header.append(card = ('BMJD_TDB', bmjd_tdb_val, "exp. midpoint value from OPENDATE and OPENTIME headers"))
            ####
            dlambda=np.copy(np.roll(wavelength, -1) - wavelength)
            #now make the last wavelength behave correctly by just making it be the same as the second to last value
            dlambda[-1]=dlambda[-2]
            output_hdu=fits.PrimaryHDU(wavelength, header=header)
            hdu1=fits.ImageHDU(flux_opt)
            standin_ones=np.ones(flux_opt.shape)
            hdu2=fits.ImageHDU(background)
            hdu3=fits.ImageHDU(noise)
            hdu4=fits.ImageHDU(dlambda)
            hdulist=fits.HDUList([output_hdu,hdu1,hdu2,hdu3,hdu4])
            hdulist.writeto(new_filename, overwrite=True)
        except IndexError:
            print("Presumably this is actually a file that was converted, but was done prior to implementing the header to track it",filename,'\n\n')
    return








specfiles=glob.glob('*.ms.fits')
#specfiles=glob.glob('*.fits')
print('\n\n',specfiles,'\n\n')
for fname in specfiles:
    print(fname)
    print(clean_filename(fname))
    convert_file(fname)
    hdu=fits.open(fname)
    data_all=hdu[0].data
    header=hdu[0].header
    plt.figure()
    try:
        flux_opt=data_all[0,0,:]
        wavelength=data_all[5,0,:]
        plt.plot(wavelength,data_all[1,0,:],label='raw extracted spectrum')
    except IndexError as error:
        print("IndexError:", error)
        print('Not compatible with plotting part of code at this point.')

    plt.plot(wavelength,flux_opt, label='optimally extracted spectrum')
    plt.title(fname)
    plt.legend()
    plt.show()
