"""
Created by Ben Kaiser (UNC-Chapel Hill) 2024-03-11.

Should open a bunch of FITS files and plot their targeted coordinates both in RA/Dec and ObsRA/ObsDec. We'll see what the heck is going on with them I guess.

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




#file_base='*LTT*besselR*.fits'
#file_base='*WDJ1500*.fits'
#file_base='*GD153*im*.fits'
#file_base='*GD71*besselR*.fits'
file_base='*EG274*.fits'


filenames=glob(file_base)
filenames=sorted(filenames)

ra_list=[]
dec_list=[]
obsra_list=[]
obsdec_list=[]

offset_list=['','','','','','','offset','offset','','','offset','','offset']
#target_coords=coords.SkyCoord(ra='08:41:30.7', dec='-32:56:04.5', unit=(u.hourangle, u.degree)) #LTT3218
#target_coords=coords.SkyCoord(ra='12:57:02.37', dec='+22:01:56.0', unit=(u.hourangle, u.degree)) #GD153
#target_coords=coords.SkyCoord(ra='05:52:27.51', dec='+15:53:16.6', unit=(u.hourangle, u.degree)) #GD71

target_coords=coords.SkyCoord(ra='16:23:33.75', dec=' -39:13:47.5', unit=(u.hourangle, u.degree)) #EG274


    
def plot_coordinates(file_list):
    ra_vals=[]
    dec_vals=[]
    obsra_vals=[]
    obsdec_vals=[]
    time_vals=[]
    rotator_vals=[]
    num_list=[]
    lastairmass=0
    lasttime=''
    for i, filename in enumerate(file_list):
        header=fits.getheader(filename)
        #coordinates=coords.SkyCoord(ra=header['ra'], dec=header['dec'], unit=(u.hourangle, u.degree))
        #obscoordinates=coords.SkyCoord(ra=header['obsra'], dec=header['obsdec'], unit=(u.hourangle, u.degree))
        #print(header['ra'],header['dec'])
        #print(coordinates.ra.arcsecond, coordinates.dec.arcsecond)
        #ra_vals.append(coordinates.ra.arcsecond)
        #dec_vals.append(coordinates.dec.arcsecond)
        print(filename[:4])
        num_list.append(filename[:4])
        ra_vals.append(header['ra'])
        dec_vals.append(header['dec'])
        obsra_vals.append(header['obsra'])
        obsdec_vals.append(header['obsdec'])
        #time_vals.append(header['BMJD_TDB'])
        lastairmass=header['airmass']
        lasttime=header['opentime']
        if header['rotator']> 180:
            rotator_val=header['rotator']-360
        else:
            rotator_val=header['rotator']
        rotator_vals.append(rotator_val)
    print('ra_vals',ra_vals)
    coordinates=coords.SkyCoord(ra=ra_vals,dec=dec_vals,unit=(u.hourangle, u.degree))
    obscoordinates=coords.SkyCoord(ra=obsra_vals, dec=obsdec_vals, unit=(u.hourangle, u.degree))
    print('coordinates.ra', coordinates.ra.degree)
    print(num_list)
    num_list=np.array(num_list)
    #num_list=np.str_(num_list)
    print('num_list',num_list)
    print(len(num_list), len(offset_list))
    plt.plot(coordinates.ra.degree, coordinates.dec.degree, label='RA/Dec',marker='o',color='b',alpha=0.3)
    plt.plot(obscoordinates.ra.degree, obscoordinates.dec.degree, label='ObsRA/ObsDec',marker='o',color='r',alpha=0.3)
    for index in range(0,len(coordinates)):
        #print(num_list[index],offset_list[index])
        #plt.text(coordinates[index].ra.degree, coordinates[index].dec.degree+(np.random.rand()-0.5)*.0002, num_list[index]+offset_list[index],color='b',alpha=0.3)
        #plt.text(obscoordinates[index].ra.degree, obscoordinates[index].dec.degree+(np.random.rand()-0.5)*.0002, num_list[index]+offset_list[index],color='r',alpha=0.3)
        plt.text(coordinates[index].ra.degree, coordinates[index].dec.degree+.00001, num_list[index],color='b',alpha=0.3)
        plt.text(obscoordinates[index].ra.degree, obscoordinates[index].dec.degree+.00001, num_list[index],color='r',alpha=0.3)
        
    plt.plot(target_coords.ra.degree, target_coords.dec.degree,label='Target List Coords', marker='*', color='k')
    plt.ylabel('Dec in decimal degrees')
    plt.xlabel('RA in decimal degrees')
    #plt.title('LTT 3218 Alignment images from 2024-03-04')
    #plt.title('GD71 Alignment images from 2021-01-09')
    #plt.title('WDJ1500-1603 On target spectra from both nights')
    #plt.title('EG274 Alignment images from 2023-02-19')
    #title_string='EG274 Alignment images from 2023-09-13'
    #title_string='EG274 Alignment images from 2023-09-14'
    title_string='EG274 Alignment images from 2023-09-13&14'
    
    plt.title(title_string+', airmass '+str(lastairmass)+', at ' +lasttime)

    plt.legend()
    plt.show()
    ra_array=np.array(ra_vals) 
    dec_array=np.array(dec_vals)
    delta_ra=ra_array-ra_array[0]
    delta_dec=dec_array-dec_array[0]
    for number, thing in enumerate(delta_ra):
        print(number, thing)
        plt.text(delta_ra[number], delta_dec[number], str(number))
    plt.plot(delta_ra, delta_dec, marker='o', color='b')
    plt.xlabel(r'$\Delta$ RA (arcseconds)')
    plt.ylabel(r'$\Delta$ Dec (arcseconds)')
    plt.show()
    return


plot_coordinates(filenames)
