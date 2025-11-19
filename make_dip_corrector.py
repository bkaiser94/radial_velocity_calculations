"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-11-18

This script should generate the dip-correction interpolation function that can then be used in another script to actually do the corrections. Hopefully without having to do the fitting every single time.

With an assist from Copilot on how to save interpolations using pickles.

Ok well copilot is pretty adamant that the more robust method is to actually construct the interpolator every time I want to use it... the problem with that is that my rebinning stuff doesn't work that way... I suppose I could just save the rebinned spectra and then load those right?


"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
import scipy.interpolate as scinterp
#import cosmics
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const
from astropy.modeling import models as asmodels
from astropy.modeling import fitting as asfitting
from astropy.table import Table, Column
#import pickle

import get_cal_params as gcp
import cal_params as cp
import spec_plot_tools as spt
import plot_spec as ps



filename='ravg_fwctb.WDJ1223m1852_400m1.fits'
input_dir='misc_spectra/'

rebinning_width=30. #width of new bins of rebinned spectrum to use for the interpolation in angstroms
bb_teff=8600. #Teff of blackbody in K that should be used for correcting the DC spectrum. This should be its effective temperature.

dip_region=[4705,5050]

input_dir=cp.ref_dir+input_dir


filename_list=['ravg_fwctb.WDJ1223m1852_400m1.fits',
              'ravg_fwctb.WDJ1218m3307_400m1.fits',
              'ravg_fwctb.WDJ1241m1335_400m1.fits'
    
    
    ]

bb_teff_list=[8780.,
              8633.,
              8254.
    
    
    ]# All taken from examining the entry for each white dwarf in MWDD. 

#########3

#The rebinned spectra seem to be placing their wavelength centers on the left edges of the bin instead of in the middle !!!!!!
###########

def make_interpolator(filename=filename,bb_teff=bb_teff):
    target_spec, header, target_noise=spt.retrieve_spec(input_dir+filename)
    #Do extinction correction
    #target_spec= spt.correct_extinction(target_spec, header, plot_all=True)

    hdu=fits.open(input_dir+filename)
    ps.plot_spectrum(target_spec, filename, header, smooth=False, norm=False, pix_width=5, kernel_type='box')
    #get lower bin centers
    low_bin=np.nanmin(target_spec[0])+rebinning_width
    high_bin=np.nanmax(target_spec[0])-rebinning_width
    extra_waves=(np.nanmax(target_spec[0])-np.nanmin(target_spec[0]))%rebinning_width
    #Now make the bin centers run up to some value
    cut_low=6600.
    #target_spec=spt.clean_spectrum(target_spec, np.nanmin(target_spec[0]),np.nanmax(target_spec[0]),[[cut_low,cut_low+extra_waves]])
    bin_centers_low=np.arange(np.nanmin(target_spec[0])+rebinning_width/2.,cut_low,rebinning_width)
    bin_centers_high=np.arange(np.nanmax(target_spec[0])-rebinning_width/2.,cut_low+rebinning_width,-1*rebinning_width)
    #print('bin_centers_high',bin_centers_high)
    #bin_centers=np.arange(np.nanmin(target_spec[0])+0.5*rebinning_width,np.nanmax(target_spec[0]),rebinning_width)
    bin_centers=np.append(bin_centers_low,bin_centers_high)
    #print(bin_centers)
    binned_spec=spt.rebin_generic_spec(target_spec,hdu[4].data,bin_centers,np.ones(bin_centers.shape)*rebinning_width)
    plt.scatter(binned_spec[0],binned_spec[1],label='original binned')
    #Clean the binned_spec
    #binned_spec=spt.clean_spectrum(binned_spec,3740, 7060,[[6850,6950]]) #trimming the first bin that is missing light for some reason and trimming out the strong O2 band in the 400M1 spectrum
    binned_spec=spt.clean_spectrum(binned_spec,3700, 7100,[[6850,6950]]) # trimming out the strong O2 band in the 400M1 spectrum
    #yakima=scinterp.Akima1DInterpolator(target_spec[0],target_spec[1])(target_spec[0])
    smoothed_spec=ps.convolve_spectrum(target_spec,header,kernel_type='box',pix_width=rebinning_width/hdu[4].data[3], kernel_width=50.)
    plt.scatter(binned_spec[0],binned_spec[1],label='trimmed binned')
    #yakima=scinterp.Akima1DInterpolator(binned_spec[0],binned_spec[1])(target_spec[0])
    akima_interpolator=scinterp.Akima1DInterpolator(binned_spec[0],binned_spec[1])
    sakima_interpolator=scinterp.Akima1DInterpolator(smoothed_spec[0],smoothed_spec[1])
    yakima=akima_interpolator(target_spec[0])
    sakima=sakima_interpolator(target_spec[0])
    ps.plot_spectrum(binned_spec, '30 A-binned spec', header, smooth=False, norm=False, pix_width=5, kernel_type='box')
    ps.plot_spectrum(smoothed_spec, '30 A-smoothed spec', header, smooth=False, norm=False, pix_width=5, kernel_type='box')
    plt.plot(target_spec[0],yakima,label='Akima Interp')
    plt.plot(target_spec[0],sakima,label='smoothed Akima Interp',marker='o')
    bb_spec=spt.blackbody_spec(target_spec[0],bb_teff)
    ####print(bb_spec)  
    ps.plot_spectrum(bb_spec, str(bb_teff)+' K BB', header, norm=True, smooth=False)
    spt.show_plot(show_legend=True, line_id='cool_wd', convert_to_air=True,actually_show=False)
    spt.show_plot(show_legend=False, line_id='h', convert_to_air=True, actually_show=True)

    norm_bb_spec=ps.norm_spectrum(bb_spec,ps.norm_range,show_norm_range=False)
    norm_interp_spec=ps.norm_spectrum(np.vstack([target_spec[0],yakima]),ps.norm_range,show_norm_range=False)

    plt.plot(norm_bb_spec[0],norm_bb_spec[1],label='normed BB spec')
    plt.plot(norm_interp_spec[0],norm_interp_spec[1],label='normed interp spec')
    plt.xlabel('Wavelength Angstroms')
    plt.legend()
    plt.show()

    div_flux=norm_bb_spec[1]/norm_interp_spec[1]
    plt.plot(target_spec[0],div_flux,label='norm bb spec/norm interp spec')
    plt.xlabel('Wavelength Angstroms')
    plt.legend()
    plt.show()


    corrected_spec=np.vstack([target_spec[0],target_spec[1]*div_flux])
    ps.plot_spectrum(corrected_spec, 'Hopefully corrected w/ leftovers', header, smooth=False, norm=True, pix_width=5, kernel_type='box')
    ps.plot_spectrum(target_spec, filename, header, smooth=False, norm=False, pix_width=5, kernel_type='box')
    #plt.show()
    spt.show_plot(show_legend=True, line_id='h', convert_to_air=True, actually_show=True)

    too_many_points=np.linspace(1000,10000,90000)
    plt.plot(too_many_points, akima_interpolator(too_many_points),label='unhinged interpolation')
    ps.plot_spectrum(norm_interp_spec, filename, header, norm=False, smooth=False)
    spt.show_plot(show_legend=True, line_id='h', convert_to_air=True, actually_show=True)
    #print('\n\n\n\n\n=================\ndiv_flux in function',div_flux)
    plt.plot(div_flux)
    plt.show()
    return target_spec[0],div_flux

#make_interpolator()
div_list=[]
for filename, bb_teff in zip(filename_list, bb_teff_list):
    new_waves,new_div=make_interpolator(filename=filename, bb_teff=bb_teff)
    #print('\n\n\nnew_div')
    #print(new_div)
    #print(type(new_div))
    div_list.append([new_waves,new_div])
    #wave_list.append(new_waves)
    
div_array=np.array(div_list)
#wave_array=np.array(wave_list)
print('div_array.shape', div_array.shape)
for num,entry in enumerate(div_array):
    plt.plot(entry[0],entry[1],label=str(num))

avg_div=np.mean(div_array,axis=0)
med_div=np.median(div_array,axis=0)
plt.plot(div_array[0,0],avg_div[1],label='avg')
plt.plot(div_array[0,0],med_div[1],label='median')
plt.legend()
plt.show()

good_inds=np.where(~np.isnan(med_div[1]))
print(med_div.shape)
sub_div=med_div[:,good_inds[0]]
for value in sub_div[0]:
    print(value)
np.savetxt(input_dir+'dip_correction_points.csv',sub_div)

new_akima=scinterp.Akima1DInterpolator(med_div[0],med_div[1])

