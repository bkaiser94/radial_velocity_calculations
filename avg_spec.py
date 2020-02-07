"""
Created by Ben Kaiser 2019-03-27 (UNC-Chapel Hill)

@author: Ben Kaiser

Should average together the counts and uncertainties of matching extracted spectra.

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
import time
start = time.time()

import spec_plot_tools as spt

#glob_string= 'wctb*'
glob_string='fwctb*'
#filenames= glob('wctb*')
filenames= glob(glob_string)

#low_index=10
low_index= len(glob_string)-1+6
high_index=-5

do_dlambda_ext= True
do_super=False
filename_matches={}
prev_core=''
filename_set=[]
#filename_set={}
print(filenames)

def get_core_name(filename):
    return  filename[low_index:high_index]
    #return filename.split('_')[1]

for i,filename in enumerate(filenames[:-1]):
    #core_name= filename[low_index:high_index]
    core_name= get_core_name(filename)
    #print('filename:', filename, 'core_name:', core_name)
    filename_set.append(filename)
    if core_name != filenames[i+1][low_index:high_index]:
        #filename_set.append(filename)
        #print('match')
        #print(core_name, prev_core)
        try:
            #print("filename_matches[core_name]", filename_matches[core_name])
            prev_list=filename_matches[core_name]
            #print('prev_list', prev_list)
            new_list=prev_list+filename_set
            filename_matches[core_name]=new_list
            #filename_matches[core_name].append(filename_set)
            #filename_matches[core_name]=filename_matches[core_name].append(filename_set)
            #print("filename_matches[core_name]", filename_matches[core_name])
        except KeyError:
            print('No previous core names matching:', core_name)
            print('previous core_names:')
            for thing in filename_matches:
                print(thing)
            filename_matches[core_name]=filename_set
            #if len(filename_set) > 1:
                #filename_matches[core_name]=filename_set
            #else:
                #filename_matches[core_name]=filename_set[0]
            #print("filename_matches[core_name]", filename_matches[core_name])
        filename_set=[]
        print("resetting filename_set")
    else:
        pass
    
    #else:
        #try:
            #print(filename_matches[prev_core])
            #filename_matches[prev_core].append(filename_set)
        #except KeyError:
            #filename_matches[prev_core]=filename_set
        #filename_matches.append([filename_set])
        #print(filename_set)
        #filename_set=[]
        #filename_set.append(filename)
    #prev_core= core_name
    #filename_set.append(filename)
    print('======')
    print(filename)
    print(core_name)
    print(filename_set)
    print('+++++++')
#filename_set.append(filename)
final_core_name= get_core_name(filenames[-1])
print('final_core_name', final_core_name)
#filename_set.append(filenames[-1]) #this assumed the last file matched the preceding file, which isn't always the case
#filename_matches.append([filename_set])
for sets in filename_matches:
    print('======')
    print(sets)
    print('+++++++')
try:
    print(filename_matches[final_core_name])
    prev_list = filename_matches[final_core_name]
    print('prev_list', prev_list)
    print('filenames[-1]', filenames[-1])
    if core_name == final_core_name:
        filename_set.append(filenames[-1])
        prev_list= prev_list+filename_set
    else:
        prev_list.append(filenames[-1])
    print('prev_list', prev_list)
    #filename_matches[final_core_name].append(filename_set)
    filename_matches[final_core_name]=prev_list
    print('filename_matches[final_core_name]', filename_matches[final_core_name])
except KeyError:
    #filename_matches[final_core_name]=filename_set
    print("last filename doesn't have a previously stored set, which probably means it's the last one of a set as used to be usual.")
    print('No previous core names matching:', final_core_name)
    print('previous core_names:')
    for thing in filename_matches:
        print(thing)
    print('filename_set', filename_set)
    filename_set.append(filenames[-1])
    filename_matches[final_core_name]=filename_set
for sets in filename_matches:
    print('======')
    print(filename_matches[sets])
    print('+++++++')
    
#print(filename_matches['eg274_400m2'])

#######################################3
#########################################
##NOw do the actual image combination with those lists of filenames that can be used for combination.



def avg_spectra(target_list):
    wave_list=[]
    flux_list=[]
    noise2_list=[]
    sky_list=[]
    if do_dlambda_ext:
        dlambda_list=[]
    else:
        pass
    print("target_list", target_list)
    for target_file in target_list:
        print('target_file', target_file)
        target_spec, header, target_noise = spt.retrieve_spec(target_file)
        wave_list.append(target_spec[0])
        flux_list.append(target_spec[1])
        noise2_list.append(target_noise[1]**2)
        hdu=fits.open(target_file)
        sky=np.copy(hdu[2].data)
        sky_list.append(sky)
        global do_dlambda_ext
        if do_dlambda_ext:
            try:
                dlambda_spec= np.copy(hdu[4].data)
                print('dlambda_spec extension (4) is present')
                dlambda_list.append(dlambda_spec)
            except IndexError as error:
                print("IndexError:", error)
                print('No dlambda_spec extension presumably, so not doing it for any other spectra for this target')
                do_dlambda_ext=False
        else:
            print('do_dlambda_ext=False')
            
        #sky_list(
        plt.plot(target_spec[1], label=target_file)
    flux_array=np.array(flux_list)
    wave_array=np.array(wave_list)
    print(flux_array.shape)
    noise2_array= np.array(noise2_list)
    sky_array=np.array(sky_list)
    std_wave=np.std(wave_array,axis=0)
    avg_std_wave=np.nanmean(std_wave)
    max_std_wave=np.nanmax(std_wave)
    print('average standard deviation of wavelength values', avg_std_wave)
    print('max standard deviation of wavelength values', max_std_wave)
    avg_wave= np.nanmean(wave_array,axis=0)
    avg_flux= np.nanmean(flux_array, axis=0)
    avg_noise2= np.sum(noise2_array, axis=0)/noise2_array.shape[0]**2
    avg_noise= np.sqrt(avg_noise2)
    avg_sky= np.nanmean(sky_array,axis=0)
    if do_dlambda_ext:
        dlambda_array= np.array(dlambda_list)
        avg_dlambda= np.nanmean(dlambda_array, axis=0)
    plt.plot(avg_flux, label='avg')
    #plt.title(target_file)
    plt.legend()
    plt.show()
    core_name=target_file[low_index:high_index]
    avg_noise=avg_noise/avg_flux
    output_name= "avg_"+ glob_string[:-1]+'.'+ core_name + ".fits"
    #output_name= "avg_"+ core_name + ".fits"
    print('output_name:', output_name)
    print('saving hopefully...')
    hdu=fits.PrimaryHDU(avg_wave, header= header)
    hdu1= fits.ImageHDU(avg_flux)
    #hdu2= fits.ImageHDU(np.ones(avg_flux.shape))
    hdu2= fits.ImageHDU(avg_sky)
    hdu3 = fits.ImageHDU(avg_noise)
    if do_dlambda_ext:
        hdu4=fits.ImageHDU(avg_dlambda)
        hdulist= fits.HDUList([hdu,hdu1,hdu2,hdu3, hdu4])
    else:
        hdulist= fits.HDUList([hdu,hdu1,hdu2,hdu3])
    #hdulist= fits.HDUList([hdu,hdu1,hdu2,hdu3])
    hdulist.writeto(output_name, overwrite= True)
    print('should have saved')
    return


def make_super_spec(target_list, save=True, plot_all=True):
    wave_list=[]
    flux_list=[]
    noise2_list=[]
    sky_list=[]
    if do_dlambda_ext:
        dlambda_list=[]
    else:
        pass
    print("target_list", target_list)
    for target_file in target_list:
        print('target_file', target_file)
        target_spec, header, target_noise = spt.retrieve_spec(target_file)
        wave_list.append(target_spec[0])
        flux_list.append(target_spec[1])
        noise2_list.append(target_noise[1]**2)
        hdu=fits.open(target_file)
        sky=np.copy(hdu[2].data)
        sky_list.append(sky)
        core_name=target_file[low_index:high_index]
        global do_dlambda_ext
        if do_dlambda_ext:
            try:
                dlambda_spec= np.copy(hdu[4].data)
                print('dlambda_spec extension (4) is present')
                dlambda_list.append(dlambda_spec)
            except IndexError as error:
                print("IndexError:", error)
                print('No dlambda_spec extension presumably, so not doing it for any other spectra for this target')
                do_dlambda_ext=False
        else:
            print('do_dlambda_ext=False')
            
        #sky_list(
        plt.plot(target_spec[0],target_spec[1], label=target_file)
    flux_array=np.array(flux_list)
    flux_array=np.copy(flux_array.ravel())
    wave_array=np.array(wave_list)
    wave_array=np.copy(wave_array.ravel())
    #print(flux_array.shape).ravel()
    noise2_array= np.array(noise2_list)
    noise2_array=np.copy(noise2_array.ravel())
    sky_array=np.array(sky_list)
    sky_array=np.copy(sky_array.ravel())
    sorted_indices= np.argsort(wave_array)
    sorted_waves=wave_array[sorted_indices]
    sorted_flux= flux_array[sorted_indices]
    sorted_sky= sky_array[sorted_indices]
    sorted_noise2= noise2_array[sorted_indices]
    sorted_noise=np.sqrt(sorted_noise2)
    if do_dlambda_ext:
        dlambda_array= np.array(dlambda_list)
        dlambda_array=np.copy(dlambda_array.ravel())
        sorted_dlambda= dlambda_array[sorted_indices]
    if plot_all:
        plt.plot(sorted_waves,sorted_flux, label='super')
        #plt.title(target_file)
        plt.legend()
        plt.show()
    else:
        pass
    header.append(card=('SUPER', True, 'This is a spectrum that just merges all of the spectra without averaging'))
    hdu=fits.PrimaryHDU(sorted_waves, header= header)
    hdu1= fits.ImageHDU(sorted_flux)
    #hdu2= fits.ImageHDU(np.ones(avg_flux.shape))
    hdu2= fits.ImageHDU(sorted_sky)
    hdu3 = fits.ImageHDU(sorted_noise)
    if do_dlambda_ext:
        hdu4=fits.ImageHDU(sorted_dlambda)
        hdulist= fits.HDUList([hdu,hdu1,hdu2,hdu3, hdu4])
    else:
        hdulist= fits.HDUList([hdu,hdu1,hdu2,hdu3])
    
    if save:
        output_name= "super_"+ glob_string[:-1]+'.'+ core_name + ".fits"
        print('output_name:', output_name)
        print('saving hopefully...')
        hdulist.writeto(output_name, overwrite= True)
        print('should have saved')
    else:
        print('save=', save)
    return hdulist
    
    
def make_rebin_avg_spec(target_list):
    """
    combine all of the spectra into the same wavelength bins using a weighted average method.
    
    This method requires knowing dlambda, so I'm just going to let it crash if there isn't a dlambda present
    
    """
    #super_hdulist= make_super_spec(target_list, save=False, plot_all=False)
    #ref objects are the ones that will be used to establish the baseline for the rebinned average spectrum
    ref_spec, ref_header, ref_noise= spt.retrieve_spec(target_list[0])
    ref_hdu=fits.open(target_list[0])
    ref_sky=np.copy(ref_hdu[2].data)
    core_name=target_list[0][low_index:high_index]
    ref_dlambda= np.copy(ref_hdu[4].data)
    rebin_fluxes=[]
    rebin_skies= []
    rebin_noise2s=[]
    for target_file in target_list:
        print('target_file', target_file)
        rebin_list= spt.rebin_spec(target_file, ref_spec[0], ref_dlambda)
        rebin_fluxes.append(rebin_list[1])
        rebin_skies.append(rebin_list[2])
        rebin_noise2s.append(rebin_list[3])
    
    avg_rebin_flux= np.nanmean(rebin_fluxes, axis=0)
    avg_rebin_sky= np.nanmean(rebin_skies, axis=0)
    avg_rebin_noise2= np.nanmean(rebin_noise2s, axis=0)**2/np.nansum(rebin_noise2s, axis=0)
    avg_rebin_noise= np.sqrt(avg_rebin_noise2)
    avg_rebin_noise=avg_rebin_noise/avg_rebin_flux
    hdu=fits.PrimaryHDU(ref_spec[0], header= ref_header)
    hdu1= fits.ImageHDU(avg_rebin_flux)
    #hdu2= fits.ImageHDU(np.ones(avg_flux.shape))
    hdu2= fits.ImageHDU(avg_rebin_sky)
    hdu3 = fits.ImageHDU(avg_rebin_noise)
    hdu4=fits.ImageHDU(ref_dlambda)
    hdulist= fits.HDUList([hdu,hdu1,hdu2,hdu3, hdu4])
    output_name= "ravg_"+ glob_string[:-1]+'.'+ core_name + ".fits"
    print('output_name:', output_name)
    print('saving hopefully...')
    hdulist.writeto(output_name, overwrite= True)
    print('should have saved')
    return

print("\n\n\n\n")
print("##################")
print("##################")
print("##################")
print("\n\n\n\n")

for sets in filename_matches:
    print('======')
    print(filename_matches[sets])
    do_dlambda_ext=True
    try:
        avg_spectra(filename_matches[sets])
    except ValueError as error:
        print(error)
        print('\nskipping the lame averaging for ', filename_matches[sets],'\n')
    try:
        make_rebin_avg_spec(filename_matches[sets])
    except IndexError as error:
        print(error)
        print('\nskipping', filename_matches[sets],'\n')
    #make_super_spec(filename_matches[sets])
    if do_super:
        make_super_spec(filename_matches[sets])
    else:
        pass
    print('+++++++')







