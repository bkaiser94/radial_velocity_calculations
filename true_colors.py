"""
Created by Benjamin C. Kaiser (UNC-Chapel Hill) 2021-06-03

Please cite me if you use this for images with something like: "made using true_colors from 
Benjamin C. Kaiser"

**** This is not a "science" code. I make no guarantees that this actually is the real color you 
would see if looking at a given star or blackbody. This code is a good faith attempt to do so,
and I am reasonably confident that the end result is roughly the color a person with normal
color vision would see. Again, it could be wrong; I did it in like 5 days.

***** I was inspired to make this code after seeing this PheT

Take an optical spectrum and generate an RGB circle that shows the "true" color of the
spectrum. This could probably be generalized to spectra that do not cover the full visible 
range such that the code infers the slope at those other parts, but initially it's just going to to 
work for my Goodman data of J1644, and it will probably actually use the stitched spectrum 
that I used to make the figure 1 of the Science paper.

The spectrum will be rescaled to the CIE 1931 XYZ system
(https://en.wikipedia.org/wiki/CIE_1931_color_space ;wikipedia link because it's plain-ish 
language, but I didn't use the math here), using linear interpolation between the grid points of the same CIE 1931 XYZ curves as in Wyman et al. 2013
(http://jcgt.org/published/0002/02/01/). CIE 1931 XYZ curves were retrieved from https://www.rit.edu/cos/colorscience/rc_useful_data.php


I'm then using skimage.color.xyz2rgb() to convert these XYZ colors to the sRGB color space, 
which I'm pretty sure is the endzone for the color transformation.

Then I just have to make a circle with that color. I'm 80% sure this last step will inexplicably be 
the most difficult for me based on my extreme difficulties with simple tasks *shrug*.

There are parts of this code that I left in deliberately to show how to edit it for contrived 
spectral inputs to make interesting colors. Basically, I would just make flat spectra at 0 with 
sections with flux=1 to see what colors I could make. 


***** These are not hard and fast requirements for the code to run, but they are the versions I used:

Python 3.6.5
astropy 3.2.3
skimage 0.14.2
numpy 1.19.3
"""
#import astropy
#import skimage
#import numpy

#print(astropy.__version__)
#print(skimage.__version__)
#print(numpy.__version__)

#from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from astropy import units as u
from astropy.table import Table, Column


import skimage.color
import skimage.data

#import spec_plot_tools as spt
#import cal_params as cp


#input_file='stitched_J1644_spectrum.fits'
#input_file='ravg_fwctb.sdssj1501m0816_400m1.fits'

#wd_name='WD J1644-0449'

default_bb_teff=5800.

#bb_teff=5800.
#bb_teff=11000.
#bb_teff=3830.

#input_file='stitched_J2356_spectrum.fits'
#wd_name='WD J2356-209'
#bb_teff=4040.

#input_file='ravg_fwctb.WISEA0615m1247_400m1.fits'
#wd_name='WISEA J0615-1247'
#bb_teff=3300.

#input_file='ravg_fwctb.sdssj1501m0816_400m1.fits'
#wd_name='SDSS J1501-0816'
#bb_teff=4500.

#visible_range=[340., 720.] #wavelengths in nm that should even be considered
#visible_range=[600., 660.] #wavelengths in nm that should even be considered
#visible_range=[380., 780.] #wavelengths in nm that should even be considered
visible_range=[375.,725.]



CIE1931_XYZ_file='RIT_CIE_XYZ_1931_colorfunctions.csv'
#CIE1931_XYZ_file=cp.true_color_dir+CIE1931_XYZ_file

CIE1931_XYZ_table=Table.read(CIE1931_XYZ_file)
#CIE1931_XYZ_table.pprint()

#default_CIE_method='single_lobe'
default_CIE_method='interp'


def x1931(wave, method=default_CIE_method):
    """
    from equation 2 of Wyman et al. 2013
    
    wave should be in nanometers, I'm 99% sure; I didn't see Wyman et al. explicitly state to do 
    so, but they discuss wavelengths in nm throughout the paper, so it's not a crazy leap.
    
    """
    if method=='single_lobe':
        return 1.065*np.exp(-0.5*((wave-595.8)/33.33)**2)+0.366 * np.exp(-0.5 * ((wave-446.8)/19.44)**2)
    elif method=='interp':
        return np.interp(wave, CIE1931_XYZ_table['Wavelength (nm)'], CIE1931_XYZ_table['xbar'])

def y1931(wave, method=default_CIE_method):
    """
    from equation 2 of Wyman et al. 2013
    
    wave should be in nanometers, I'm 99% sure; I didn't see Wyman et al. explicitly state to do 
    so, but they discuss wavelengths in nm throughout the paper, so it's not a crazy leap.
    
    """
    if method=='single_lobe':
        return  1.014 * np.exp(-0.5*((np.log(wave)-np.log(556.3))/0.075)**2)
    elif method =='interp':
        return np.interp(wave, CIE1931_XYZ_table['Wavelength (nm)'], CIE1931_XYZ_table['ybar'])



def z1931(wave, method=default_CIE_method):
    """
    from equation 2 of Wyman et al. 2013
    
    wave should be in nanometers, I'm 99% sure; I didn't see Wyman et al. explicitly state to do 
    so, but they discuss wavelengths in nm throughout the paper, so it's not a crazy leap.
    
    """
    if method=='single_lobe':
        return 1.839 * np.exp(-0.5*((np.log(wave)-np.log(449.8))/0.051)**2)
    elif method=='interp':
        return np.interp(wave, CIE1931_XYZ_table['Wavelength (nm)'], CIE1931_XYZ_table['zbar'])

def get_CIE_val(spec, dwave, CIEfunc, method=default_CIE_method):
    match_vals=CIEfunc(spec[0], method=method)
    radiance_vals=spec[1]*np.pi
    product_vals= match_vals*radiance_vals*dwave
    return np.sum(product_vals)


def get_CIE1931_XYZ(spec, dwave, method=default_CIE_method):
    X=get_CIE_val(spec, dwave, x1931, method=method)
    Y=get_CIE_val(spec, dwave,y1931, method=method)
    Z=get_CIE_val(spec,dwave,z1931, method=method)
    
    return [X, Y, Z]


### I'm going to assume that the input spectrum will be in units of f_lambda because that's the spectrum I'll be using in the near-term.

#the stellar true-color PHet simulation at the link on the next line gives the spectral radiance conversion from spectral power density, which is apparently just to multiply by pi, so not ground-breaking.
#https://github.com/phetsims/blackbody-spectrum/blob/4d8b70c78ee85210d7b79661c65c5ba8732c8273/js/blackbody-spectrum/model/BlackbodyBodyModel.js

def planck_function(wave, temperature):
    """
    Just using the form given in the Phet-stored values because it's easier
    
    returns Watts per square meter per micron with wavelength input in nanometers
    
    """
    constA= 3.74192e-16
    constB=1.438770e7
    return constA / ((wave**5)*  (np.exp(constB/(wave * temperature))-1))



#input_spec,header, noise=spt.retrieve_spec(input_file)
#hdu=fits.open(input_file)
#delta_wave=np.copy(hdu[4].data)

#input_spec[0]=input_spec[0]-650. #shifting wavelength values to try to line up absorption with color.

#######


def trim_spec(input_spec, min_wave, max_wave):
    lower_indices = np.where(input_spec[0]< max_wave)
    trimmed_waves= input_spec[0][lower_indices]
    trimmed_flux= input_spec[1][lower_indices]
    upper_indices= np.where(trimmed_waves > min_wave)
    trimmed_waves= trimmed_waves[upper_indices]
    trimmed_flux = trimmed_flux[upper_indices]
    trimmed_spec= np.vstack([trimmed_waves, trimmed_flux])
    #print trimmed_spec.shape
    return trimmed_spec

def remove_range(input_spec, bound_list):
    """
    Removes wavelengths and flux values from the array that fall in the range specified by bound_list
    """
    wave_array= input_spec[0]
    other_array= input_spec[1]
    lower_bound = bound_list[0]
    upper_bound= bound_list[1]
    low_mask = np.where(wave_array < lower_bound)
    high_mask= np.where(wave_array > upper_bound)
    low_waves= wave_array[low_mask]
    high_waves= wave_array[high_mask]
    low_other = other_array[low_mask]
    high_other= other_array[high_mask]
    merge_waves= np.append(low_waves, high_waves)
    merge_other= np.append(low_other, high_other)
    return np.vstack([merge_waves, merge_other])


def sort_spectrum(input_spec):
    """
    Fixes the ordering of the wavelengths (and associated fluxes) so that the go from least to greatest, thereby removing errant lines through the plotted spectra and allowing for correct interpolation in model-fitting.
    """
    sort_indices = np.argsort(input_spec[0])
    sorted_waves= input_spec[0][sort_indices]
    sorted_flux= input_spec[1][sort_indices]
    sorted_spectrum= np.vstack([sorted_waves, sorted_flux])
    
    return sorted_spectrum



def clean_spectrum(input_spec, min_wave, max_wave, mask_list):
    """
    input_spec should be a vstack of wavelengths and the flux (or error)
    """
    clean_spec= trim_spec(input_spec, min_wave, max_wave)
    for mask in mask_list:
        clean_spec= remove_range(clean_spec, mask)
    clean_spec= sort_spectrum(clean_spec)
    return np.copy(clean_spec)

def get_true_colors(input_file='', bb_teff=default_bb_teff):
    """
    Inputs:
    
    input_file - string of a csv file (comma-delimited) with columns containing in order wavelengths (in angstroms), flux (in f_lambda; ideally erg/cm^2/s/angstrom, but only is vital to be f_lambda), pixel widths in angstroms. The pixel widths in angstroms should probably be smaller than 5 angstroms for the sake of accuracy (meaning the "wavelengths" should also be within 5 angstroms of each other more or less.
    
    
    bb_teff - temperature in Kelvin of the blackbody curve to also be 
    
    
    Returns:
    
    star_rgb - a list of the sRGB colors of the integrated spectrum with values between 0 and 1
    
    ***** This code does not detect gaps in your spectrum nor does it care if you don't fully cover the visible range 375-720 nm, but if your spectrum does not cover that range or has gaps, your true color output WILL NOT BE CORRECT (most likely)! It's possible by some fluke that a missing portion of spectrum 
    
    ******** The input spectrum MUST BE FLUX-CALIBRATED if you want this to at all be correct. The relative fluxes are what make colors.
    
    
    
    """
    
    wd_name=input_file
    if ((input_file =='') or (input_file==' ')):
        print('\n\nNo input_file specificied', input_file,'\n\n')
        print('Input spectrum will be replaced in future steps with a flat spectrum covering the visible range as a result')
        input_wavelengths=np.arange(visible_range[0]-10.,visible_range[1]+10.,1.0)
        delta_wave=np.ones(input_wavelengths.shape)
        input_flux=np.ones(input_wavelengths.shape)
        input_spec=np.vstack([input_wavelengths, input_flux])
        
    else:
        input_array=np.genfromtxt(input_file, delimiter=',').T
        input_spec=np.vstack([input_array[0],input_array[1]])
        delta_wave=np.copy(input_array[2])
    input_spec[0]=input_spec[0]*0.1 #converting angstroms to nm
    delta_wave=delta_wave*0.1 ##converting angstroms to nm
    input_spec[1]=input_spec[1]*10. #converting to nm^-1

    delta_wave_spec= np.vstack([np.copy(input_spec[0]),delta_wave])
    #vis_delta_wave_spec=spt.clean_spectrum(delta_wave_spec,visible_range[0],visible_range[1],[])
    #vis_input_spec=spt.clean_spectrum(input_spec, visible_range[0],visible_range[1],[])
    vis_delta_wave_spec=clean_spectrum(delta_wave_spec,visible_range[0],visible_range[1],[])
    vis_input_spec=clean_spectrum(input_spec, visible_range[0],visible_range[1],[])

    #def get_CIE_val(spec, dwave, CIEfunc):
        #match_vals=CIEfunc(spec[0])
        #radiance_vals=spec[1]*np.pi
        #product_vals= match_vals*radiance_vals*dwave
        #return np.sum(product_vals)

    Xval=get_CIE_val(vis_input_spec,vis_delta_wave_spec[1],x1931)
    Yval=get_CIE_val(vis_input_spec,vis_delta_wave_spec[1],y1931)
    Zval=get_CIE_val(vis_input_spec,vis_delta_wave_spec[1],z1931)

    print('X', Xval)
    print('Y', Yval)
    print('Z', Zval)
    star_XYZ=np.array([Xval, Yval, Zval])
    renormed_XYZ=star_XYZ/np.max(star_XYZ)
    print('renormed_XYZ',renormed_XYZ)

    interp_XYZ=np.array(get_CIE1931_XYZ(vis_input_spec,vis_delta_wave_spec[1], method='interp'))
    renormed_interp_XYZ=interp_XYZ/np.max(interp_XYZ)

    print('interp XYZ', interp_XYZ)
    print('renormed_interp_XYZ', renormed_interp_XYZ)



    #blackbody_3800K=planck_function(vis_input_spec[0], 3800.)
    #blackbody_3800K=planck_function(vis_input_spec[0], 10000.)
    blackbody_3800K=planck_function(vis_input_spec[0], bb_teff)
    
    ### section where I messed around with contrived spectra to make neat colors #####

    #blackbody_3800K[:]=1.
    #blackbody_3800K[:]=0.

    #blackbody_3800K[np.where((vis_input_spec[0]>400.) & (vis_input_spec[0] < 500.))]=1.
    #blackbody_3800K[np.where((vis_input_spec[0]>550.) & (vis_input_spec[0] < 660.))]=1.

    #blackbody_3800K[np.where((vis_input_spec[0]>400.) & (vis_input_spec[0] < 600.))]=1.

    #blackbody_3800K[np.where((vis_input_spec[0]>350.) & (vis_input_spec[0] < 450.))]=1.

    ###### end section for messing around with contrived spectra (except for later places where display labels would need to be changed) ########


    bb_spec=np.vstack([vis_input_spec[0], blackbody_3800K])

    bb_Xval=get_CIE_val(bb_spec,vis_delta_wave_spec[1],x1931)
    bb_Yval=get_CIE_val(bb_spec,vis_delta_wave_spec[1],y1931)
    bb_Zval=get_CIE_val(bb_spec,vis_delta_wave_spec[1],z1931)

    bb_XYZ=np.array([bb_Xval,bb_Yval,bb_Zval])
    renormed_bb_XYZ=bb_XYZ/np.max(bb_XYZ)
    print('renormed_bb_XYZ', renormed_bb_XYZ)

    bb_rgb=skimage.color.xyz2rgb([[renormed_bb_XYZ]])


    #astronaut=skimage.data.astronaut()
    #print(astronaut.shape)
    #print(astronaut)
    star_rgb=skimage.color.xyz2rgb([[renormed_XYZ]])

    print('Input spectrum RGB', star_rgb)
    print(star_rgb.shape)
    print('Input spectrum RGB 255', star_rgb*255)
    #print('RGB 130', star_rgb*130)


    print('BB RGB', bb_rgb)
    print('BB RGB 255', bb_rgb*255)
    #print('BB RGB 130', bb_rgb*130)

    norm_blackbody_3800K=blackbody_3800K/np.max(blackbody_3800K)
    norm_vis_input=vis_input_spec[1]/np.max(vis_input_spec[1])

    plot_waves=np.arange(visible_range[0],visible_range[1],1.)
    #plt.plot(plot_waves,x1931(plot_waves, method='single_lobe'), label='x1931 single_lobe', color='r')
    #plt.plot(plot_waves,y1931(plot_waves, method='single_lobe'), label='y1931 single_lobe', color='g')
    #plt.plot(plot_waves,z1931(plot_waves, method='single_lobe'), label='z1931 single_lobe', color='b')
    plt.plot(plot_waves,x1931(plot_waves, method='interp'), label='x1931 interp', color='r', linestyle=':')
    plt.plot(plot_waves,y1931(plot_waves,method='interp'), label='y1931 interp',color='g', linestyle=':')
    plt.plot(plot_waves,z1931(plot_waves, method='interp'), label='z1931 interp', color='b', linestyle=':')

    plt.plot(vis_input_spec[0],vis_input_spec[1]/np.max(vis_input_spec[1]),label=wd_name)
    #plt.plot(vis_input_spec[0],vis_input_spec[1]/np.max(vis_input_spec[1]),label=wd_name+  'with a -650 angstrom shift')

    plt.plot(vis_input_spec[0], norm_blackbody_3800K, label=str(bb_teff)+' K Blackbody')
    #plt.plot(vis_input_spec[0], norm_blackbody_3800K, label='contrived purple spectrum')
    plt.scatter(500, 2, c=star_rgb[0],edgecolors='k',s=160)
    plt.scatter(500, 1.8, c=bb_rgb[0],edgecolors='k',s=160)
    plt.text(510,1.8, str(bb_teff)+' K Blackbody True Color')
    #plt.text(510,1.8, 'Contrived True Color')

    plt.text(510,2,wd_name+' True Color')
    plt.xlabel('Wavelength (nm)')
    
    #plt.text(510,2,wd_name+'  Color for -650 angstrom shift')


    #for single_val in plot_waves:
        #single_spec=np.vstack([[single_val],[1]])
        #this_X=get_CIE_val(single_spec,vis_delta_wave_spec[1][0],x1931,method='single_lobe')
        #this_Y=get_CIE_val(single_spec,vis_delta_wave_spec[1][0],y1931,method='single_lobe')
        #this_Z=get_CIE_val(single_spec,vis_delta_wave_spec[1][0],z1931,method='single_lobe')
        #this_XYZ=np.array([this_X,this_Y,this_Z])
        #renormed_this_XYZ=this_XYZ/np.max(this_XYZ)
        #this_rgb=skimage.color.xyz2rgb([[renormed_this_XYZ]])
        #plt.scatter(single_val, -0.25, c=this_rgb[0])
        
    #plt.text(720,-0.25,'single_lobe')

    for single_val in plot_waves:
        single_spec=np.vstack([[single_val],[1]])
        this_X=get_CIE_val(single_spec,vis_delta_wave_spec[1][0],x1931,method='interp')
        this_Y=get_CIE_val(single_spec,vis_delta_wave_spec[1][0],y1931, method='interp')
        this_Z=get_CIE_val(single_spec,vis_delta_wave_spec[1][0],z1931, method='interp')
        this_XYZ=np.array([this_X,this_Y,this_Z])
        renormed_this_XYZ=this_XYZ/np.max(this_XYZ)
        this_rgb=skimage.color.xyz2rgb([[renormed_this_XYZ]])
        plt.scatter(single_val, -0.40, c=this_rgb[0])

    plt.text(720,-0.40,'interp')

    plt.legend()
    plt.show()

    #plt.plot(vis_input_spec[0],x1931(vis_input_spec[0])*norm_vis_input*vis_delta_wave_spec[1],label='x1931 cross spec')
    #plt.plot(vis_input_spec[0],y1931(vis_input_spec[0])*norm_vis_input*vis_delta_wave_spec[1],label='y1931 cross spec')
    #plt.plot(vis_input_spec[0],z1931(vis_input_spec[0])*norm_vis_input*vis_delta_wave_spec[1],label='z1931 cross spec')
    ##plt.legend()
    ##plt.show()

    #plt.plot(vis_input_spec[0],x1931(vis_input_spec[0])*1*vis_delta_wave_spec[1],label='x1931 cross 1')
    #plt.plot(vis_input_spec[0],y1931(vis_input_spec[0])*1*vis_delta_wave_spec[1],label='y1931 cross 1')
    #plt.plot(vis_input_spec[0],z1931(vis_input_spec[0])*1*vis_delta_wave_spec[1],label='z1931 cross 1')

    #plt.plot(vis_input_spec[0],x1931(vis_input_spec[0])*norm_blackbody_3800K*vis_delta_wave_spec[1],label='x1931 cross norm BB '+ str(bb_teff)+' K')
    #plt.plot(vis_input_spec[0],y1931(vis_input_spec[0])*norm_blackbody_3800K*vis_delta_wave_spec[1],label='y1931 cross norm BB '+ str(bb_teff)+' K')
    #plt.plot(vis_input_spec[0],z1931(vis_input_spec[0])*norm_blackbody_3800K*vis_delta_wave_spec[1],label='z1931 cross norm BB '+ str(bb_teff)+' K')
    #plt.legend()
    #plt.show()
    
    return star_rgb





#############################3
#############################


if __name__ =='__main__':
    
    input_filename=sys.argv[1]
    print('input_filename:', input_filename)
    try:
        input_bb_teff= float(sys.argv[2])
    except IndexError:
        print("No Blackbody Teff provided as second command line argument.")
        print("Using default Blackbody Teff", default_bb_teff)
        input_bb_teff=default_bb_teff
        
    get_true_colors(input_file=input_filename, bb_teff=input_bb_teff)
    
    








