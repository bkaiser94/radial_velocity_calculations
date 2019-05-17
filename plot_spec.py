"""
Created 2019-01-29 by Ben Kaiser (UNC-Chapel Hill)

@author: Ben Kaiser

This should just be able to plot a given spectrum and various things about it.

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



#print start
import wdatmos
import spec_plot_tools as spt

slit_width = 1.0 #arcseconds
pixel_scale = 0.3 #arcseconds per pixel_scale
slit_width = slit_width/pixel_scale #slit width in pixels
test_wavelength = 4686
test_width = 40
test_side = test_width/2

sdss_pix_width = 10

#wavelength_offset=60
wavelength_offset=0

#filenames = glob(sys.argv[1])
#filenames= glob('wctb*')
#filenames= glob('wctb*SDSS*')
sdss_path = '/Users/BenKaiser/Desktop/SDSS_speclib/'
#print(filenames)
plot_wavelength=True
plot_400m2_tell= False
#norm_range=[1240,1280]
#norm_range=[1560,1590]
#norm_range=[40,80]

#norm_range=[7042,7046]
#norm_range=[7490,7510]
norm_range=[7517,7556]
#norm_range=[8074,8140]
#norm_range=[8058,8231]
#norm_range=[5100,5400]

#norm_range=np.array(norm_range)+wavelength_offset

#file_setting='all_avg'
#file_setting='command'
#file_setting='all_wctb'
file_setting= 'compare_SDSS'

file_setting= 'all_SDSS'



if file_setting=='all_avg':
    print(file_setting)
    #filenames=glob('avg_*')
    filenames=glob('avg_fwctb*fits')
    #single_iterate=False
    #double_iterate=True
    single_iterate=True
    double_iterate=False

elif file_setting=='all_wctb':
    print(file_setting)
    filenames=glob('wctb*')
    single_iterate=True
    double_iterate=False

elif file_setting =='command':
    filename1=sys.argv[1]
    filename2=sys.argv[2]
    #filename1= glob(sys.argv[1])
    #filename2=glob(sys.argv[2])
    print(sys.argv[1])
    print(filename1)
    print(sys.argv[2])
    print(filename2)
    single_iterate=False
    double_iterate=False
    
elif file_setting=='compare_SDSS':
    filename=sys.argv[1]
    sdss_names = glob(sdss_path+'*.fits')
    print('sdss_names:',sdss_names)
    single_iterate=False
    double_iterate=False
    
elif file_setting=='all_SDSS':
    filenames= glob('*sdss*')
    sdss_names = glob(sdss_path+'*.fits')
    single_iterate=False
    double_iterate=False
    
else:
    print('\n\nno file_setting specificied\n\n')
    


def norm_spectrum(input_spec, norm_range, wave_range=True):
    if wave_range:
        print('\n\nwavelength-based norm range selected',norm_range)
        print('\n\n')
        out_flux= input_spec[1]/np.nanmean(input_spec[1][np.where((input_spec[0]>norm_range[0])&(input_spec[0]<norm_range[1]))])
        #out_flux= input_spec[1]/(np.nanmax(input_spec[1][np.where((input_spec[0]>norm_range[0])&(input_spec[0]<norm_range[1]))])-np.nanmin(input_spec[1][np.where((input_spec[0]>norm_range[0])&(input_spec[0]<norm_range[1]))]))
    else:
        print('\n\npixel-based norm range selected',norm_range)
        print('\n\n')
        pass
    out_spec=np.vstack([input_spec[0],out_flux])
    return out_spec


def convolve_spectrum(target_spec, header, kernel_type='gaussian', pix_width=3):
    #pix_width =3
    #pix_width =20
    fluxes= np.copy(target_spec[1])
    wavelengths = np.copy(target_spec[0])
    if kernel_type=='gaussian':
        see_sig = float(header['SEE_SIG']) #sigma value of gaussian fit to do the 
        see_kernel = conv.Gaussian1DKernel(see_sig, x_size = int(slit_width), mode = 'oversample')
        see_kernel.normalize()
        spec_conv = conv.convolve(fluxes, see_kernel)
    elif kernel_type== 'box':
        pix_kernel = conv.Box1DKernel(width = int(pix_width), mode = 'oversample')
        pix_kernel.normalize()
        spec_conv = conv.convolve(fluxes, pix_kernel)
    else:
        pass
    spec_out = np.vstack([wavelengths, spec_conv])
    return spec_out


def plot_spectrum(spec, filename, header, smooth=False, kernel_type='gaussian', norm=False, forced_title='', pix_width=3):
    title_string=filename
    if smooth:
        spec= convolve_spectrum(spec, header, kernel_type=kernel_type, pix_width=pix_width)
        title_string=title_string+ ' smoothed'
        #plt.title(filename+ ' smoothed')
    else:
        #plt.title(filename)
        pass
    if norm:
        #spec[1]= spec[1]/np.nanmean(spec[1])
        #spec[1]=spec[1]/np.nanmean(spec[1][1240:1280])
        #spec[1]=spec[1]/np.nanmean(spec[1][norm_range[0]:norm_range[1]])
        spec=norm_spectrum(spec, norm_range)
        #spec[1]=spec[1]/np.nanmax(spec[1])
        #spec[1]=spec[1]/np.nanmean(spec[1][1560:1590])
        #spec[1]=spec[1]/np.nanmean(spec[1][1020:1029])
        title_string= title_string+ ' normed'+str(norm_range)
        plt.ylabel('Flux (normed)')
    else:
        plt.ylabel('Flux (ergs/cm/cm/s/A 10**-16)')
        pass
    if plot_wavelength:
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.plot(spec[0], spec[1], label=filename)
    else:
        plt.xlabel('Pixel')
        plt.plot(spec[1], label=filename)
    print(title_string)
    plt.title(title_string)
    plt.xlim(np.nanmin(spec[0]), np.nanmax(spec[0]))
    plt.ylim(bottom=0)
    #plt.ylabel('Flux')
    #plt.title(filename)
    #plt.plot(spec[0], spec[1])
    #plt.show()
    return

def plot_dwavelength(spec):
    plt.ylabel(r'delta Wavelength ($\AA$)')
    plt.xlabel('Wavelength ($\AA$)')
    plt.title(filename)
    dlambda= spec[0][1:]-spec[0][:-1]
    plt.plot(spec[0][:-1], dlambda)
    plt.show()
    
    #plt.ylabel(r'Wavelength ($\AA$)')
    #plt.xlabel('pixel')
    #plt.title(filename)
    #plt.plot(spec[0])
    #plt.show()
    #return


def plot_SNR(spec, noise, filename):
    center_pixel = np.argmin(np.abs(spec[0]-test_wavelength))
    measured_std = np.std(spec[1][center_pixel-test_side:center_pixel+test_side])
    print("sigma in " + str(test_width) + " pixel range around " + str(test_wavelength)+ " angstroms", measured_std)
    sigma_range = noise[1][center_pixel-test_side:center_pixel+test_side]
    print("Predicted sigmas of " + str(test_width) +" pixel range around " + str(test_wavelength)+ " angstroms", "min:" + str(np.min(sigma_range)), "mean:" + str(np.mean(sigma_range)), "max:" + str(np.max(sigma_range)))
    print("Mean S/N:", np.mean(spec[1]/noise[1]))
    plt.xlabel('Noise')
    plt.title(filename)
    plt.hist(sigma_range)
    plt.show()
    
    if plot_wavelength:
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.vlines([spec[0][center_pixel-test_side],spec[0][center_pixel+test_side]], np.min(spec[1]/noise[1]), np.max(spec[1]/noise[1]))
        plt.plot(spec[0], spec[1]/noise[1], color = 'r')
    else:
        plt.xlabel('pixel')
        plt.vlines([center_pixel-test_side,center_pixel+test_side], np.min(spec[1]/noise[1]), np.max(spec[1]/noise[1]))
        plt.plot( spec[1]/noise[1], color = 'r')
    plt.ylabel('Signal/Noise')
    #plt.vlines([spec[0][center_pixel-test_side],spec[0][center_pixel+test_side]], np.min(spec[1]/noise[1]), np.max(spec[1]/noise[1]))
    plt.title(filename)
    #plt.plot(spec[0], spec[1]/noise[1], color = 'r')
    plt.show()
    return
    
def plot_sky(filename):
    hdu=fits.open(filename)
    sky=hdu[2].data
    if plot_wavelength:
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.plot(hdu[0].data, sky)
    else:
        plt.plot(sky)
        plt.xlabel('pixel')
    plt.title(filename+' sky')
    plt.show()
    return

def plot_diff_spec(spec1, spec2, filename1, filename2, header, smooth=False, kernel_type='gaussian', norm=False):
    plot_spectrum(spec1, filename1, header, norm=norm, smooth=smooth, kernel_type=kernel_type)
    plot_spectrum(spec2, filename2, header, norm=norm, smooth=smooth, kernel_type=kernel_type)
    title_string= filename1+' - ' + filename2
    if norm:
        #spec1[1]=spec1[1]/np.nanmean(spec1[1][norm_range[0]:norm_range[1]])
        #spec2[1]=spec2[1]/np.nanmean(spec2[1][norm_range[0]:norm_range[1]])
        spec1=norm_spectrum(spec1, norm_range)
        spec2=norm_spectrum(spec2, norm_range)
        title_string= title_string+ ' normed'+str(norm_range)
        plt.ylabel('Flux (normed)')
        #spec1[1]=spec1[1]/np.nanmax(spec1[1])
        #spec2[1]=spec2[1]/np.nanmax(spec2[1])
    else:
        pass
    diff_flux=spec1[1]-spec2[1]
    diff_spec=np.vstack([spec1[0], diff_flux])
    #plot_spectrum(diff_spec, filename1+' - ' + filename2, header, smooth=smooth, kernel_type=kernel_type)     
    plot_spectrum(diff_spec, title_string, header, smooth=smooth, kernel_type=kernel_type)
    plt.ylim(top=np.nanmax([spec1[1],spec2[1]])+0.5)
    plt.ylim(bottom=np.nanmin(diff_flux)-0.5)
    return


if file_setting== 'compare_SDSS':
    target_spec1, header1, target_noise1= spt.retrieve_spec(filename)
    target_spec1[0]=target_spec1[0]+wavelength_offset
    target_spec1= norm_spectrum(target_spec1, norm_range)
    for filename2 in sdss_names:
        target_spec2, header2, target_noise2= spt.retrieve_sdss_spec(filename2)
        target_spec2= spt.clean_spectrum(target_spec2, np.nanmin(target_spec1[0]), np.nanmax(target_spec1[0]), [])
        target_spec2=norm_spectrum(target_spec2, norm_range)
        plt.ylim(top=np.nanmax(np.hstack([target_spec2[1], target_spec1[1]]))+0.5)
        plot_spectrum(target_spec1, filename, header1, norm=True, smooth=True, kernel_type='box')
        plot_spectrum(target_spec2, filename2, header2, norm=True, smooth=True, kernel_type='box', pix_width=sdss_pix_width)
        #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=True, norm=True)
        #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=True)
        #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=False)
        plt.axhline(y=0, linestyle='--', color='k')
        plt.legend()
        plt.title(filename+ ' & '+ filename2.split('/')[-1])
        plt.show()
        
if file_setting=='all_SDSS':
    for filename1 in filenames:
        target_spec1, header1, target_noise1= spt.retrieve_sdss_spec(filename1)
        target_spec1=norm_spectrum(target_spec1, norm_range)
        for filename2 in sdss_names:
            target_spec2, header2, target_noise2= spt.retrieve_sdss_spec(filename2)
            target_spec2= spt.clean_spectrum(target_spec2, np.nanmin(target_spec1[0]), np.nanmax(target_spec1[0]), [])
            target_spec2=norm_spectrum(target_spec2, norm_range)
            plt.ylim(top=np.nanpercentile(np.hstack([target_spec2[1], target_spec1[1]]), 99.9)+0.5)
            plot_spectrum(target_spec1, filename1, header1, norm=False, smooth=True, kernel_type='box', pix_width= sdss_pix_width)
            plot_spectrum(target_spec2, filename2, header2, norm=True, smooth=True, kernel_type='box', pix_width=sdss_pix_width)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=True, norm=True)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=True)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=False)        plt.axhline(y=0, linestyle='--', color='k')
            plt.legend()
            plt.title(filename1+ ' & '+ filename2.split('/')[-1])
            plt.show()

if file_setting=='command':
        target_spec1, header1, target_noise1= spt.retrieve_spec(filename1)
        target_spec2, header2, target_noise2= spt.retrieve_spec(filename2)
        #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=True, norm=True)
        #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=True, norm=False)
        #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=True)
        plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=False)
        plt.axhline(y=0, linestyle='--', color='k')
        plt.legend()
        plt.show()

if single_iterate:
    for filename in filenames:
        target_spec, header, target_noise= spt.retrieve_spec(filename)
        #conv_spec= convolve_spectrum(target_spec, header)
        #plot_spectrum(target_spec, filename, header, smooth=True)
        #plot_spectrum(target_spec, filename, header, smooth=True, kernel_type='box', norm=True)
        #plot_spectrum(target_spec, filename, header)
        #plot_spectrum(target_spec, filename, header, norm=True)
        plot_spectrum(target_spec, filename, header, norm=True, smooth=True)
        #plot_spectrum(target_spec, filename, header, smooth=True)
        #plot_sky(filename)
        #plot_SNR(target_spec, target_noise, filename)
        #plot_dwavelength(target_spec)
        #plt.legend()
        #plt.show()
    plt.legend()
    plt.show()
else:
    pass


if double_iterate:
    for filename1 in filenames:
        target_spec1, header1, target_noise1= spt.retrieve_spec(filename1)
        for filename2 in filenames:
            target_spec2, header2, target_noise2= spt.retrieve_spec(filename2)
            plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=True, norm=True)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=True)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=False)
            plt.axhline(y=0, linestyle='--', color='k')
            plt.legend()
            plt.show()
else:
    pass

