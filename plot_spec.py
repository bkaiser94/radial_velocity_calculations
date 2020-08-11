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
from astropy.table import Table, Column
import scipy.interpolate as scinterp
import time
start = time.time()


#plt.rc('font', size =18)

#print start
import wdatmos
import spec_plot_tools as spt
import cal_params as cp

slit_width = 3.0 #arcseconds
pixel_scale = 0.3 #arcseconds per pixel_scale
slit_width = slit_width/pixel_scale #slit width in pixels
test_wavelength = 4686
test_width = 40
test_side = test_width/2

pix_width=3
sdss_pix_width = 10
#sdss_scale_factor=20.6 #BOSS scaling
sdss_scale_factor= 1.467 #SDSS spectrograph scaling
sdss_seeing=0.7 #arcsec seeing
sdss_see_sig=sdss_seeing/2.355/ pixel_scale

#wavelength_offset=60
#wavelength_offset=20
#wavelength_offset=15
#wavelength_offset=27
#wavelength_offset=35
#wavelength_offset=-30
wavelength_offset=0

#filenames = glob(sys.argv[1])
#filenames= glob('wctb*')
#filenames= glob('wctb*SDSS*')
sdss_path = '/Users/BenKaiser/Desktop/SDSS_speclib/'
#sdss_path= sdss_path+'G0_K5/'
#sdss_path = '/Users/BenKaiser/Desktop/SDSS_speclib/G0_K5/'

tell_filename='LBL_A30_s0_w005_R0060000_T.fits'
#tell_filename='LBL_A30_s0_w200_R0060000_T.fits'
#print(filenames)
plot_wavelength=True
plot_400m2_tell= False
#norm_range=[1240,1280]
#norm_range=[1560,1590]
#norm_range=[40,80]

#norm_range=[7042,7046]
#norm_range=[7490,7510] #outside telluric
#norm_range=[7470, 7530]
#norm_range=[7517,7556] #20190528
#norm_range=[7860.,8050.]
#norm_range=[8074,8140]
#norm_range=[8058,8231]
#norm_range=[5100,5400]
#norm_range=[6090,6240]
#norm_range=[6360,6420]
#norm_range=[6570,6620]
norm_range=[6640,6670]#20190530 400M1 norm range
#norm_range=[6630,6690]#wider double norm range
#norm_range=[5740,5850]
#norm_range=[5270,5560]
#norm_range=[3800, 4000]
####norm_range=np.array(norm_range)+wavelength_offset

#file_setting='all_avg'
#file_setting='command' #this is essentially the version for comparing 2 goodman spectra to each other
#file_setting='all_wctb'
#file_setting='all_fwctb'
file_setting= 'compare_SDSS'
#file_setting= 'compare_only_SDSS' #this should compare the spectra beginning with 'sdss' to other objects
#file_setting= 'all_SDSS'
#file_setting= 'two_arm'
#file_setting= 'all_super'
#file_setting= 'two_arm_compare_SDSS'
#file_setting='null' #option if you want to call this script in another script. It prevents anything from actually being executed.

single_iterate= False
double_iterate= False #file_settings change these in their little sections ahead if they should be changed


if file_setting=='all_avg':
    print(file_setting)
    filenames=glob('ravg_fwctb*fits')
    #filenames=glob('*ravg*2356*fits')
    #filenames=glob('ravg_fwctb*1150*fits')
    #filenames=glob('ravg_fwctb*2126*fits')
    #filenames=glob('ravg_fwctb*SDSS*n*fits')
    #filenames=glob('*avg_fwctb*DQpec*fits')
    #filenames=glob('ravg_wctb*fits')
    #filenames=glob('*avg_fwctb*eg274*')
    #filenames=glob('ravg_fwctb*Gaia*1644*')
    #filenames=glob('avg_fwctb*eg274*fits')
    #filenames=glob('ravg_fwctb*aia*1644*fits')
    #filenames=glob('avg_wctb*fits')
    #filenames=glob('avg_fwctb*fits')
    #filenames=glob('avg_fwctb*SDSSJ1252*')
    #filenames=glob('avg_*')
    single_iterate=True
    double_iterate=False
    #single_iterate=False
    #double_iterate=True

elif file_setting=='all_wctb':
    print(file_setting)
    filenames=glob('wctb*')
    #filenames=glob('wctb*Feige110_*')
    #filenames=glob('wctb*aia*2320*')
    single_iterate=True
    double_iterate=False
    

elif file_setting=='all_fwctb':
    print(file_setting)
    filenames=glob('fwctb*')
    #filenames=glob('fwctb*LTT*')
    #filenames=glob('fwctb*SDSS*n*')
    #filenames=glob('fwctb*2356*')
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
    #filename=glob('ravg_fwctb*')
    print('filename:', filename)
    #sdss_names = glob(sdss_path+'*.fits')
    #sdss_names = glob(sdss_path+'*M*.fits')
    #sdss_names = glob(sdss_path+'*K*.fits')
    sdss_names = glob(sdss_path+'SDSS*.fits')
    #sdss_names = glob(sdss_path+'*WDpec*.fits')
    #sdss_names = glob(sdss_path+'sdss*.fits')
    print('sdss_names:',sdss_names)
    single_iterate=True
    double_iterate=False
    
elif file_setting=='all_SDSS':
    #sdss_names= glob(sdss_path+'*sdss*')
    #sdss_names = glob(sdss_path+'*.fits')
    sdss_names = glob(sdss_path+'SDSSJ1636*.fits')
    filenames=sdss_names
    single_iterate=False
    double_iterate=False
    
elif file_setting== 'two_arm':
    #m1_names =glob('avg_fwctb*400m1*fits')
    #m2_names= glob('avg_fwctb*400m2*fits')
    m1_names =glob('ravg_fwctb*400m1*fits')
    m2_names= glob('ravg_fwctb*400m2*fits')
    #m1_names =glob('ravg_fwctb*1644*400m1*fits')
    #m2_names= glob('ravg_fwctb*1644*400m2*fits')
    #m1_names =glob('super_fwctb*400m1*fits')
    #m2_names= glob('super_fwctb*400m2*fits')
    #m1_names =glob('avg_wctb*400m1*fits')
    #m2_names= glob('avg_wctb*400m2*fits')
    filenames= zip(m1_names, m2_names)
    
elif file_setting=='all_super':
    print(file_setting)
    #filenames=glob('super_fwctb*1644*fits')
    filenames=glob('super_fwctb*fits')
    single_iterate=True
    double_iterate=False
    
elif file_setting =='two_arm_compare_SDSS':
    print(file_setting)
    filename1=sys.argv[1]
    filename2=sys.argv[2]
    sdss_names = glob(sdss_path+'*SDSS*.fits')
    #sdss_names = glob(sdss_path+'*M*.fits')
    single_iterate=False
    double_iterate=False
    
else:
    print('\n\nno file_setting specificied\n\n')
    


def norm_spectrum(input_spec, norm_range, wave_range=True, show_norm_range=True):
    if wave_range:
        print('\n\nwavelength-based norm range selected',norm_range)
        print('\n\n')
        out_flux= input_spec[1]/np.nanmean(input_spec[1][np.where((input_spec[0]>norm_range[0])&(input_spec[0]<norm_range[1]))])
        #out_flux= input_spec[1]/(np.nanmax(input_spec[1][np.where((input_spec[0]>norm_range[0])&(input_spec[0]<norm_range[1]))])-np.nanmin(input_spec[1][np.where((input_spec[0]>norm_range[0])&(input_spec[0]<norm_range[1]))]))
        if show_norm_range:
            plt.axvspan(norm_range[0], norm_range[1], alpha=0.1, color='r')
        else:
            pass
    else:
        print('\n\npixel-based norm range selected',norm_range)
        print('\n\n')
        pass
    
    out_spec=np.vstack([input_spec[0],out_flux])
    return out_spec


def convolve_spectrum(target_spec, header, kernel_type='gaussian', pix_width=pix_width, kernel_width=slit_width):
    #pix_width =3
    #pix_width =20
    fluxes= np.copy(target_spec[1])
    wavelengths = np.copy(target_spec[0])
    if kernel_type=='gaussian':
        #see_sig = float(header['SEE_SIG']) #sigma value of gaussian fit to do the 
        #see_sig= sdss_scale_factor*see_sig
        see_sig=pix_width
        #see_kernel = conv.Gaussian1DKernel(see_sig, x_size = int(slit_width), mode = 'oversample')
        if int(kernel_width)%2 ==0:
            kernel_width=kernel_width+1
        else:
            pass
        see_kernel = conv.Gaussian1DKernel(see_sig, x_size = int(kernel_width), mode = 'oversample')
        see_kernel.normalize()
        spec_conv = conv.convolve(fluxes, see_kernel)
    elif kernel_type=='sdss_match':
        see_sig = float(header['SEE_SIG']) #sigma value of gaussian fit to do the 
        conv_see_sig= np.sqrt(see_sig**2- sdss_see_sig**2)
        print('sdss_see_sig', sdss_see_sig)
        print('see_sig', see_sig)
        print('conv_see_sig', conv_see_sig)
        #see_sig= sdss_scale_factor*see_sig
        conv_see_sig=conv_see_sig*sdss_scale_factor
        
        #see_kernel = conv.Gaussian1DKernel(see_sig, mode = 'oversample')
        see_kernel=conv.Gaussian1DKernel(conv_see_sig, mode='oversample')
        see_kernel.normalize()
        spec_conv = conv.convolve(fluxes, see_kernel)
        pix_kernel = conv.Box1DKernel(width = int(sdss_scale_factor*pix_width), mode = 'oversample')
        pix_kernel.normalize()
        spec_conv = conv.convolve(spec_conv, pix_kernel)
    elif kernel_type== 'box':
        pix_kernel = conv.Box1DKernel(width = int(pix_width), mode = 'oversample')
        pix_kernel.normalize()
        spec_conv = conv.convolve(fluxes, pix_kernel)
    else:
        print('something wrong with convolution attempt')
        pass
    spec_out = np.vstack([wavelengths, spec_conv])
    return spec_out


def plot_spectrum(spec, filename, header, smooth=False, kernel_type='gaussian', norm=False, forced_title='', pix_width=pix_width, offset=0, color='None', kernel_width=slit_width, norm_range=norm_range):
    title_string=filename
    label_string= filename
    try:
        #label_string= label_string+str(header['airofavg'])[:4]
        pass
    except KeyError:
        pass
    if smooth:
        spec= convolve_spectrum(spec, header, kernel_type=kernel_type, pix_width=pix_width, kernel_width=kernel_width)
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
        #plt.ylabel('Flux (normed)')
        plt.ylabel(r'$f_{\lambda}$ (arbitrary units)')
    else:
        #plt.ylabel('Flux (ergs/cm/cm/s/A 10**-16)')
        try:
            plt.ylabel(header['units'])
        except TypeError:
            pass
        except KeyError:
            pass
    if plot_wavelength:
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.xlim(np.nanmin(spec[0]), np.nanmax(spec[0]))
        #plt.xlim(np.nanmin(spec[0]), np.nanmax(spec[0]))
        #color='b'
        #if 'eg274' in filename.lower():
            #color='r'
        #elif 'gd153' in filename.lower():
            #color='g'
        #elif 'wd2356' in filename.lower():
            #color='magenta'
        #elif 'sdssj1555' in filename.lower():
            #color= 'k'
        #else:
            #pass
        #plt.plot(spec[0], spec[1]+offset, label=filename, color=color)
        if color != 'None':
            plt.plot(spec[0], spec[1]+offset, label=label_string, color=color)
        else:
            plt.plot(spec[0], spec[1]+offset, label=label_string)
    else:
        plt.xlabel('Pixel')
        plt.plot(spec[1]+offset, label=label_string)
    print(title_string)
    plt.title(title_string)
    #plt.xlim(np.nanmin(spec[0]), np.nanmax(spec[0]))
    #plt.ylim(bottom=0)
    #plt.ylabel('Flux')
    #plt.title(filename)
    #plt.plot(spec[0], spec[1])
    #plt.show()
    return


def plot_telluric_spectrum(wave_range, pix_width=pix_width, smooth=False, kernel_type='gaussian', kernel_width= 300., tell_filename=tell_filename):
    tell_spec= spt.retrieve_telluric_model(tell_filename, wave_range)
    #plt.plot(tell_spec[0], tell_spec[1], label='Telluric absorption')
    plot_spectrum(tell_spec, 'Telluric Absorption', '', smooth=smooth, kernel_type=kernel_type, pix_width=pix_width, kernel_width=kernel_width)
    return


def plot_dwavelength(spec, filename, read_in=True):
    plt.ylabel(r'delta Wavelength ($\AA$)')
    plt.xlabel('Wavelength ($\AA$)')
    if read_in:
        hdu=fits.open(filename)
        dlambda= hdu[4].data
    
        plt.title(filename)
        plt.plot(spec[0], label=filename)
    else:
        dlambda= spec[0][1:]-spec[0][:-1]
        plt.plot(spec[0][:-1], dlambda, label=filename)
    #plt.plot(spec[0], dlambda, label=filename)
    #plt.show()
    
    #plt.ylabel(r'Wavelength ($\AA$)')
    #plt.xlabel('pixel')
    #plt.title(filename)
    #plt.plot(spec[0])
    #plt.show()
    return
    
def plot_pix_shifts(file_list):
    for filename in file_list:
        header= fits.getheader(filename)
        #print(filename, header['ROTATOR'],header['airofavg'], header['BMJD_TDB'])
        print(filename, header['ROTATOR'],header['CAM_ANG'], header['CAM_TARG'])
        #plt.errorbar(header['BMJD_TDB'], header['airofavg'], yerr=header['airofstd'], color='b')
        plt.errorbar(header['ROTATOR'], header['airofavg'], yerr=header['airofstd'], color='b', marker='o', markersize=4)
        #plt.scatter( header['BMJD_TDB'],header['ROTATOR'], color='b')
        #plt.errorbar(header['AIRMASS'], header['airofavg'], yerr=header['airofstd'], color='b')
    #plt.xlabel('BMDJ_TDB')
    plt.xlabel('ROTATOR')
    #plt.ylabel('ROTATOR')
    plt.ylabel('pixel offset of air lines')
    plt.show()
    return


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

def plot_SNR_from_file(filename):
    hdu=fits.open(filename)
    wavelengths=np.copy(hdu[0].data)
    scaled_noise= np.copy(hdu[3].data)
    plt.plot(wavelengths, 1/scaled_noise, label=filename)
    plt.title(filename)
    plt.ylabel('S/N')
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.xlim(np.nanmin(wavelengths), np.nanmax(wavelengths))
    return

def get_median_dlambda(input_spec):
    
    return np.nanmedian(input_spec[0]-np.roll(input_spec[0],1))

    
def plot_sky(filename, offset=0, line_labels=True, convolve=False):
    hdu=fits.open(filename)
    wavelengths=np.copy(hdu[0].data)
    sky=np.copy(hdu[2].data)
    sky_spec=np.vstack([wavelengths, sky])
    if convolve:
        sky_spec= convolve_spectrum(sky_spec, 'dummy_header', kernel_type='box', pix_width=pix_width)
    else:
        pass
    sky= sky_spec[1]
    sky=sky/np.nanmax(sky)
    sky=sky+offset
    #color='b'
    #if 'eg274' in filename.lower():
        #color='r'
    #elif 'gd153' in filename.lower():
        #color='g'
    #elif 'wd2356' in filename.lower():
        #color='magenta'
    #elif 'sdssj1555' in filename.lower():
        #color= 'k'
    #else:
        #pass
    airline_array= Table.read(cp.line_list_dir+cp.airline_name, format='ascii.tab')
    #print(airline_array)
    use_array=np.int_(airline_array['use'])
    #good_airlines= np.copy(airline_array[np.where(use_array==1)])
    good_airlines= np.copy(airline_array)
    air_waves = np.float_(good_airlines['User'])
    #air_names= good_airlines['Name']+good_airlines['Name2']
    if line_labels:
        for air_wave, name, name2 in zip(air_waves, good_airlines['Name'], good_airlines['Name2']):
            #print(name+name2, type(name))
            air_name=name+name2
            plt.axvline(x=air_wave, linestyle='--', color=cp.airline_color[air_name[:2]])
            plt.text(air_wave, np.nanmax(sky), air_name, color=cp.airline_color[air_name[:2]], rotation=90)
    else:
        pass
    print("need to put dlambda into this part again since we're about to move around wavelengths in the future.")
    if plot_wavelength:
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.plot(hdu[0].data, sky, label=filename+' sky')
        #plt.plot(hdu[0].data, sky, label=filename, color=color)
        #plt.plot(hdu[0].data, sky, label=hdu[0].header['airmass'])
    else:
        plt.plot(sky)
        plt.xlabel('pixel')
    plt.title(filename+' sky')
    #plt.show()
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
    #plot_spectrum(diff_spec, title_string, header, smooth=smooth, kernel_type=kernel_type)
    plt.ylim(top=np.nanmax([spec1[1],spec2[1]])+0.5)
    plt.ylim(bottom=np.nanmin(diff_flux)-0.5)
    return

if __name__ == '__main__':

    if file_setting== 'compare_SDSS':
        target_spec1, header1, target_noise1= spt.retrieve_spec(filename)
        #target_spec1, header1, target_noise1= spt.retrieve_sdss_spec(filename)
        target_spec1[0]=target_spec1[0]+wavelength_offset
        target_spec1= norm_spectrum(target_spec1, norm_range)
        for filename2 in sdss_names:
            target_spec2, header2, target_noise2= spt.retrieve_sdss_spec(filename2)
            #target_spec2= spt.clean_spectrum(target_spec2, np.nanmin(target_spec1[0]), np.nanmax(target_spec1[0]), [])
            target_spec2=norm_spectrum(target_spec2, norm_range)
            #plt.ylim(top=np.nanmax(np.hstack([target_spec2[1], target_spec1[1]]))+0.5)
            plt.ylim(top=np.nanpercentile(np.hstack([target_spec2[1], target_spec1[1]]), 99.9)+0.5)
            plot_spectrum(target_spec1, filename, header1, norm=True, smooth=True, kernel_type='box')
            #plot_spectrum(target_spec1, filename, header1, norm=True, smooth=True, kernel_type='gaussian')
            #plot_spectrum(target_spec1, filename, header1, norm=False, smooth=True, kernel_type='box')
            #plot_spectrum(target_spec1, filename, header1, norm=True, smooth=True, kernel_type='box', pix_width=sdss_pix_width)
            #plot_spectrum(target_spec2, filename2.split('/')[-1], header2, norm=True, smooth=True, kernel_type='box', pix_width=sdss_pix_width)
            plot_spectrum(target_spec2, filename2.split('/')[-1], header1, norm=True, smooth=True, kernel_type='sdss_match', pix_width=pix_width)
            #plot_spectrum(target_spec2, filename2, header2, norm=False, smooth=True, kernel_type='box', pix_width=sdss_pix_width)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=True, norm=True)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=True)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=False)
            plt.axhline(y=0, linestyle='--', color='k')
            #plt.legend()
            plt.title(filename+ ' & '+ filename2.split('/')[-1])
            #plt.show()
            spt.show_plot(line_id='alkali')
            
    if file_setting=='all_SDSS':
        for filename1 in filenames:
            target_spec1, header1, target_noise1= spt.retrieve_sdss_spec(filename1)
            target_spec1=norm_spectrum(target_spec1, norm_range)
            for filename2 in sdss_names:
                target_spec2, header2, target_noise2= spt.retrieve_sdss_spec(filename2)
                target_spec2= spt.clean_spectrum(target_spec2, np.nanmin(target_spec1[0]), np.nanmax(target_spec1[0]), [])
                target_spec2=norm_spectrum(target_spec2, norm_range)
                plt.ylim(top=np.nanpercentile(np.hstack([target_spec2[1], target_spec1[1]]), 99.9)+0.5)
                plot_spectrum(target_spec1, filename1.split('/')[-1], header1, norm=True, smooth=True, kernel_type='box', pix_width= sdss_pix_width)
                plot_spectrum(target_spec2, filename2.split('/')[-1], header2, norm=True, smooth=True, kernel_type='box', pix_width=sdss_pix_width)
                tell_filename='LBL_A30_s0_w005_R0060000_T.fits'
                #tell_filename='LBL_A30_s0_w200_R0060000_T.fits'
                plot_telluric_spectrum([3700, 9000], smooth=True, pix_width=30, tell_filename=tell_filename)
                tell_filename='LBL_A30_s0_w200_R0060000_T.fits'
                plot_telluric_spectrum([3700, 9000], smooth=True, pix_width=30,tell_filename=tell_filename)
                #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=True, norm=True)
                #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=True)
                #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=False)        plt.axhline(y=0, linestyle='--', color='k')
                #plt.legend()
                plt.title(filename1+ ' & '+ filename2.split('/')[-1])
                #plt.show()
                spt.show_plot(line_id='alkali', convert_to_air=False)
        #plt.legend()
        #plt.show()

    if file_setting=='command':
            target_spec1, header1, target_noise1= spt.retrieve_spec(filename1)
            target_spec1[0]=target_spec1[0]+wavelength_offset
            target_spec2, header2, target_noise2= spt.retrieve_spec(filename2)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=True, norm=True, kernel_type='box')
            plot_spectrum(target_spec1, filename1, header1, norm=True, smooth=False, kernel_type='box', pix_width=pix_width)
            #plot_spectrum(target_spec1, filename1, header1, norm=True, smooth=True, kernel_type='box')
            #plot_spectrum(target_spec1, filename1, header1, norm=True, smooth=True, kernel_type='box', pix_width=pix_width*2)
            plot_spectrum(target_spec2, filename2, header2, norm=True, smooth=False, kernel_type='box', pix_width=pix_width)
            #plot_spectrum(target_spec1, filename1, header1, norm=True, smooth=True, kernel_type='box', pix_width=10)
            #plot_spectrum(target_spec2, filename2, header2, norm=True, smooth=True, kernel_type='box', pix_width=10)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=True, norm=False)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=True)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=False)
            plt.axhline(y=0, linestyle='--', color='k')
            #plt.legend()
            #plt.show()
            spt.show_plot()
            
    if file_setting=='two_arm':
        #sdss_names = glob(sdss_path+'*1555*.fits')
        #sdss_names = glob(sdss_path+'*1159*.fits')
        #sdss_spec, sdssheader, sdss_noise= spt.retrieve_sdss_spec(sdss_names[0])
        #plot_spectrum(sdss_spec, sdss_names[0], sdssheader, norm=False, smooth=True, kernel_type='box', pix_width=sdss_pix_width)
        #plt.xlim(3700,9000)
        #spt.show_plot() 
        for m1_name, m2_name in filenames:
            target_spec1, header1, target_noise1= spt.retrieve_spec(m1_name)
            target_spec2, header2, target_noise2= spt.retrieve_spec(m2_name)
            
            ##sdss_spec= spt.clean_spectrum(target_spec2, np.nanmin(target_spec1[0]), np.nanmax(target_spec1[0]), [])
            #target_spec1[0]=target_spec1[0]+wavelength_offset
            #plot_sky(m1_name, offset=0)
            #plot_sky(m2_name, offset=0)
            
            plot_spectrum(target_spec1, m1_name, header1, norm=True, smooth=False, kernel_type='box')
            plot_spectrum(target_spec2, m2_name, header2, norm=True, smooth=False, kernel_type='box')
            
            #plt.scatter(target_spec1[0],target_spec1[1], color='b')
            #plt.scatter(target_spec2[0], target_spec2[1], color='r')
            #plt.show()
            #dl1=target_spec1[0]-np.roll(target_spec1[0],1)
            #dl2=target_spec2[0]-np.roll(target_spec2[0],1)
            #plt.scatter(target_spec1[0][1:],dl1[1:], color='b', label='400M1')
            #plt.scatter(target_spec2[0][1:],dl2[1:] , color='r', label='400M2')
            #plt.ylabel(r'$\Delta \lambda (\AA)$')
            #plt.xlabel(r'Wavelength ($\AA$)')
            #plt.title(header1['date-obs']+'\t'+header2['date-obs'])
            #plt.legend()
            #plt.show()
            
            
            #target_spec2=norm_spectrum(target_spec2, norm_range)
            #target_spec1=norm_spectrum(target_spec1, norm_range)
            #plot_spectrum(target_spec1, m1_name, header1, norm=True, smooth=True, kernel_type='box')
            #plot_spectrum(target_spec2, m2_name, header2, norm=True, smooth=True, kernel_type='box')
            
            #target_spec2=norm_spectrum(target_spec2, norm_range)
            #target_spec1=norm_spectrum(target_spec1, norm_range)
            #plot_spectrum(target_spec1, m1_name, header1, norm=True, smooth=False, kernel_type='gaussian')
            #plot_spectrum(target_spec2, m2_name, header2, norm=True, smooth=False, kernel_type='gaussian')
            
            #plot_dwavelength(target_spec1, m1_name)
            #plot_dwavelength(target_spec2,m2_name)
            
            #plot_spectrum(target_spec1, m1_name, header1, norm=False, smooth=True, kernel_type='box')
            #plot_spectrum(target_spec2, m2_name, header2, norm=False, smooth=True, kernel_type='box')
            
            #try:
                #hdu1= fits.open(m1_name)
                #dlambda1=  hdu1[4].data
                #hdu2= fits.open(m2_name)
                #dlambda2= hdu2[4].data
                #nu_spec1= spt.flambda_to_fnu(target_spec1, dlambda1)
                #nu_spec2=spt.flambda_to_fnu(target_spec2, dlambda2)
                #plot_spectrum(nu_spec2, 'fnu2', header2, smooth=False, norm=False, kernel_type='box')
                #plot_spectrum(nu_spec1, 'fnu1', header1,smooth=False, norm=False, kernel_type='box')
                #plt.ylabel(r'$f_{\nu}$ $10^{-28} erg cm^{-2} s}^{-1} Hz^{-1}$')
                #plt.title(m1_name)
            #except IndexError as error:
                #print(error)
            
            #plot_spectrum(sdss_spec, sdss_names[0], sdssheader, norm=False, smooth=True, kernel_type='box', pix_width=sdss_pix_width)
            
            #plt.ylim(top=np.percentile(np.hstack([target_spec1[1], target_spec2[1]]),99.9)*1.1)
            plt.xlim(3700,9000)
            spt.show_plot(line_id='alkali', convert_to_air=True)
        spt.show_plot()
            
    if file_setting=='two_arm_compare_SDSS':
        target_spec1, header1, target_noise1= spt.retrieve_spec(filename1)
        target_spec2, header2, target_noise2= spt.retrieve_spec(filename2)
        
        
        target_spec1[0]=spt.air_to_vac(target_spec1[0])
        target_spec2[0]=spt.air_to_vac(target_spec2[0])
        
        #target_spec1= spt.flambda_to_fnu(target_spec1)
        #target_spec2=spt.flambda_to_fnu(target_spec2)
        
        target_spec1=norm_spectrum(target_spec1, norm_range)
        target_spec2= norm_spectrum(target_spec2, norm_range)
        
        for sdss_filename in sdss_names:
            
            sdss_spec, sdss_header, sdss_noise2= spt.retrieve_sdss_spec(sdss_filename, wave_medium='vac')
            
            #sdss_spec[0]=spt.vac_to_air(sdss_spec[0])
            
            sdss_spec= spt.clean_spectrum(sdss_spec, 3700, 9000, [])
            
            
            #sdss_spec=spt.flambda_to_fnu(sdss_spec)
            
            sdss_spec= norm_spectrum(sdss_spec, norm_range)
            
            
            print('header1 see_sig', header1['SEE_SIG'])
            print('header2 see_sig', header2['SEE_SIG'])
            print('median dlambda spec1:', get_median_dlambda(target_spec1))
            print('median dlambda spec2:', get_median_dlambda(target_spec2))
            print('median dlambda sdss spec:', get_median_dlambda(sdss_spec))
            print('scaling factor:', get_median_dlambda(target_spec1)/get_median_dlambda(sdss_spec))
            print('scaling factor:', get_median_dlambda(target_spec2)/get_median_dlambda(sdss_spec))
            plt.ylim(top=np.percentile(np.hstack([target_spec2[1], target_spec1[1], sdss_spec[1]]), 99.9)+0.5)
            
            plot_spectrum(sdss_spec, sdss_filename.split('/')[-1], sdss_header, norm=True, smooth=True, kernel_type='box', pix_width=sdss_pix_width, color='r')
            #plot_spectrum(sdss_spec, sdss_filename.split('/')[-1], header1, norm=True, smooth=True, kernel_type='sdss_match', pix_width=pix_width, color='r')
            plot_spectrum(target_spec1, filename1, header1, norm=True, smooth=True, kernel_type='box', pix_width=pix_width, color='g')
            plot_spectrum(target_spec2, filename2, header2, norm=True, smooth=True, kernel_type='box', pix_width=pix_width, color='b')
            #plot_dwavelength(target_spec1, filename1, read_in=False)
            #plot_dwavelength(target_spec2, filename2, read_in=False)
            #plot_dwavelength(sdss_spec, sdss_filename.split('/')[-1], read_in=False)
            plt.xlim(3700,9000)
            #plt.ylabel(r'$F_{\nu}$ (normalized)')
            spt.show_plot(line_id='alkali', convert_to_air=True)
            
    if single_iterate:
        counter=0
        for filename in filenames:
            target_spec, header, target_noise= spt.retrieve_spec(filename)
            hdu= fits.open(filename)
            #dlambda= hdu[4].data
            #target_spec[0]=target_spec[0]+wavelength_offset
            #print(filename, 'mean: ', np.nanmean(target_spec[1]))
            #nu_spec= spt.flambda_to_fnu(target_spec, dlambda)
            #conv_spec= convolve_spectrum(target_spec, header)
            #plot_spectrum(target_spec, filename, header, smooth=True, norm=False, kernel_type='box')
            #plt.scatter(header['BMJD_TDB'], np.sum(target_spec[1]*dlambda))
            #plt.errorbar(target_spec[0], target_spec[1], yerr=target_noise[1], label=filename, marker='o')
            #plot_spectrum(nu_spec, 'fnu', header, smooth=True, norm=False, kernel_type='box')
            plot_spectrum(target_spec, filename, header, smooth=True, norm=True, pix_width=header['see_sig'], kernel_type='gaussian')
            #plot_spectrum(target_spec, filename, header, smooth=True, norm=False,  kernel_type='gaussian',pix_width=header['SEE_SIG'])
            #target_spec[1]=header['airmass']
            #plot_spectrum(target_spec, filename, header, norm=False, smooth=True, kernel_type='box', pix_width=10)
            #plot_spectrum(target_spec, str(header['airmass']), header, norm=False, smooth=True, kernel_type='box', pix_width=10)
            #plot_spectrum(target_spec, filename, header, norm=True, smooth=False,offset=counter*-0.2, kernel_type='box')
            #plot_spectrum(target_spec, filename, header, norm=True, smooth=True, kernel_type='box', offset=counter)
            #plot_spectrum(nu_spec, filename, header, norm=False, smooth=True, kernel_type='box')
            #plt.plot(target_spec[0], dlambda,  label=filename, marker='o', markersize=10-counter)
            #plot_spectrum(target_spec, filename, header, smooth=True, kernel_type='gaussian', norm=True)
            #plot_sky(filename, offset=0, line_labels=False, convolve=False)
            #if header['airmass']<1.5:
                #plot_sky(filename, offset=0)
            #else:
                #pass
            #plot_SNR(target_spec, target_noise, filename)
            #plot_SNR_from_file(filename)
            #plot_dwavelength(target_spec, filename)
            #plot_telluric_spectrum([3700, 9000], smooth=True, pix_width=30)
            #plot_telluric_spectrum([3700,9000], smooth=True, pix_width=30, tell_filename='LBL_A30_s0_w200_R0060000_T.fits')
            #spt.show_plot(show_telluric=False, show_legend=False)
            #spt.show_plot(show_legend=True, line_id='alkali', convert_to_air=True)
            #plt.legend()
            #plt.show()
            
            counter+=1
            #try:
                #plt.title(header['airoftyp'])
            #except KeyError:
                #pass
        #plt.xlim(3700,9000)
        #plt.ylabel(r'$f_{\nu}$ (arbitrary units)')
        #plt.legend(loc='best')
        #plt.show()
        #plot_telluric_spectrum([3700, 9000], smooth=True, pix_width=30)
        plot_telluric_spectrum([3700, 9000], smooth=True, pix_width=30)
        plot_telluric_spectrum([3700,9000], smooth=True, pix_width=30, tell_filename='LBL_A30_s0_w200_R0060000_T.fits')
        spt.show_plot(show_legend=True, line_id='alkali', convert_to_air=True)
        #plt.ylabel('Integrated Flux (10^-16 erg/cm^2/s)')
        #plt.xlabel('BMJD_TDB')
        #plt.show()
        #spt.show_plot(show_legend=True, show_telluric=True)
        #spt.show_plot(show_telluric=True, )
        #plt.legend()
        #plt.show()
        plot_pix_shifts(filenames)
    else:
        pass


    if double_iterate:
        for filename1 in filenames:
            target_spec1, header1, target_noise1= spt.retrieve_spec(filename1)
            for filename2 in filenames:
                target_spec2, header2, target_noise2= spt.retrieve_spec(filename2)
                plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=True, kernel_type='box')
                #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=True)
                #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=False)
                #plot_spectrum(target_spec1, filename1, header1, norm=False, smooth=True, kernel_type='box')
                #plot_spectrum(target_spec2, filename2, header2, norm=False, smooth=True, kernel_type='box')
                plt.axhline(y=0, linestyle='--', color='k')
                #plt.legend()
                #plt.show()
                spt.show_plot()
    else:
        pass

