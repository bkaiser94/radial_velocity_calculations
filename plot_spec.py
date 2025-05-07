"""
Created 2019-01-29 by Ben Kaiser (UNC-Chapel Hill)

@author: Ben Kaiser

This should just be able to plot a given spectrum and various things about it.

"""
from __future__ import print_function


#import matplotlib
#matplotlib.use('pdf')


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

plt.rc('lines',linewidth=0.5)
#plt.rc('font', size =18)

#print start
#import wdatmos
import spec_plot_tools as spt
import cal_params as cp

slit_width = 3.2 #arcseconds
pixel_scale = 0.3 #arcseconds per pixel_scale
slit_width = slit_width/pixel_scale #slit width in pixels
test_wavelength = 4875
test_width = 40
test_side = test_width/2

pix_width=10
sdss_pix_width = 10
boss_scale_factor=20.6 #BOSS scaling or really probably the drizzled combined spectra from Kesseli et al. 2017 because they up-resolved their class-integrated spectra
sdss_scale_factor= 1.467 #SDSS spectrograph scaling
#sdss_scale_factor=boss_scale_factor#BOSS scaling

sdss_seeing=0.7 #arcsec seeing
sdss_see_sig=sdss_seeing/2.355/ pixel_scale

sdss_flux_scale_factor=0.5 #factor by which to multiply normalized SDSS spectrum before subtracting in compare_SDSS

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

#output_dir='/Users/BenKaiser/Desktop/white_dwarfs_posing_as_Kstars_PDFs/'
output_dir='/Users/BenKaiser/Desktop/white_dwarfs_posing_as_Kstars_PDFs_fluxes/'

tell_filename='LBL_A30_s0_w005_R0060000_T.fits'
#tell_filename='LBL_A15_s0_w015_R0060000_T.fits'
#tell_filename='LBL_A30_s0_w200_R0060000_T.fits'

#print(filenames)
plot_wavelength=True
plot_400m2_tell= False
default_norm_show=True
save_plots=False

#########################
bb_teff=9000. #K


#mwdd_filename='WDJ0254p4255_DXP_Kilic2025_fnu_MWDD_spec.csv'
#mwdd_filename='GJ3866_DC_Heatm_Teff5000_Bergeron2001_fnu_MWDD_spec.csv'
#mwdd_filename='WD1208p076_DA_Teff5150_subsavage2009_fnu_MWDD_spec.csv'
#mwdd_filename='WDJ1753m8405_DCsus_Teff5020K_OBrien2024_flambda_MWDD_spec.csv'
mwdd_filename='WDJ1838m5357_DAsus_Teff5150_OBrien2024_flambda_MWDD_spec.csv'
#mwdd_filename='WDJ1828m6537_DAsus_Teff5050_OBrien2024_flambda_MWDD_spec.csv'
#mwdd_filename='WDJ1018m3438_DAsus_Teff5090_OBrien2024_flambda_MWDD_spec.csv'
#mwdd_filename='WDJ2236m5548_DZsus_Teff5132_OBrien2024_flambda_MWDD_spec.csv'
#mwdd_filename='WDJ1604p2310_DAsus_Teff5100_Kilic2025_fnu_MWDD_spec.csv'
#mwdd_filename='LSPMJ1951p4209_DAsus_Teff5100_Tremblay2020_fnu_MWDD_spec.csv'
#mwdd_filename='Lawd34_DC_Teff7096K_OBrien2024_flambda_MWDD_spec.csv'
#mwdd_filename='WDJ1018m3438_DAsus_Teff5090_OBrien2024_flambda_MWDD_spec.csv'
#mwdd_filename='GJ1221_DXP_Giammichele2012_fnu_MWDD_spec.csv'
#mwdd_filename='GJ1228_DXP_Bergeron2011_fnu_MWDD_spec.csv'

mwdd_list=[
    'WDJ1753m8405_DCsus_Teff5020K_OBrien2024_flambda_MWDD_spec.csv',
    'WDJ1838m5357_DAsus_Teff5150_OBrien2024_flambda_MWDD_spec.csv',
    'WDJ1828m6537_DAsus_Teff5050_OBrien2024_flambda_MWDD_spec.csv',
    'WDJ1018m3438_DAsus_Teff5090_OBrien2024_flambda_MWDD_spec.csv'
    #'WDJ2236m5548_DZsus_Teff5132_OBrien2024_flambda_MWDD_spec.csv'
    
    ]

##########################

#norm_range=[1240,1280]
#norm_range=[1560,1590]
#norm_range=[40,80]

#norm_range=[7042,7046]
#norm_range=[7490,7510] #outside telluric
#norm_range=[7470, 7530]
#norm_range=[7440, 7550]
#norm_range=[7517,7556] #20190528
#norm_range=[7860.,8050.]
#norm_range=[8074,8140]
#norm_range=[8058,8231]
#norm_range=[5100,5400]
#norm_range=[4420,4500]
#norm_range=[4800.,5000.]
#norm_range=[6090,6240]
norm_range=[6360,6420]
#norm_range=[6570,6620]
#norm_range=[6640,6670]#20190530 400M1 norm range
#norm_range=[6630,6690]#wider double norm range
#norm_range=[6773, 6817]#red side of Li I line short of telluric
#norm_range=[6630,6690]
#norm_range=[5740,5850]
#norm_range=[5310,5490]
#norm_range=[4600, 4700]
#norm_range=[8640,8790]
#norm_range=[4975,5050]
#norm_range=[6745, 6815] #DQpec normalization range
####norm_range=np.array(norm_range)+wavelength_offset

file_setting='all_avg'
#file_setting='command' #this is essentially the version for comparing 2 goodman spectra to each other
#file_setting='all_wctb'
#file_setting='all_fwctb'
#file_setting= 'compare_SDSS'
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
    #filenames=glob('ravg_wctb*2320*bkg*.fits')
    filenames=glob('ravg_wctb*2320*off*6_*.fits')
    #filenames=glob('fravg_wctb*2320*.fits')
    #filenames=glob('*ravg_wctb*fits')
    #filenames=glob('ravg_fwctb*2124*fits')
    #filenames=glob('ravg_fwctb*0212*fits')
    #filenames=glob('ravg_fwctb*WDJ1426*fits')
    #filenames=glob('ravg_fwctb*1824*other2_*fits')
    #filenames=glob('ravg_fwctb*J0850*fits')
    #filenames=glob('ravg_fwctb*other*fits')
    #filenames=glob('ravg_fwctb*20220325*fits')
    #filenames=glob('ravg_zfwctb*fits')
    #filenames=glob('*ravg_fwctb*J0400*fits')
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
    
    filenames= sorted(filenames)
    single_iterate=True
    double_iterate=False
    #single_iterate=False
    #double_iterate=True

elif file_setting=='all_wctb':
    print(file_setting)
    filenames=glob('wctb*2320*6637_*')
    #filenames=glob('wctb*J0850*')
    #filenames=glob('wctb*4414other_*')
    #filenames=glob('wctb*1824*1213_*')
    #filenames=glob('wctb*Feige110_*')
    #filenames=glob('wctb*aia*2320*')
    
    filenames= sorted(filenames)
    single_iterate=True
    double_iterate=False
    

elif file_setting=='all_fwctb':
    print(file_setting)
    filenames=glob('fwctb*')
    #filenames=glob('fwctb*0850p195*')
    #filenames=glob('fwctb*WDJ2124*')
    #filenames=glob('fwctb*SDSSJ1312*')
    #filenames=glob('fwctb*1430*')
    #filenames=glob('fwctb*WISE*')
    #filenames=glob('fwctb*Feige*')
    #filenames=glob('fwctb*SDSS*n*')
    #filenames=glob('fwctb*J0834*')
    
    filenames= sorted(filenames)
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
    #sdss_names = glob(sdss_path+'*Dwarf*.fits')
    #sdss_names = glob(sdss_path+'G0_K5/*K*.fits')
    #sdss_names = glob(sdss_path+'*1651*.fits')
    #sdss_names = glob(sdss_path+'*M*.fits')
    #sdss_names = glob(sdss_path+'*K*.fits')
    #sdss_names = glob(sdss_path+'SDSS*1636*.fits')
    #sdss_names = glob(sdss_path+'LHS*.fits')
    #sdss_names = glob(sdss_path+'*DQpec*.fits')
    #sdss_names = glob(sdss_path+'*sdss*.fits')
    #sdss_names=glob(sdss_path+'*J0804*.fits')
    #sdss_names=glob(sdss_path+'*Dox*.fits')
    #sdss_names=glob(sdss_path+'Kstars_mistaken_for_WDs/*.fits')
    #sdss_names=glob(sdss_path+'Kstars_mistaken_for_WDs/*SDSSJ0738p4114*.fits')
    #sdss_names=glob(sdss_path+'Kstars_mistaken_for_WDs/*SDSSJ1312*.fits')
    #sdss_names=glob(sdss_path+'Kstars_mistaken_for_WDs/*SDSSJ2348*.fits')
    sdss_names=glob(sdss_path+'WD_pec/*.fits')
    #sdss_names=glob(sdss_path+'DC/*.fits')

    sdss_names= sorted(sdss_names)
    #print('sdss_names:',sdss_names)
    single_iterate=False
    double_iterate=False
    
elif file_setting=='all_SDSS':
    #sdss_names= glob(sdss_path+'*sdss*')
    #sdss_names = glob(sdss_path+'*.fits')
    #sdss_names = glob(sdss_path+'*DQpec*.fits')
    #sdss_names = glob(sdss_path+'SDSSJ1636*.fits')
    #sdss_names = glob(sdss_path+'SDSSJ1330*.fits')
    sdss_names=sorted(glob(sdss_path+'G0_K5/*K*.fits'))
    #sdss_names=sorted(glob(sdss_path+'*Dwarf*.fits'))
    filenames = sorted(glob(sdss_path+'Kstars_mistaken_for_WDs/*SDSSJ2348*.fits'))
    #filenames=sdss_names
    single_iterate=False
    double_iterate=False
    
elif file_setting== 'two_arm':
    #m1_names =glob('avg_fwctb*400m1*fits')
    #m2_names= glob('avg_fwctb*400m2*fits')
    #m1_names =glob('ravg_fwctb*WISE*400m1*fits')
    #m2_names= glob('ravg_fwctb*WISE*400m2*fits')
    m1_names =glob('ravg_fwctb*WDJ1948*m1*fits')
    m2_names= glob('ravg_fwctb*WDJ1948*m2*fits')
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
    #sdss_names = glob(sdss_path+'*Dwarf*.fits')
    #sdss_names = glob(sdss_path+'*sdss*.fits')
    #sdss_names = glob(sdss_path+'*SDSS*.fits')
    #sdss_names = glob(sdss_path+'*Dox*.fits')
    sdss_names = glob(sdss_path+'*WD_pec/*.fits')
    #sdss_names = glob(sdss_path+'*DQpec*.fits')
    #sdss_names = glob(sdss_path+'*M*.fits')
    #sdss_names = glob(sdss_path+'*1636*.fits')
    #sdss_names = glob(sdss_path+'G0_K5/*G*.fits')
    #sdss_names = glob(sdss_path+'DC/*.fits')
    
    sdss_names= sorted(sdss_names)
    single_iterate=False
    double_iterate=False
    
else:
    print('\n\nno file_setting specificied\n\n')
    


def norm_spectrum(input_spec, norm_range, wave_range=True, show_norm_range=default_norm_show):
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


def plot_spectrum(spec, filename, header, smooth=False, kernel_type='gaussian', norm=False, forced_title='', pix_width=pix_width, offset=0, color='None', kernel_width=slit_width, norm_range=norm_range,alpha=1):
    title_string=filename
    label_string= filename
    try:
        #label_string= label_string+str(header['airofavg'])[:4]
        #label_string= label_string+','+str(header['BMJD_TDB'])[:7]+','+str(header['see_FWHM'])[:3]
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
        plt.xlabel(r'Wavelength ($\mathrm{\AA}$)')
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
            plt.plot(spec[0], spec[1]+offset, label=label_string, color=color,alpha=alpha)
        else:
            plt.plot(spec[0], spec[1]+offset, label=label_string,alpha=alpha)
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


def plot_telluric_spectrum(wave_range, pix_width=pix_width, smooth=False, kernel_type='gaussian', kernel_width= 300., tell_filename=tell_filename, color= 'None',show_filename=False):
    label= 'Telluric Absorption'
    if show_filename:
        label=label+', '+tell_filename
    else:
        pass
    tell_spec= spt.retrieve_telluric_model(tell_filename, wave_range)
    #plt.plot(tell_spec[0], tell_spec[1], label='Telluric absorption')
    if color == 'None':
        plot_spectrum(tell_spec,label, '', smooth=smooth, kernel_type=kernel_type, pix_width=pix_width, kernel_width=kernel_width)
    else:
        plot_spectrum(tell_spec,label, '', smooth=smooth, kernel_type=kernel_type, pix_width=pix_width, kernel_width=kernel_width, color=color)
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
        #plt.scatter( header['AIRMASS'],header['ROTATOR'], color='b')
        #plt.errorbar(header['AIRMASS'], header['airofavg'], yerr=header['airofstd'], color='b')
    #plt.xlabel('BMDJ_TDB')
    plt.xlabel('ROTATOR')
    #plt.ylabel('ROTATOR')
    plt.ylabel('pixel offset of air lines')
    plt.show()
    return


def plot_rot_time(file_list):
    for filename in file_list:
        header= fits.getheader(filename)
        #print(filename, header['ROTATOR'],header['airofavg'], header['BMJD_TDB'])
        #print(filename, header['ROTATOR'],header['CAM_ANG'], header['CAM_TARG'])
        #plt.errorbar(header['BMJD_TDB'], header['airofavg'], yerr=header['airofstd'], color='b')
        #plt.errorbar(header['ROTATOR'], header['airofavg'], yerr=header['airofstd'], color='b', marker='o', markersize=4)
        plt.scatter( header['BMJD_TDB'],header['ROTATOR'], color='b')
        #plt.scatter( header['AIRMASS'],header['ROTATOR'], color='b')
        #plt.errorbar(header['AIRMASS'], header['airofavg'], yerr=header['airofstd'], color='b')
    plt.xlabel('BMDJ_TDB')
    #plt.xlabel('ROTATOR')
    plt.ylabel('ROTATOR')
    #plt.ylabel('pixel offset of air lines')
    plt.show()
    return



def plot_head_2_head(file_list, header1, header2,offset=10.):
    """
    Plot 2 arbitrary header values from a given set of files.
    
    """
    for filename in file_list:
        header= fits.getheader(filename)
        #print(filename, header['ROTATOR'],header['airofavg'], header['BMJD_TDB'])
        #print(filename, header['ROTATOR'],header['CAM_ANG'], header['CAM_TARG'])
        #plt.errorbar(header['BMJD_TDB'], header['airofavg'], yerr=header['airofstd'], color='b')
        #plt.errorbar(header['ROTATOR'], header['airofavg'], yerr=header['airofstd'], color='b', marker='o', markersize=4)
        #if header2.lower()=='cam_ang':
            #plt.scatter(header[header1], header['cam_targ']-header['cam_ang'],color='b')
            #plt.ylabel('cam_targ - cam_ang')
            #plt.xlabel(header1)
        #elif header1.lower()=='cam_ang':
            #plt.scatter( header['cam_targ']-header['cam_ang'], header[header2], color='b')
            #plt.xlabel('cam_targ - cam_ang')
        if 'a'=='b':
            pass
        elif ((header1.lower()=='rotator') and (header[header1] > 180)):
            plt.scatter( header[header1]-360,header[header2], color='b')
            plt.ylabel(header2)
            plt.xlabel(header1)
        elif ((header2.lower()=='rotator') and (header[header2] > 180)):
            plt.scatter(header[header1], header[header2]-360, color='b')
            plt.ylabel(header2)
            plt.xlabel(header1)
        #if header1.lower()=='mount_el':
            #plt.scatter(1./np.tan(header[header1]/180.*np.pi), header[header2],color='b')
        #if header2.lower()=='atm_refraction':
            #refraction_val=spt.calculate_atmospheric_refraction(header)
            #plt.scatter(header[header1], refraction_val, color='b')
            #plt.xlabel(header1)
            #plt.ylabel(header2+ ' (arcseconds)')
        #if header2.lower()=='delta_atm_refraction':
            #refraction_val=spt.calculate_atmospheric_refraction(header)
            #other_refraction_val=spt.calculate_atmospheric_refraction(header, offset=offset)
            #plt.scatter(header[header1], refraction_val-other_refraction_val, color='b')
            #plt.xlabel(header1)
            #plt.ylabel(header2+ ' (arcseconds) for '+ str(offset)+ ' arcminute offset')
        elif ((header1.lower()=='seeing' ) and (header2.lower()=='see_fwhm')):
            plt.scatter( header[header1],header[header2]*0.3, color='b')
            plt.ylabel(header2)
            plt.xlabel(header1)
        else:
            plt.scatter( header[header1],header[header2], color='b')
            plt.ylabel(header2)
            plt.xlabel(header1)
        #plt.scatter( header['AIRMASS'],header['ROTATOR'], color='b')
        #plt.errorbar(header['AIRMASS'], header['airofavg'], yerr=header['airofstd'], color='b')
    #plt.xlabel('ROTATOR')
    #plt.ylabel('pixel offset of air lines')
    #plt.grid()
    plt.show()
    
    
    return 


def plot_dheaddhead_2_dheaddhead(file_list,header1, dheader1, header2, dheader2):
    """
    Plot 2 discrete derivative type things for 2 arbitrary header values versus 2 other headers from a given set of files.
    
    """
    for number, filename in enumerate(file_list[:-1]):
        header= fits.getheader(filename)
        next_header=fits.getheader(file_list[number+1])
        show_delta_vals=False
        if header1.lower()=='rotator':
            if (next_header[header1]-header[header1])<-100:
                print('crossed 0 for rotator so slipping to underside')
                header1_numerator=((next_header[header1])-(header[header1]-360))
                show_delta_vals=True
            else:
                header1_numerator=(next_header[header1]-header[header1])
        else:
            pass
        delta_header1=header1_numerator/(next_header[dheader1]-header[dheader1])
        delta_header2=(next_header[header2]-header[header2])/(next_header[dheader2]-header[dheader2])
        if show_delta_vals:
            print(filename,'d '+header1+' /'+'d '+dheader1, delta_header1)
            print(filename,'d '+header2+' /'+'d '+dheader2, delta_header2)
        else:
            pass
        #print(filename, header['ROTATOR'],header['airofavg'], header['BMJD_TDB'])
        #print(filename, header['ROTATOR'],header['CAM_ANG'], header['CAM_TARG'])
        #plt.errorbar(header['BMJD_TDB'], header['airofavg'], yerr=header['airofstd'], color='b')
        #plt.errorbar(header['ROTATOR'], header['airofavg'], yerr=header['airofstd'], color='b', marker='o', markersize=4)
        #plt.scatter( header[header1],header[header2], color='b')
        plt.scatter(delta_header1,delta_header2, color='b')
        #plt.scatter( header['AIRMASS'],header['ROTATOR'], color='b')
        #plt.errorbar(header['AIRMASS'], header['airofavg'], yerr=header['airofstd'], color='b')
    plt.xlabel('d '+ header1+'/'+ 'd '+dheader1)
    #plt.xlabel('ROTATOR')
    plt.ylabel('d '+ header2+'/'+ 'd '+dheader2)
    #plt.ylabel(header2)
    #plt.ylabel('pixel offset of air lines')
    plt.show()
    
    
    
    return


def plot_dheaddhead_2_head(file_list,header1, header2,dheader2='null'):
    """
    Plot 2 discrete derivative type things for 2 arbitrary header values versus 2 other headers from a given set of files.
    
    """
    if dheader2=='null':
        dheader2=header1
    else:
        pass
    for number, filename in enumerate(file_list[:-1]):
        header= fits.getheader(filename)
        next_header=fits.getheader(file_list[number+1])
        show_delta_vals=False
        if header2.lower()=='rotator':
            if (next_header[header2]-header[header2])<-100:
                print('crossed 0 for rotator so slipping to underside')
                header2_numerator=((next_header[header2])-(header[header2]-360))
                show_delta_vals=True
            else:
                header2_numerator=(next_header[header2]-header[header2])
        else:
            header2_numerator=(next_header[header2]-header[header2])
        delta_header2=header2_numerator/(next_header[dheader2]-header[dheader2])
        if dheader2.lower()=='bmjd_tdb':
            delta_header2=delta_header2/u.day
            delta_header2=delta_header2.to(1/u.second).value
        else:
            pass
        if show_delta_vals:
            print(filename,'d '+header2+' /'+'d '+header1, delta_header2)
        else:
            pass
        #print(filename, header['ROTATOR'],header['airofavg'], header['BMJD_TDB'])
        #print(filename, header['ROTATOR'],header['CAM_ANG'], header['CAM_TARG'])
        #plt.errorbar(header['BMJD_TDB'], header['airofavg'], yerr=header['airofstd'], color='b')
        #plt.errorbar(header['ROTATOR'], header['airofavg'], yerr=header['airofstd'], color='b', marker='o', markersize=4)
        #plt.scatter( header[header1],header[header2], color='b')
        if ((header1.lower()=='rotator') and (header[header1])>180):
            plt.scatter(header[header1]-360,delta_header2, color='b')
        else:
            plt.scatter(header[header1],delta_header2, color='b')
        #plt.scatter( header['AIRMASS'],header['ROTATOR'], color='b')
        #plt.errorbar(header['AIRMASS'], header['airofavg'], yerr=header['airofstd'], color='b')
    plt.xlabel(header1)
    #plt.xlabel('ROTATOR')
    plt.ylabel('d '+ header2+'/'+ 'd '+dheader2)
    #plt.ylabel(header2)
    #plt.ylabel('pixel offset of air lines')
    plt.show()
    
    
    
    return


def plot_coordinates(file_list):
    ra_vals=[]
    dec_vals=[]
    time_vals=[]
    rotator_vals=[]
    for filename in file_list:
        header=fits.getheader(filename)
        coordinates=coords.SkyCoord(ra=header['ra'], dec=header['dec'], unit=(u.hourangle, u.degree))
        print(header['ra'],header['dec'])
        print(coordinates.ra.arcsecond, coordinates.dec.arcsecond)
        ra_vals.append(coordinates.ra.arcsecond)
        dec_vals.append(coordinates.dec.arcsecond)
        time_vals.append(header['BMJD_TDB'])
        if header['rotator']> 180:
            rotator_val=header['rotator']-360
        else:
            rotator_val=header['rotator']
        rotator_vals.append(rotator_val)
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



def plot_SNR(spec, noise, filename):
    header=fits.getheader(filename)
    center_pixel = np.argmin(np.abs(spec[0]-test_wavelength))
    measured_std = np.std(spec[1][int(center_pixel-test_side):int(center_pixel+test_side)])
    print("sigma in " + str(test_width) + " pixel range around " + str(test_wavelength)+ " angstroms", measured_std)
    sigma_range = noise[1][int(center_pixel-test_side):int(center_pixel+test_side)]
    print("Predicted sigmas of " + str(test_width) +" pixel range around " + str(test_wavelength)+ " angstroms", "min:" + str(np.min(sigma_range)), "mean:" + str(np.mean(sigma_range)), "max:" + str(np.max(sigma_range)))
    print("Mean S/N (per pixel) from Noise spec:", np.mean(spec[1]/noise[1]))
    print('Mean S/N (per pixel) from measured standard dev in target '+ str(test_wavelength)+'A range:',np.mean(spec[1][int(center_pixel-test_side):int(center_pixel+test_side)])/measured_std)
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

    
def plot_sky(filename, offset=0, line_labels=True, convolve=False, color='None',norm=False, norm_range=norm_range,divide_width=False,show_date=False):
    hdu=fits.open(filename)
    wavelengths=np.copy(hdu[0].data)
    label=filename+' sky'
    sky=np.copy(hdu[2].data)
    if divide_width:
        sky=sky/hdu[0].header['width']
        #label=label+' per pixel (width='+str(hdu[0].header['width'])+')'
        label=label+' per pixel'
    sky_spec=np.vstack([wavelengths, sky])
    if convolve:
        sky_spec= convolve_spectrum(sky_spec, 'dummy_header', kernel_type='box', pix_width=pix_width)
    else:
        pass
    if norm:
        #spec[1]= spec[1]/np.nanmean(spec[1])
        #spec[1]=spec[1]/np.nanmean(spec[1][1240:1280])
        #spec[1]=spec[1]/np.nanmean(spec[1][norm_range[0]:norm_range[1]])
        sky_spec=norm_spectrum(sky_spec, norm_range)
        sky= sky_spec[1]
        label=label+' normed'
    else:
        pass
        #sky= sky_spec[1]
        #sky=sky/np.nanmax(sky)
        #sky=sky+offset
    #sky= sky_spec[1]
    #sky=sky/np.nanmax(sky)
    #sky=sky+offset
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
    if show_date:
        label=label+' '+hdu[0].header['date-obs']
    else:
        pass
    if plot_wavelength:
        plt.xlabel(r'Wavelength ($\AA$)')
        if color== 'None':
            plt.plot(hdu[0].data, sky, label=label)
        else:
            plt.plot(hdu[0].data, sky, label=label, color=color)
        #plt.plot(hdu[0].data, sky, label=filename, color=color)
        #plt.plot(hdu[0].data, sky, label=hdu[0].header['airmass'])
    else:
        plt.plot(sky)
        plt.xlabel('pixel')
    plt.title(filename+' sky')
    #plt.show()
    return

def plot_diff_spec(spec1, spec2, filename1, filename2, header, smooth=False, kernel_type='gaussian', norm=False):
    #plot_spectrum(spec1, filename1, header, norm=norm, smooth=smooth, kernel_type=kernel_type)
    #plot_spectrum(spec2, filename2, header, norm=norm, smooth=smooth, kernel_type=kernel_type)
    title_string= filename1+' - ' + filename2
    if norm:
        #spec1[1]=spec1[1]/np.nanmean(spec1[1][norm_range[0]:norm_range[1]])
        #spec2[1]=spec2[1]/np.nanmean(spec2[1][norm_range[0]:norm_range[1]])
        spec1=norm_spectrum(spec1, norm_range)
        spec2=norm_spectrum(spec2, norm_range)
        spec2[1]=spec2[1]
        title_string= title_string+ ' normed'+str(norm_range)
        plt.ylabel('Flux (normed)')
        #spec1[1]=spec1[1]/np.nanmax(spec1[1])
        #spec2[1]=spec2[1]/np.nanmax(spec2[1])
    else:
        pass
    try:
        diff_flux=spec1[1]-spec2[1]
    except ValueError as error:
        print("ValueError:", error)
        print("Probably tried to subtract an SDSS spectrum from a Goodman spectrum, so we need to rebin")
        greater_spec_index=np.greater(spec1[1].shape[0], spec2[1].shape[0])
        print(spec1[1].shape[0], spec2[1].shape[0],np.greater(spec1[1].shape[0], spec2[1].shape[0]))
        if greater_spec_index==False:
            spec1_hdu=fits.open(filename1)
            spec1_dlambda=spec1_hdu[4].data
            spec2_dlambda=np.copy(np.roll(spec2[0], -1) - spec2[0])
            spec2=spt.rebin_generic_spec(spec2[:,:-1], spec2_dlambda[:-1], spec1[0], spec1_dlambda)
            diff_flux=spec1[1]-spec2[1]
        else:
            spec2_hdu=fits.open(filename2)
            spec2_dlambda=spec2_hdu[4].data
            spec1_dlambda=np.copy(np.roll(spec1[0], -1) - spec1[0])
            spec1=spt.rebin_generic_spec(spec1[:,:-1], spec1_dlambda[:-1], spec2[0], spec2_dlambda)
            diff_flux=spec1[1]-spec2[1]
    diff_spec=np.vstack([spec1[0], diff_flux])
    plot_spectrum(diff_spec, filename1+' - ' + filename2, header, smooth=smooth, kernel_type=kernel_type)     
    #plot_spectrum(diff_spec, title_string, header, smooth=smooth, kernel_type=kernel_type)
    plt.ylim(top=np.nanmax([spec1[1],spec2[1]])+0.5)
    plt.ylim(bottom=np.nanmin(diff_flux)-0.5)
    return


def plot_white_lightcurve(filenames,header_name='BMJD_TDB',wave_range=[0.,30000.]):
    """
    Produce a white-light curve from spectra by integrating flux and plotting versus time.
    
    
    """
    times=[]
    white_fluxes=[]
    for filename in filenames:
        hdu=fits.open(filename)
        header=fits.getheader(filename)
        fluxes=np.copy(hdu[1].data)
        dlambda=np.copy(hdu[4].data)
        ###
        inbounds=np.where((hdu[0].data>wave_range[0]) & (hdu[0].data < wave_range[1]))
        fluxes=fluxes[inbounds]
        dlambda=dlambda[inbounds]
        
        ####
        white_light=np.sum(fluxes*dlambda)
        white_fluxes.append(white_light)
        print('white_light', white_light)
        if ((header_name.lower()=='rotator') and (header[header_name])>180):
            times.append(header[header_name]-360)
            xlabel=header_name
        elif header_name.lower()=='mount_el':
           times.append(1./np.tan(header[header_name]/180.*np.pi))
           xlabel='cot('+header_name+')'
        elif header_name.lower()=='cam_ang':
            #plt.scatter( header['cam_targ']-header['cam_ang'], header[header2], color='b')
            #plt.xlabel('cam_targ - cam_ang')
            times.append(header['cam_targ']-header['cam_ang'])
            xlabel=header_name
        else:
            xlabel=header_name
            times.append(header[header_name])
    
    plt.plot(times, white_fluxes, marker='o',linestyle='None')
    plt.xlabel(xlabel)
    plt.ylabel('White Light Fluxes (erg/cm/cm/s)')
    plt.show()
    return




if __name__ == '__main__':

    if file_setting== 'compare_SDSS':
        fbadpix_list=[]
        print('\nCompare SDSS selected\n\n')
        target_spec1, header1, target_noise1= spt.retrieve_spec(filename)
        filename1=filename
        #target_spec1, header1, target_noise1= spt.retrieve_sdss_spec(filename)
        target_spec1[0]=target_spec1[0]+wavelength_offset
        mwdd_spec=spt.retrieve_mwdd_spec(mwdd_filename,convert_to_flambda=False)
        #mwdd_spec=spt.retrieve_mwdd_spec(mwdd_filename, convert_to_flambda=True)
        print('retrieved MWDD spec')
        #target_spec1= norm_spectrum(target_spec1, norm_range)
        #for filename2 in sdss_names:
        print('sdss_names',sdss_names)
        for  index, filename2 in enumerate(sdss_names):
            if save_plots:
                plt.figure(figsize=(1036./540*6.,6.))
            else:
                pass
            target_spec2, header2, target_noise2= spt.retrieve_sdss_spec(filename2, wave_medium='air')
            
            
            #fbadpix_list.append(fits.getheader(filename2)['fbadpix'])
            #print('fbadpix_list',fbadpix_list)
            
            #target_spec2= spt.clean_spectrum(target_spec2, np.nanmin(target_spec1[0]), np.nanmax(target_spec1[0]), [])
            #target_spec2=norm_spectrum(target_spec2, norm_range)
            #target_spec2[1]=target_spec2[1]*sdss_flux_scale_factor
            #plt.ylim(top=np.nanmax(np.hstack([target_spec2[1], target_spec1[1]]))+0.5)
            #plt.ylim(top=np.nanpercentile(np.hstack([target_spec2[1], target_spec1[1]]), 99.9)+0.5)
            #plot_spectrum(target_spec1, filename, header1, norm=True, smooth=True, kernel_type='gaussian')
            #plot_spectrum(target_spec1, filename, header1, norm=False, smooth=True, kernel_type='box')
            #plot_spectrum(target_spec1, filename, header1, norm=True, smooth=True, kernel_type='box', pix_width=sdss_pix_width)
            
            #plot_spectrum(target_spec2, filename2.split('/')[-1], header2, norm=True, smooth=True, kernel_type='box', pix_width=sdss_pix_width, color='k')
            
            #plot_telluric_spectrum([3700, 9000], smooth=True, pix_width=30, color='r')
            
            #plot_spectrum(target_spec2, filename2, header1, norm=True, smooth=True, kernel_type='gaussian', pix_width=8.*header1['see_sig'],kernel_width=10.*slit_width)
            
            #plot_spectrum(target_spec2, filename2, header1, norm=True, smooth=True, kernel_type='gaussian', pix_width=8.*header1['see_sig'],kernel_width=10.*slit_width)
            
            
            
            #plot_spectrum(target_spec2, filename2.split('/')[-1], header1, norm=True, smooth=True, kernel_type='box', pix_width=sdss_pix_width,color='k')
            
            
            #plot_spectrum(target_spec1, filename, header1, norm=False, smooth=True, kernel_type='box', pix_width=5,alpha=1.0,color='k')
            
            #plot_spectrum(target_spec2, filename2.split('/')[-1], header1, norm=True, smooth=True, kernel_type='box', pix_width=sdss_pix_width*sdss_scale_factor)
            
            #plot_spectrum(target_spec2, filename2.split('/')[-1], header1, norm=True, smooth=True, kernel_type='box', pix_width=5*sdss_scale_factor)
            #for mwdd_filename in mwdd_list:
                        #mwdd_spec=spt.retrieve_mwdd_spec(mwdd_filename,convert_to_flambda=False)
                        #plot_spectrum(mwdd_spec, mwdd_filename, header1, norm=True, smooth=True, kernel_type='box',pix_width=10)
            plot_spectrum(mwdd_spec, mwdd_filename, header1, norm=True, smooth=True, kernel_type='box',pix_width=10)
            plot_spectrum(target_spec1, filename, header1, norm=True, smooth=True, kernel_type='box', pix_width=5)
            bb_spec=spt.blackbody_spec(target_spec1[0],bb_teff)
            print(bb_spec)
            plot_spectrum(bb_spec, str(bb_teff)+' K BB', header1, norm=True, smooth=False)
                    
            model_file_list=sorted(glob(cp.balmer_model_dir+'*.txt'))
            for model_filename in model_file_list:
                model_spec=spt.retrieve_blouin_model(model_filename)
                trim_model_spec=spt.clean_spectrum(model_spec,3600,7000,[])
                norm_model_spec=norm_spectrum(trim_model_spec,norm_range)
                label=model_filename.split('/')[-1]
                plt.plot(norm_model_spec[0],norm_model_spec[1],label=model_filename.split('/')[-1])
            spt.plot_wien_max_wave(np.arange(4500,6600, 500))
            
            #plot_spectrum(spt.flambda_to_fnu(target_spec2), filename2.split('/')[-1], header1, norm=True, smooth=True, kernel_type='box', pix_width=5*sdss_scale_factor)
            ##plot_spectrum(spt.flambda_to_fnu(mwdd_spec), mwdd_filename, header1, norm=True, smooth=True, kernel_type='box',pix_width=5)
            #plot_spectrum(spt.flambda_to_fnu(target_spec1), filename, header1, norm=True, smooth=True, kernel_type='box', pix_width=5)
            #plt.ylabel('fnu normalized')

            #plot_spectrum(target_spec2, filename2.split('/')[-1], header1, norm=False, smooth=True, kernel_type='sdss_match', pix_width=pix_width,color='r',alpha=0.8)
            
            
            ####### plotting individual SDSS spectra frames too#####
            
            
            #sdss_spec_list, sdss_header_list, sdss_noise_list=spt.retrieve_allframes_sdss_spec(filename2, wave_medium='air')
            #extension=4
            #total_extensions=len(sdss_spec_list)
            #try:
                #for sdss_spec,sdss_header in zip(sdss_spec_list, sdss_header_list):
                    #try:
                        #date=sdss_header['date-obs']+' '+sdss_header['taihms']
                    #except KeyError as error:
                        #print("KeyError:",error)
                        #date=str(0)
                    ##datetime=Time(sdss_header['tai'],format='tai')
                    ##testtai=Time(0,format='tai')
                    ##print('testtai', testtai.tai)
                    ##print('testtai', testtai.utc)
                    #if extension<total_extensions/2+4:
                        #color='b'
                    ##if extension >total_extensions/2+4:
                    #else:
                        #color='r'
                    ##if ((extension==7) or (extension==11)):
                        ##color='g'
                    #plot_spectrum(sdss_spec, filename2.split('/')[-1]+str(extension)+' '+date, header1, norm=False, smooth=True, kernel_type='sdss_match', pix_width=pix_width, alpha=0.4,color=color)
                    ##plot_spectrum(sdss_spec, filename2.split('/')[-1]+str(extension)+' '+str(sdss_header['tai-beg']), header1, norm=False, smooth=True, kernel_type='sdss_match', pix_width=pix_width, alpha=0.8,color=color)
                    #extension+=1
            #except IndexError as error:
                    #print("IndexError:", error)
                    #print("\n\nFile that lacked the extensions:",filename2,"\n\n")
            
            ###### end plotting individual sdss frames#########
            
            #plot_spectrum(target_spec2, filename2.split('/')[-1]+' coadd', header1, norm=False, smooth=True, kernel_type='sdss_match', pix_width=pix_width,color='k',alpha=1.0)
            
            #filename3=sdss_names[index+1]
            #filename3=filename2
            #target_spec3, header3, target_noise3=spt.retrieve_sdss_spec(filename3,wave_medium='air')
            #plot_spectrum(target_spec3, filename3.split('/')[-1], header1, norm=True, smooth=False, kernel_type='sdss_match', pix_width=pix_width,color='b',alpha=0.4)
            
            
            
            #target_spec2= norm_spectrum(target_spec2, norm_range)
            #target_spec2[1]=target_spec2[1]+0.5
            #plt.axhline(y=0.5, color='b',linestyle=':')
            #plot_spectrum(target_spec2, filename2.split('/')[-1]+' +0.5 flat flux', header1, norm=True, smooth=True, kernel_type='sdss_match', pix_width=pix_width,color='r')
            
            #plot_sky(filename, offset=0, line_labels=False, convolve=False)
            #plot_spectrum(target_spec2, filename2, header2, norm=False, smooth=True, kernel_type='box', pix_width=sdss_pix_width)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=False)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=True)
            #plot_diff_spec(target_spec2, target_spec1, filename2, filename1, header1, smooth=True, norm=True)
            #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=False)
            plt.axhline(y=0, linestyle='--', color='k')
            #plt.legend()
            plt.title(filename+ ' & '+ filename2.split('/')[-1])
            #plt.xlim(np.nanmin(target_spec1[0]), np.nanmax(target_spec1[0]))
            #plt.show()
            spt.show_plot(line_id='h',convert_to_air=True,actually_show=False)
            spt.show_plot(line_id='cool_wd',convert_to_air=True,actually_show=True)
            plt.show()
            #spt.show_plot(line_id='h',convert_to_air=True, actually_show=True)
            #spt.show_plot(line_id='molecules',convert_to_air=False,label_pos=1.05)
            
            #plt.ylim(0,np.percentile(np.append(target_spec1[1],target_spec2[1],axis=0),0.9)*1.05)
            if save_plots:
                output_name=filename2.split('/')[-1]
                output_name=output_name.split('.')[0]
                output_name=output_name+'.pdf'
                plt.savefig(output_dir+output_name)
                plt.clf()
            else:
                pass
        #plt.title('FBADPIX header values for all spectra')
        #plt.hist(fbadpix_list,bins=np.arange(0,0.1,0.001))
        #plt.show()
            
    if file_setting=='all_SDSS':
        for filename1 in filenames:
            print('doing first file')
            target_spec1, header1, target_noise1= spt.retrieve_sdss_spec(filename1)
            target_spec1=norm_spectrum(target_spec1, norm_range)
            for filename2 in sdss_names:
                print('doing second file',filename2)
                target_spec2, header2, target_noise2= spt.retrieve_sdss_spec(filename2)
                #target_spec2= spt.clean_spectrum(target_spec2, np.nanmin(target_spec1[0]), np.nanmax(target_spec1[0]), [])
                target_spec2=norm_spectrum(target_spec2, norm_range)
                #plt.ylim(top=np.nanpercentile(np.hstack([target_spec2[1], target_spec1[1]]), 99.9)+0.5)
                target_dlambda=target_spec1[0][400]-target_spec1[0][399]
                template_dlambda=target_spec2[0][400]-sdss_spec2[0][399]
                print('target_dlambda',target_dlambda)
                print('template_dlambda',template_dlambda)
                print('target kernel',main_kernel)
                template_kernel=sdss_pix_width*(target_dlambda/template_dlambda)
                print('template_kernel',template_kernel)
                
                plot_spectrum(target_spec1, filename1.split('/')[-1], header1, norm=True, smooth=True, kernel_type='box', pix_width= sdss_pix_width)
                plot_spectrum(target_spec2, filename2.split('/')[-1], header2, norm=True, smooth=True, kernel_type='box', pix_width=template_kernel,alpha=0.7)
            
                #tell_filename='LBL_A30_s0_w005_R0060000_T.fits'
                
                
                
                #smooth_spec1=convolve_spectrum(target_spec1, header1, kernel_type='box', pix_width=sdss_pix_width, kernel_width=slit_width)
                #smooth_spec2=convolve_spectrum(target_spec2, header2, kernel_type='box', pix_width=sdss_pix_width, kernel_width=slit_width)
                
                
                
                
                #plt.plot((const.c/(target_spec1[0]*u.angstrom)).to(u.Hz),target_spec1[1],label=filename1.split('/')[-1])
                #plt.plot((const.c/(target_spec2[0]*u.angstrom)).to(u.Hz),target_spec2[1],label=filename2.split('/')[-1])
                #plt.plot((const.c/(target_spec1[0]*u.angstrom)).to(u.Hz),spt.flambda_to_fnu(smooth_spec1)[1],label=filename1.split('/')[-1])
                #plt.plot((const.c/(target_spec2[0]*u.angstrom)).to(u.Hz),spt.flambda_to_fnu(smooth_spec2)[1],label=filename2.split('/')[-1])
                
                #plt.plot(target_spec1[0],spt.flambda_to_fnu(smooth_spec1)[1],label=filename1.split('/')[-1])
                #plt.plot(target_spec2[0],spt.flambda_to_fnu(smooth_spec2)[1],label=filename2.split('/')[-1])
                
                #rv=1000
                #target_spec2[0]=spt.get_doppler_shifted(target_spec2[0],rv)
                #plt.plot(target_spec2[0],spt.flambda_to_fnu(smooth_spec2)[1],label=filename2.split('/')[-1]+' w/ RV of'+str(rv)+' km/s')
                
                #wave_shift=-151
                #plt.plot(target_spec2[0]+wave_shift,spt.flambda_to_fnu(smooth_spec2)[1],label=filename2.split('/')[-1]+' w/ shift of ' + str(wave_shift) + ' A')
                
                ###tell_filename='LBL_A30_s0_w200_R0060000_T.fits'
                ###plot_telluric_spectrum([3700, 9000], smooth=True, pix_width=30, tell_filename=tell_filename)
                
                tell_filename='LBL_A30_s0_w200_R0060000_T.fits'
                #plot_telluric_spectrum([3700, 9000], smooth=True, pix_width=30,tell_filename=tell_filename)
                #plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=True, norm=True)
                
                ####plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=True)
                ####plot_diff_spec(target_spec1, target_spec2, filename1, filename2, header1, smooth=False, norm=False)        plt.axhline(y=0, linestyle='--', color='k')
                ####plt.legend()
                
                
                #plt.title(filename1+ ' & '+ filename2.split('/')[-1])
                #plt.show()
                #plt.legend()
                #plt.show()
                #spt.show_plot(line_id='alkali', convert_to_air=False)
                #spt.show_plot(line_id='mystery', convert_to_air=False)
                #spt.show_plot(line_id='cyclotron3900', convert_to_air=False)
                #spt.show_plot(line_id='cool_wd', convert_to_air=False)
                #spt.show_plot(line_id='C2_bands',convert_to_air=False)
                plt.ylim(0,np.percentile(np.append(target_spec1[1],target_spec2[1],axis=0),99.5)*1.05)
                spt.show_plot(line_id='cool_wd',convert_to_air=True,actually_show=False)
                spt.show_plot(line_id='h',convert_to_air=True, actually_show=True)
                
        #plt.legend()
        plt.show()

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
            plot_sky(m1_name, offset=0,line_labels=False)
            plot_sky(m2_name, offset=0,line_labels=False)
            
            #plot_spectrum(target_spec1, m1_name, header1, norm=True, smooth=False, kernel_type='box')
            plot_spectrum(target_spec2, m2_name, header2, norm=True, smooth=True, kernel_type='gaussian',pix_width=header1['see_sig'])
            plot_spectrum(target_spec1, m1_name, header1, norm=True, smooth=True, kernel_type='gaussian', pix_width=header2['see_sig'])
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
        #spt.show_plot(line_id='alkali', convert_to_air=True)
        spt.show_plot(show_legend=True, line_id='cool_wd', convert_to_air=True)
        spt.show_plot()
            
    if file_setting=='two_arm_compare_SDSS':
        target_spec1, header1, target_noise1= spt.retrieve_spec(filename1)
        target_spec2, header2, target_noise2= spt.retrieve_spec(filename2)
        mwdd_spec=spt.retrieve_mwdd_spec(mwdd_filename)
        
        target_spec1[0]=spt.air_to_vac(target_spec1[0])
        target_spec2[0]=spt.air_to_vac(target_spec2[0])
        
        #gaia_spec, gaia_noise=spt.retrieve_gaia_xp_spec('LHS2534')
        
        #target_spec1= spt.flambda_to_fnu(target_spec1)
        #target_spec2=spt.flambda_to_fnu(target_spec2)
        
        #target_spec1=norm_spectrum(target_spec1, norm_range)
        #target_spec2= norm_spectrum(target_spec2, norm_range)
        
        for sdss_filename in sdss_names:
            
            sdss_spec, sdss_header, sdss_noise2= spt.retrieve_sdss_spec(sdss_filename, wave_medium='vac')
            
            #sdss_spec[0]=spt.vac_to_air(sdss_spec[0])
            
            sdss_spec= spt.clean_spectrum(sdss_spec, 3700, 9000, [])
            
            
            #sdss_spec=spt.flambda_to_fnu(sdss_spec)
            
            #sdss_spec= norm_spectrum(sdss_spec, norm_range)
            
            
            print('header1 see_sig', header1['SEE_SIG'])
            print('header2 see_sig', header2['SEE_SIG'])
            print('median dlambda spec1:', get_median_dlambda(target_spec1))
            print('median dlambda spec2:', get_median_dlambda(target_spec2))
            print('median dlambda sdss spec:', get_median_dlambda(sdss_spec))
            print('scaling factor:', get_median_dlambda(target_spec1)/get_median_dlambda(sdss_spec))
            print('scaling factor:', get_median_dlambda(target_spec2)/get_median_dlambda(sdss_spec))
            plt.ylim(top=np.percentile(np.hstack([target_spec2[1], target_spec1[1], sdss_spec[1]]), 99.9)+0.5)
            
            #plot_spectrum(sdss_spec, sdss_filename.split('/')[-1], header2, norm=True, smooth=True, kernel_type='gaussian', pix_width=pix_width, color='r')
            plot_spectrum(sdss_spec, sdss_filename.split('/')[-1], header1, norm=True, smooth=True, kernel_type='sdss_match', pix_width=pix_width, color='r')
            #plot_spectrum(target_spec1, filename1, header1, pix_width=header1['see_sig'],norm=True, smooth=True, kernel_type='gaussian',  color='g')
            #plot_spectrum(target_spec2, filename2, header2, norm=True, smooth=True, pix_width=header2['see_sig'], kernel_type='gaussian',  color='b')
            
            plot_spectrum(target_spec1, filename1, header1, pix_width=5,norm=True, smooth=True, kernel_type='box',  color='g')
            plot_spectrum(target_spec2, filename2, header2, norm=True, smooth=True, pix_width=5, kernel_type='box',  color='b')
            
            plot_spectrum(mwdd_spec, mwdd_filename, header1, norm=True, smooth=False, color='k')

            #plot_spectrum(gaia_spec,'LHS2534 GaiaDR3 XP',header1, smooth=False, norm=False)
            
            #plot_dwavelength(target_spec1, filename1, read_in=False)
            #plot_dwavelength(target_spec2, filename2, read_in=False)
            #plot_dwavelength(sdss_spec, sdss_filename.split('/')[-1], read_in=False)
            #plot_telluric_spectrum([3700, 9000], smooth=True, pix_width=30)
            plt.xlim(3700,9000)
            #plt.ylabel(r'$F_{\nu}$ (normalized)')
            spt.show_plot(line_id='cool_wd', convert_to_air=False)
            #spt.show_plot(line_id='mystery', convert_to_air=False)
            
    if single_iterate:
        counter=0
        ew_list=[]
        ew_noise_list=[]
        bmjd_list=[]
        sig_list=[]
        snr_list=[]
        print('\n\n',filenames,'\n\n')
        for filename in filenames:
            target_spec, header, target_noise= spt.retrieve_spec(filename)
            hdu= fits.open(filename)
            sig_list.append(np.std(target_spec[1]))
            #print(filename)
            #gaia_spec, gaia_noise=spt.retrieve_gaia_xp_spec('LTT3218')
            ####Na-D line for J0850
            ew,ew_noise=spt.get_ew(filename,[5728.,6068.], cont_method='poly1',cont_width=120.,noise_method='rms', plot_all=False )
            
            #H-beta for GD71
            #ew,ew_noise=spt.get_ew(filename,[4718.,5008.], cont_method='poly1',cont_width=120.,noise_method='rms', plot_all=True )

            #H-alpha for GD71
            #ew,ew_noise=spt.get_ew(filename,[6448.,6678.], cont_method='poly1',cont_width=120.,noise_method='rms', plot_all=True )
            
            #ew,ew_noise=spt.get_ew(filename,[4150.,4302.], cont_method='poly1',cont_width=60.,noise_method='rms', plot_all=False )
            
            #ew,ew_noise=spt.get_ew(filename,[4190.,4262.], cont_method='poly1',cont_width=60.,noise_method='rms', plot_all=False )
            
            #Li I line for J0850
            #ew,ew_noise=spt.get_ew(filename,[6678,6738], cont_method='poly1',cont_width=120.,noise_method='rms', plot_all=False )
            
              ##Ca I line
            #ew,ew_noise=spt.get_ew(filename,[4200,4258], cont_method='poly1',cont_width=120.,noise_method='rms', plot_all=True )
            
            ####H-alpha
            #ew,ew_noise=spt.get_ew(filename,[6513.,6613.], cont_method='poly1',cont_width=120.,noise_method='rms', plot_all=False )
            
            #mystery emission
            #ew,ew_noise=spt.get_ew(filename,[6255.,6295.], cont_method='poly1',cont_width=120.,noise_method='rms', plot_all=False )
            
            ew_list.append(ew)
            ew_noise_list.append(ew_noise)
            bmjd_list.append(header['bmjd_tdb'])
            noise_inds=np.where((target_spec[0]>5166.)&(target_spec[0]<5532.))
            this_noise=np.std(target_spec[1][noise_inds])
            this_snr=np.mean(target_spec[1][noise_inds])/this_noise
            snr_list.append(this_snr)
            #print(header['senscurv'], header['seeing'])
            #plot_telluric_spectrum([3700, 9000], smooth=True, pix_width=30, color='r')
            print('seeing from trace:',header['see_FWHM'])
            print('\n',header['OPENDATE'])
            print('see_FWHM*0.3"', header['see_FWHM']*0.3)
            print('seeing actual header',header['seeing'])
            print('see_sig',header['see_sig'])
            print('extraction_width',header['width'])
            print('width in units of see_sig', header['width']/header['see_sig'])
            #print('bkgshift',header['bkgshift'], 'bkgwidth',header['bkgwidth'])
            #print('see_FWHM*0.3/DIMM_seeing"', header['see_FWHM']*0.3/header['seeing'])
            #print('DIMM seeing', header['seeing'],'\n')
            dlambda= hdu[4].data
            #target_spec[0]=target_spec[0]+wavelength_offset
            #print(filename, 'mean: ', np.nanmean(target_spec[1]))
            nu_spec= spt.flambda_to_fnu(target_spec, dlambda)
            #conv_spec= convolve_spectrum(target_spec, header)
            #plot_spectrum(target_spec, filename, header, smooth=True, norm=False, kernel_type='box')
            #plt.scatter(header['BMJD_TDB'], np.sum(target_spec[1]*dlambda))
            #plt.errorbar(target_spec[0], target_spec[1], yerr=target_noise[1], label=filename, marker='o')
            #plot_spectrum(nu_spec, 'fnu', header, smooth=True, norm=True, kernel_type='box')
            #plot_spectrum(target_spec, filename, header, smooth=True, norm=True, pix_width=header['see_sig'], kernel_type='gaussian')
            
            #plot_spectrum(target_spec, filename, header, smooth=True , norm=True, pix_width=5., kernel_type='box',alpha=1.0)
            #plot_spectrum(target_spec, filename+' smoothed', header, smooth=True , norm=False, pix_width=5., kernel_type='box',alpha=0.7)
            #plot_spectrum(gaia_spec,'LTT3218 GaiaDR3 XP',header, smooth=False, norm=False)
            
            #plot_spectrum(target_spec, filename, header, smooth=True, norm=True, pix_width=5., kernel_type='box',alpha=1.0,offset=(counter%6.)*1.)
            
            #plot_spectrum(spt.flambda_to_fnu(target_spec), filename, header, smooth=True, norm=True, pix_width=0.5*header['see_sig'], kernel_type='gaussian',alpha=1.0)
            
            #plot_spectrum(target_spec, filename, header, smooth=True, norm=True, pix_width=1., kernel_type='gaussian',alpha=1.0)
            
            #Arbitrary Resolution matching smoothing method below, but you have to specify what resolution is the target in pix_width as the number whose square has the file's see_sig's square subtracted from it.
            #plot_spectrum
            
            ######################################################
            #################################################
            #####Poorly placed trial run of tweaking of the flux calibration based on deviations of warm DCs from a blackbody curve.

            ####plot_spectrum(target_spec, filename, header, smooth=False, norm=False, pix_width=5, kernel_type='box')
            ####binned_spec=spt.rebin_generic_spec(target_spec,hdu[4].data,np.arange(np.nanmin(target_spec[0]),np.nanmax(target_spec[0]),30.),np.ones(np.arange(np.nanmin(target_spec[0]),np.nanmax(target_spec[0]),30.).shape)*30.)
            #####yakima=scinterp.Akima1DInterpolator(target_spec[0],target_spec[1])(target_spec[0])
            ####yakima=scinterp.Akima1DInterpolator(binned_spec[0],binned_spec[1])(target_spec[0])
            ####plot_spectrum(binned_spec, '30 A-binned spec', header, smooth=False, norm=False, pix_width=5, kernel_type='box')
            ####plt.plot(target_spec[0],yakima,label='Akima Interp')
            ####bb_spec=spt.blackbody_spec(target_spec[0],bb_teff)
            ####print(bb_spec)  
            ####plot_spectrum(bb_spec, str(bb_teff)+' K BB', header, norm=True, smooth=False)
            ####spt.show_plot(show_legend=True, line_id='cool_wd', convert_to_air=True,actually_show=False)
            ####spt.show_plot(show_legend=False, line_id='h', convert_to_air=True, actually_show=True)
            
            ####norm_bb_spec=norm_spectrum(bb_spec,norm_range,show_norm_range=False)
            ####norm_interp_spec=norm_spectrum(np.vstack([target_spec[0],yakima]),norm_range,show_norm_range=False)
            
            ####plt.plot(norm_bb_spec[0],norm_bb_spec[1],label='normed BB spec')
            ####plt.plot(norm_interp_spec[0],norm_interp_spec[1],label='normed interp spec')
            ####plt.xlabel('Wavelength Angstroms')
            ####plt.legend()
            ####plt.show()
            
            ####div_flux=norm_bb_spec[1]/norm_interp_spec[1]
            ####plt.plot(target_spec[0],div_flux,label='norm bb spec/norm interp spec')
            ####plt.xlabel('Wavelength Angstroms')
            ####plt.legend()
            ####plt.show()
            
            ####plot_spectrum(target_spec, filename, header, smooth=False, norm=True, pix_width=5, kernel_type='box')
            ####binned_spec=spt.rebin_generic_spec(target_spec,hdu[4].data,np.arange(np.nanmin(target_spec[0]),np.nanmax(target_spec[0]),30.),np.ones(np.arange(np.nanmin(target_spec[0]),np.nanmax(target_spec[0]),30.).shape)*30.)
            #####yakima=scinterp.Akima1DInterpolator(target_spec[0],target_spec[1])(target_spec[0])
            ####yakima=scinterp.Akima1DInterpolator(binned_spec[0],binned_spec[1])(target_spec[0])
            ####plot_spectrum(binned_spec, '30 A-binned spec', header, smooth=False, norm=True, pix_width=5, kernel_type='box')
            ####plot_spectrum(np.vstack([target_spec[0],yakima]), 'Interp Spectrum', header, smooth=False, norm=True, pix_width=5, kernel_type='box')
            #####plt.plot(target_spec[0],yakima,label='Akima Interp')
            ####bb_spec=spt.blackbody_spec(target_spec[0],bb_teff)
            ####print(bb_spec)  
            ####plot_spectrum(bb_spec, str(bb_teff)+' K BB', header, norm=True, smooth=False)
            ####corrected_spec=np.vstack([target_spec[0],target_spec[1]*div_flux])
            ####plot_spectrum(corrected_spec, 'Hopefully corrected w/ leftovers', header, smooth=False, norm=True, pix_width=5, kernel_type='box')
            #####plt.plot(target_spec[0],target_spec[1]*div_flux,label='Hopefully corrected w/ leftovers')
            ####spt.show_plot(show_legend=True, line_id='cool_wd', convert_to_air=True,actually_show=False)
            ####spt.show_plot(show_legend=False, line_id='h', convert_to_air=True, actually_show=True)
            
            ####plot_spectrum(target_spec, filename, header, smooth=False, norm=False, pix_width=5, kernel_type='box')
            ####plot_spectrum(corrected_spec, 'Hopefully corrected w/ leftovers', header, smooth=False, norm=False, pix_width=5, kernel_type='box')
            ####spt.show_plot(show_legend=True, line_id='cool_wd', convert_to_air=True,actually_show=False)
            ####spt.show_plot(show_legend=False, line_id='h', convert_to_air=True, actually_show=True)
            
            ####other_filename='ravg_fwctb.WDJ1828m6357_400m1.fits'
            ####for other_filename in sorted(glob('ravg_fwctb.*fits')):
                ####other_spec, other_header, other_noise= spt.retrieve_spec(other_filename)
                ####correct_other_spec=np.vstack([other_spec[0],other_spec[1]*(np.interp(other_spec[0],target_spec[0],div_flux))])
                ####plot_spectrum(other_spec, other_filename, other_header, smooth=True, norm=False, pix_width=5, kernel_type='box')
                ####plot_spectrum(correct_other_spec, other_filename+' corrected with previous method', other_header, smooth=True, norm=False, pix_width=5, kernel_type='box')
                ####mwdd_spec=spt.retrieve_mwdd_spec(mwdd_filename,convert_to_flambda=False)
                ####plot_spectrum(mwdd_spec, mwdd_filename, header, norm=False, smooth=True, kernel_type='box',pix_width=10)
                ####spt.show_plot(show_legend=True, line_id='cool_wd', convert_to_air=True,actually_show=False)
                ####spt.show_plot(show_legend=False, line_id='h', convert_to_air=True, actually_show=True)
                
            #####End of poorly positioned flux correction code that needs to be spun off
            ####################################################
            ####################################################
            
            
            plot_spectrum(target_spec, filename, header, smooth=True, norm=False, pix_width=5, kernel_type='box')
            #plot_sky(filename, offset=0, line_labels=False, convolve=False,norm=True)
            #spt.show_plot(show_legend=True, line_id='h')
            
            #plot_spectrum(target_spec, filename, header, smooth=True, norm=True, pix_width=0.5*header['see_sig'], kernel_type='gaussian', offset=counter*0.1)
            #plot_spectrum(target_spec, filename+', seeing '+str(header["see_FWHM"]*header["pix_scal"])[:4]+ '" '+str(header['width']), header, smooth=True, norm=False, pix_width=5, kernel_type='box')
            #try:
                #hold_spec[1]=hold_spec[1]-target_spec[1]
                #plot_spectrum(hold_spec, filename+', seeing '+str(header["see_FWHM"]*header["pix_scal"])[:4]+ ' " ', header, smooth=True, norm=False, pix_width=5, kernel_type='box')
                #hold_spec=target_spec
            #except NameError as error:
                #print("NameError:", error)
                #hold_spec=target_spec
            #plot_spectrum(nu_spec, filename, header, smooth=True, norm=True, pix_width=15, kernel_type='box')

            #plot_spectrum(target_spec, filename, header, smooth=True, norm=False,  kernel_type='gaussian',pix_width=header['SEE_SIG'])
            #target_spec[1]=header['airmass']
            #plot_spectrum(target_spec, filename, header, norm=False, smooth=True, kernel_type='box', pix_width=10)
            #plot_spectrum(target_spec, str(header['airmass']), header, norm=False, smooth=True, kernel_type='box', pix_width=10)
            #plot_spectrum(target_spec, filename, header, norm=True, smooth=True, pix_width=15,offset=counter*-0.2, kernel_type='box')
            #plot_spectrum(target_spec, filename, header, norm=True, smooth=True, kernel_type='box', offset=counter)
            #plot_spectrum(nu_spec, filename, header, norm=False, smooth=True, kernel_type='box')
            #plt.plot(target_spec[0], dlambda,  label=filename, marker='o', markersize=10-counter)
            ##plot_spectrum(target_spec, filename, header, smooth=True, kernel_type='gaussian', norm=True)
            #plot_sky(filename, offset=0, line_labels=False, convolve=False,norm=False)
            #if header['bmjd_tdb']>   60560.:
                #plot_sky(filename, offset=0, line_labels=False, convolve=False,norm=False, divide_width=True, show_date=True)
            #plot_sky(filename, offset=0, line_labels=False, convolve=False,norm=True, divide_width=False, show_date=True)
            #if header['airmass']<1.5:
                #plot_sky(filename, offset=0)
            #else:
                #pass
            #plot_SNR(target_spec, target_noise, filename)
            #plot_SNR_from_file(filename)
            #plot_dwavelength(target_spec, filename)
            #plt.axvline(x=6563., linestyle=':',color='k')
            #plot_telluric_spectrum([3700, 9000], smooth=True, pix_width=30)
            #plot_telluric_spectrum([3700,9000], smooth=True, pix_width=30, tell_filename='LBL_A30_s0_w200_R0060000_T.fits')
            #spt.show_plot(show_telluric=False, show_legend=False)
            
            #spt.show_plot(show_legend=True, line_id='cool_wd', convert_to_air=True,actually_show=False)
            #spt.show_plot(show_legend=False, line_id='h', convert_to_air=True, actually_show=True)
            #spt.show_plot(show_legend=True, line_id='cool_wd', convert_to_air=True)
            
            #spt.show_plot(show_legend=True, line_id='C2_bands', convert_to_air=True)

            #spt.show_plot(show_legend=True, line_id='cyclotron3800', convert_to_air=True)
            #spt.show_plot(show_legend=True)
            #plt.legend()
            #plt.show()
            #spt.show_plot(show_legend=True, line_id='ca', convert_to_air=True)
            counter+=1
            #try:
                #plt.title(header['airoftyp'])
            #except KeyError:
                #pass
        #plt.xlim(3700,9000)
        #plt.ylabel(r'$f_{\nu}$ (arbitrary units)')
        #plt.legend(loc='best')
        #plt.show()
        #plot_telluric_spectrum([3700, 9000], smooth=False, pix_width=30)
        #plot_telluric_spectrum([3700, 9000], smooth=True, pix_width=30)
        #plot_telluric_spectrum([3700,9000], smooth=False, pix_width=30, tell_filename='LBL_A30_s0_w200_R0060000_T.fits')
        #spt.show_plot(show_legend=True, line_id='alkali', convert_to_air=True)
        
        #####Adding plotting of model spectra for Balmer lines
        
        #model_file_list=sorted(glob(cp.balmer_model_dir+'*.txt'))
        #for model_filename in model_file_list:
            #model_spec=spt.retrieve_blouin_model(model_filename)
            #trim_model_spec=spt.clean_spectrum(model_spec,3600,7000,[])
            #norm_model_spec=norm_spectrum(trim_model_spec,norm_range)
            #label=model_filename.split('/')[-1]
            #plt.plot(norm_model_spec[0],norm_model_spec[1],label=model_filename.split('/')[-1])
        
        ###### End of model spectra plotting
        print('\n\nSNRs')
        print(snr_list)
        print('\n\n')
        #print('ew_list:',ew_list)
        plt.axhline(y=0, linestyle=':', color='k')
        print('\n\nstandard deviations',sig_list)
        #spt.plot_wien_max_wave(np.arange(4500,6600, 500))
        #spt.show_plot(show_legend=False, line_id='h', convert_to_air=True, actually_show=False)
        #spt.show_plot(show_legend=True, line_id='cool_wd', convert_to_air=True)
        plt.ylim(-0.1,0.6)
        spt.show_plot(show_legend=True, line_id='h', convert_to_air=True,actually_show=True)
        #spt.show_plot(show_legend=True, line_id='C2_bands', convert_to_air=True)
        #spt.show_plot(show_legend=True, line_id='cyclotron3800', convert_to_air=True)
        #spt.show_plot(show_telluric=False, show_legend=True, line_id='')
        #spt.show_plot(show_telluric=False, show_legend=True, line_id='J0212', convert_to_air=True)
        spt.show_plot(show_legend=True, line_id='cool_wd', convert_to_air=True)
        
        plt.plot(bmjd_list,ew_list, marker='o')
        plt.ylabel(r'ew_values ($\AA$)')
        plt.show()
        
        plt.errorbar(bmjd_list,ew_list, yerr=ew_noise_list, marker='o')
        plt.ylabel(r'ew_values ($\AA$)')
        #plt.yscale('log')
        #plt.xscale('log')
        plt.show()
        
        #plt.plot(np.array(ew_list)/np.array(ew_noise_list),marker='o')
        #plt.ylabel(r'ew/ew_noise')
        #plt.show()
        #plt.ylabel('Integrated Flux (10^-16 erg/cm^2/s)')
        #plt.xlabel('BMJD_TDB')
        #plt.show()
        #spt.show_plot(show_legend=True, show_telluric=True)
        #spt.show_plot(show_telluric=True, )
        #plt.legend()
        #plt.show()
        plot_pix_shifts(filenames)
        #plot_rot_time(filenames)
        #plot_head_2_head(filenames,'rotator','airmass')
        #plot_head_2_head(filenames,'ROTATOR','POSANGLE')
        #plot_head_2_head(filenames,'BMJD_TDB','POSANGLE')
        plot_white_lightcurve(filenames,header_name='airofavg')
        plot_head_2_head(filenames,'seeing','see_fwhm')
        plot_head_2_head(filenames,'BMJD_TDB','SEE_FWHM')
        plot_head_2_head(filenames,'BMJD_TDB','seeing')
        plot_head_2_head(filenames,'BMJD_TDB','airmass')
        plot_head_2_head(filenames,'BMJD_TDB','rotator')
        plot_head_2_head(filenames,'BMJD_TDB','cam_ang')
        plot_head_2_head(filenames,'rotator','cam_ang')
        plot_head_2_head(filenames,'rotator','airofavg')
        #plot_head_2_head(filenames,'BMJD_TDB','airmass')
        #plot_head_2_head(filenames,'BMJD_TDB','airofavg')
        #plot_head_2_head(filenames,'BMJD_TDB','seeing')
        #plot_head_2_head(filenames,'BMJD_TDB','POSANGLE')
        print('r-band approximation')
        plt.title('r-band approximation')
        plot_white_lightcurve(filenames, wave_range=[5500., 6700.])
        
        plot_white_lightcurve(filenames)
        plot_white_lightcurve(filenames,header_name='AIRMASS')
        plot_white_lightcurve(filenames,header_name='rotator')
        #plot_coordinates(filenames)
        #plot_white_lightcurve(filenames,header_name='mount_el')
        #plot_head_2_head(filenames, 'mount_el', 'delta_atm_refraction', offset=5.)
        #plot_head_2_head(filenames, 'mount_el', 'delta_atm_refraction', offset=-5.)
        #plot_white_lightcurve(filenames,header_name='rotator')
        #plot_white_lightcurve(filenames,header_name='cam_ang')
        plot_dheaddhead_2_head(filenames,'BMJD_TDB','rotator')
        plot_dheaddhead_2_head(filenames,'BMJD_TDB','airofavg')
        #plot_dheaddhead_2_head(filenames,'airofavg','rotator')
        #plot_dheaddhead_2_head(filenames,'rotator','airofavg',dheader2='BMJD_TDB')
        #plot_dheaddhead_2_dheaddhead(filenames, 'rotator', 'BMJD_TDB', 'airofavg', 'BMJD_TDB')
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

