"""
Created by Ben Kaiser (UNC-Chapel Hill) (date not known of original creation.)
This should produce the sensitivity curve that is needed for calibrate_flux.py.

Uses the outputs of actual_reduction.py.

Have to manually change targets that are being used... probably could automate that to actually load all of the
necessary options after at least changing targets... doesn't seem really worth it at this point.

Step 5 of Reduction
"""
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
import cosmics
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
import scipy.interpolate as scinterp
import os


import spec_plot_tools as spt
import cal_params as cp
import get_cal_params as gcp
import model_manipulation as mm

from plot_spec import plot_telluric_spectrum, norm_spectrum


cwd= os.getcwd()
sub_dir = cwd.split('/')[-1]


#poly_degree=3
#poly_degree= 5 #order of polynomials from before 20190506
#poly_degree=7
#poly_degree=7
#model_poly_degree= 5
#sens_fit_method='poly/poly' #400M2 method generally, gets sens curve by dividing the polynomial observed by polynomial model
#sens_fit_method='empirical' #400M1 method generally, gets sens curve by dividing the obs flux directly by the model flux and then fitting a polynomial
division_extra_deg = 2

norm_range=[6630,6690]

header_char=':'
header_delim= '\t'

#ext_corr=True #correct extinction
use_fnu=False


standard_directory= cp.standard_dir



#output_filename = "GD108_sensitivity_curve.txt"
#output_filename = 'Feige67_sensitivity_curve.txt'
#output_filename= 'LTT6248_sensitivity_curve.txt'
#output_filename= 'EG274_sensitivity_curve.txt'
#output_filename = 'LTT3218_sensitivity_curve.txt'

#standard_file = standard_directory+standard_file

##standard_name = "GD108"
#standard_name = 'Feige67'
#standard_name = 'LTT6248'
#standard_name='EG274'
#standard_name = 'GD153'
#standard_name= 'LTT3218'
#standard_name='Feige110'
#standard_name= 'LTT7987'
#standard_name='GD71'
standard_name='EG131'

observed_file='ravg_wctb.EG131_gemini_600B.fits'

##observed_file = "wcmtb.GD108930blue.fits"
##observed_file = 'wcmtb.feige67930blue.fits'
#observed_file  = 'wcmtb.ltt6248930blue.fits'
##observed_file = 'wcmtb.eg274930blue.fits'
##observed_file= 'wcmtb.ltt3218930blue.fits'

#observed_file='ravg_wctb.GD71_400m2.fits'
#observed_file='ravg_wctb.GD71_400m1.fits'

#observed_file='ravg_wctb.EG274_400m1_fix.fits'
#observed_file='ravg_wctb.LTT7987_400m1.fits'
#observed_file='ravg_wctb.LTT7987_400m2.fits'
#observed_file='ravg_wctb.EG274_400m2.fits'
#observed_file='avg_wctb.EG274_400m2.fits'


#observed_file = "wcmtb.GD108930blue.fits"
#observed_file = 'wctb.0232_feige67_930_blue.fits'
#observed_file  = 'wctb.0244_LTT6248_930_blue.fits'
#observed_file  = 'wctb.0269_LTT6248_930_blue_1arcsec.fits'
#observed_file = 'wctb.0272_eg274_930_blue.fits'
#observed_file='ravg_wctb.LTT3218_400m1_moretime.fits'
#observed_file= 'wcmtb.ltt3218930blue.fits'

#observed_file='avg_wctb.LTT7987_400m2.fits'
#observed_file='avg_wctb.LTT7987second_400m2.fits'
#observed_file='avg_wctb.LTT7987_400m1.fits'

#observed_file='avg_wctb.EG274_400m1.fits'
#observed_file='avg_wctb.EG274_400m2.fits'
#observed_file='avg_wctb.EG274_400m2_fix.fits'
#observed_file='avg_wctb.Feige110_400m2.fits'
#observed_file='avg_wctb.Feige110stand_400m2.fits'
#observed_file='avg_wctb.Feige110second_400m2.fits'
#observed_file='avg_wctb.Feige110_400m1.fits'

#observed_file='avg_wctb.Feige110_400m2.fits'

#observed_file='avg_wctb.GD153_400m2.fits'
#observed_file='avg_wctb.GD153_400m1.fits'

#observed_file='avg_wctb.EG274_400m1.fits'
#observed_file='avg_wctb.EG274_400m1_fix.fits'
#observed_file='avg_wctb.EG274_400m2.fits'
#observed_file='avg_wctb.eg274_930_blue.fits'
#observed_file='avg_wctb.eg274am104_400m2.fits'
#observed_file='avg_wctb.eg274am126_400m2.fits'
#observed_file='avg_wctb.eg274am176_400m2.fits'
#observed_file='avg_wctb.GD153_400m1.fits'
#observed_file='avg_wctb.EG274_400m1_am13.fits'
#observed_file='avg_wctb.EG274_400m2_am13.fits'
##############################
##############################

def get_star_info(starname):
    standard_dict= cp.standard_dict[starname.lower()]
    standard_dict['filename']=standard_directory+standard_dict['filename']
    return standard_dict

def get_output_header(header):
    airmass= header['AIRMASS']
    obs_time = header['OPENTIME']
    obs_date = header['OPENDATE']
    obs_time = obs_date+'T'+obs_time
    obs_time = Time(obs_time, format = 'isot', scale = 'utc').mjd
    output_header_list= ['Airmass'+header_char+str(airmass), 'MJD'+header_char+str(obs_time)]
    for in_header, out_header in zip(cp.in_headers, cp.out_headers):
        value= header[in_header]
        new_entry= out_header+header_char+str(value)
        output_header_list.append(new_entry)
    output_header= header_delim.join(output_header_list)
    return output_header

def degrade_model(model_vals, obs_vals, header):
    """
    Convolve model spectrum with the seeing and rebin it (using flux-conservative method) to the pixel-scale of
    the observation.
    
    INPUTS:
        model_vals - [wavelengths, fluxes, wavelength bin widths]
        obs_vals - [wavelengths, fluxes, wavelength bin widths]
        header - observation header, which should include the pixel width, slit width, seeing, etc.
        
    OUTPUTS:
        model_spec 
    
    """
    rebinned_model_spec= spt.rebin_generic_spec(model_vals[:2], model_vals[2], obs_vals[0], obs_vals[2])
    slit_width = spt.get_slit_width(header)
    rebinned_model_spec= mm.convolve_model_new(rebinned_model_spec, header, slit_width=slit_width)
    output_spec= rebinned_model_spec
    return output_spec

##############################








############################
##########################



core_name= observed_file.split('.')[1] #get the part of the filename that follows the first period and exclude the extension
output_filename='sens_curv_' +sub_dir+'_'+ core_name+'.txt'

if use_fnu:
    output_filename=core_name+'_fnu_sensitivity_curve.txt'
else:
    pass


obs_fits = fits.open(observed_file)
header = fits.getheader(observed_file)
obs_waves1= obs_fits[0].data
obs_flux1 = obs_fits[1].data
obs_dlambda= obs_fits[4].data
airmass = header['AIRMASS']
obs_time = header['OPENTIME']
obs_date = header['OPENDATE']
obs_time = obs_date+'T'+obs_time
obs_time = Time(obs_time, format = 'isot', scale = 'utc').mjd
exptime = header['EXPTIME']

obs_spec= np.vstack([obs_waves1, obs_flux1])
try:
    dlambda= obs_fits[4].data
    obs_spec= np.copy(spt.counts_to_flambda(obs_spec, dlambda))
    print('Observed spectrum in units of erg/s/cm^2/angstrom')
except IndexError:
    print('No dlambda extension in observed FITS file. You need to redo wave_cal.py to include that extension.')
    print('This file is most likely generated before 2019-07-16 when this change was implemented')
    sys.exit()
    
if use_fnu:
    obs_spec[1]=obs_spec[1]*1e16
    obs_spec= np.copy(spt.flambda_to_fnu(obs_spec, dlambda))
    print('Observed spectrum in units of 10**-28 erg/s/cm^2/Hz')
else:
    pass

if 'goodman' in header['INSTRUME'].lower():
    instrument='goodman'
    ext_corr=True
elif 'gmos-n' in header['INSTRUME'].lower():
    instrument='gmos-n'
    ext_corr=False
else:
    print("No instrument recognized for deciding extinction correction.")
    print('header["INSTRUME"]', header['INSTRUME'])


if ext_corr:
    obs_spec= spt.correct_extinction(obs_spec, header, plot_all=True)
else:
    pass
obs_waves1=obs_spec[0]
obs_flux1=obs_spec[1]



setup_dict= gcp.get_cal_params(header)
setup_name=setup_dict['setupname']
sens_fit_method= cp.flux_cal_dict['sens_fit_method'][setup_name]
poly_degree=cp.flux_cal_dict['obs_poly_degree'][setup_name]
model_poly_degree=cp.flux_cal_dict['model_poly_degree'][setup_name]

#standard_file = standard_directory+standard_file
standard_info = get_star_info(standard_name)
print(type(standard_info['balmer_masks']))
print(type(standard_info['other_masks']))
wavelength_masks=standard_info['balmer_masks']+standard_info['other_masks']
print("wavelength_masks:", wavelength_masks)

stand_array = np.genfromtxt(glob(standard_info['filename'])[0]).T
#output_filename= standard_dict['sens_filename']

stand_waves1 = stand_array[0]


stand_flux1 = stand_array[1]  #ergs/cm/cm/s/A (That's exactly how it's written in the README for X-shooter)


def rescale_flux(stand_flux1):
    """
    Make sure the flux is in the correct units for your purposes, which should be 1e-16
    
    """
    standard_type= standard_info['filename'].split('/')[-2]
    print('standard_type:', standard_type)
    if ((standard_type == 'xshooter_standards') or (standard_type=='gemini_north_standards')) :
        print('should be resetting fluxes')
        stand_flux1= stand_flux1*1e16 #converts to 10**-16 flux vals hopefully
    else:
        stand_flux1= stand_flux1
    return stand_flux1

if use_fnu:
    stand_flux1= rescale_flux(stand_flux1)
    stand_spec= np.vstack([stand_waves1, stand_flux1])
    stand_spec=spt.flambda_to_fnu(stand_spec)
    stand_flux1=stand_spec[1]
    stand_waves1=stand_spec[0]
else:
    pass



plt.title('model versus observed')

norm_stand_spec=norm_spectrum([stand_waves1, stand_flux1], norm_range)
norm_obs_spec=norm_spectrum([obs_waves1, obs_flux1], norm_range)
#plt.plot(stand_waves1, stand_flux1,label = 'model')
#plt.plot(obs_waves1, obs_flux1, label = 'observed')
plt.plot(norm_stand_spec[0], norm_stand_spec[1],label = 'model')
plt.plot(norm_obs_spec[0], norm_obs_spec[1], label = 'observed')
plot_telluric_spectrum([3700,9000], smooth=True, pix_width=30)
#plt.scatter(stand_waves1, stand_flux1,label = 'model')
#plt.scatter(obs_waves1, obs_flux1, label = 'observed', color='r')
#plt.legend()
#plt.show()
if use_fnu:
    plt.ylabel('Flux(10**-28 erg/s/cm^2/Hz)')
else:
    plt.ylabel('Flux(erg/s/cm^2/A)')
plt.xlabel('Wavelengths (Angstroms)')
spt.show_plot()




unshifted_waves= np.copy(stand_waves1)
#plt.plot(stand_waves1, stand_flux1, label='unshifted')
stand_waves1= spt.barycentric_vel_uncorr(header, stand_waves1)
#plt.plot(stand_waves1, stand_flux1, label='shifted to Earth')
#plt.title('With and without Barycentric RV added')
#plt.legend(loc='best')
#plt.show()

#plt.plot(unshifted_waves,unshifted_waves-stand_waves1)
#plt.title('shift in angstroms by angstrom')
#plt.show()

#plt.plot(np.roll(unshifted_waves, 1)[1:]-unshifted_waves[1:], label= 'unshifted delta lambda')
#plt.plot(np.roll(stand_waves1, 1)[1:]- stand_waves1[1:], label='shifted delta lambda')
#plt.legend()
#plt.show()


do_offset= bool(raw_input("Do you need to do a wavelength offset?(Enter nothing to skip; Enter anything to do it.)>>>"))
if do_offset:
    print("Enter the approximate wavelength for the same feature in the model and observed spectra for offset")
    model_wavelength = float(raw_input("Model spec wavelength>>>"))
    obs_wavelength= float(raw_input("Observed spec wavelength>>>"))
    #dotted_pixel=0
    #emission_pixel=0
    offset = model_wavelength-obs_wavelength
    obs_waves1=obs_waves1+offset
    plt.title('model versus observed')
    plt.plot(stand_waves1, stand_flux1/np.nanmedian(stand_flux1),label = 'model')
    plt.plot(obs_waves1, obs_flux1/np.nanmedian(obs_flux1), label = 'observed')
    plt.xlabel('Wavelengths (Angstroms)')
    if use_fnu:
        plt.ylabel('Flux(10**-28 erg/s/cm^2/Hz)')
    else:
        plt.ylabel('Flux(erg/s/cm^2/A)')
    #plt.ylabel('Flux(erg/s/cm^2/A)')
    #plt.legend()
    #plt.show()
    spt.show_plot()
else:
    print("Skipping offsetting")
    

#plt.plot(stand_waves1, stand_flux1, label='unbinned')


standard_type= standard_info['filename'].split('/')[-2]
print('standard_type:', standard_type)
if ((standard_type == 'xshooter_standards') or (standard_type=='gemini_north_standards')) :
    print('no bin widths provided for model because X-Shooter... or Gemini-North')
    model_bin_widths= np.copy(np.roll(stand_waves1, -1) - stand_waves1)
    print(stand_waves1.shape, stand_flux1.shape, model_bin_widths.shape)
    stand_waves1= stand_waves1[:-1]
    stand_flux1=stand_flux1[:-1]
    model_bin_widths= model_bin_widths[:-1] #need to remove the first two pixels because they're going to be weird... or mayb it's the last two...
    plt.plot(stand_waves1, model_bin_widths)
    plt.title('model bin widths')
    plt.show()
    
    
print('degrading model spec')
rebin_model_spec= degrade_model(np.vstack([stand_waves1, stand_flux1, model_bin_widths]),  np.vstack([obs_waves1, obs_flux1, obs_dlambda]), header)


#plt.plot(stand_waves1, stand_flux1, label='unbinned model')
plt.scatter(stand_waves1, stand_flux1, label='unbinned model')
plt.plot(rebin_model_spec[0], rebin_model_spec[1], label='rebinned model_spec')
plt.legend()
plt.show()

original_stand_waves=np.copy(stand_waves1)
original_stand_flux=np.copy(stand_flux1)

stand_waves1=rebin_model_spec[0]
stand_flux1= rebin_model_spec[1]


min_wave = np.nanmin(obs_waves1)
max_wave = np.nanmax(obs_waves1)


plt.title('interpolated model versus standard model')
#interp_model_flux = np.interp(obs_waves1, stand_waves1, stand_flux1) #
interpolator = scinterp.CubicSpline(stand_waves1, stand_flux1)
interp_model_flux= interpolator(obs_waves1)

other_interpolator=scinterp.CubicSpline(original_stand_waves,original_stand_flux)
other_interp_flux=other_interpolator(obs_waves1)
plt.plot(original_stand_waves,original_stand_flux, label='original model points', linestyle='None', marker='o',color='magenta')
plt.plot(obs_waves1, interp_model_flux, label = 'interpolated')
plt.plot(stand_waves1, stand_flux1, label = 'model')
plt.plot(obs_waves1, other_interp_flux, label='original model interpolated')
plt.legend(loc='best')
plt.show()

obs_spec = np.vstack([obs_waves1, obs_flux1])
#stand_spec= np.vstack([stand_waves1, stand_flux1])
stand_spec= np.vstack([obs_waves1, interp_model_flux ])


obs_spec = spt.clean_spectrum(obs_spec, min_wave, max_wave, wavelength_masks)
unmasked_stand_spec=spt.clean_spectrum(stand_spec, min_wave, max_wave,standard_info['balmer_masks'])
stand_spec =spt.clean_spectrum(stand_spec, min_wave, max_wave, wavelength_masks)
stand_waves = stand_spec[0]
stand_flux= stand_spec[1]
obs_waves = obs_spec[0]
obs_flux = obs_spec[1]



plt.title('model versus observed')
plt.plot(stand_waves, stand_flux,label = 'model')
plt.plot(obs_waves, obs_flux, label = 'observed')
#plt.legend()
#plt.show()
spt.show_plot()

#plt.title('wavelength values')
#plt.plot(obs_waves,np.ones(obs_waves.shape), label = 'observed', marker='o')
#plt.plot(stand_waves, np.zeros(stand_waves.shape), label = 'model', marker= 'o')
#plt.legend()
#plt.show()


#interp_obs_flux = np.interp(stand_waves, obs_waves, obs_flux)

#sens_curve_points = obs_flux/interp_model_flux
#sens_curve_fit = np.polyfit(obs_waves, sens_curve_points, 5)
#sens_curve_points= interp_obs_flux/stand_flux
#sens_curve_fit= np.polyfit(stand_waves, sens_curve_points,5)

obs_curve= np.polyfit(obs_waves, obs_flux, poly_degree)
#obs_curve= np.polyfit(obs_waves, obs_flux, cp.flux_cal_dict['obs_poly_degree'][setup_name])
#model_curve = np.polyfit(stand_waves, stand_flux, poly_degree)
#model_curve = np.polyfit(stand_waves, stand_flux, model_poly_degree)
model_curve = np.polyfit(unmasked_stand_spec[0], unmasked_stand_spec[1], model_poly_degree)
#model_curve = np.polyfit(unmasked_stand_spec[0], unmasked_stand_spec[1], cp.flux_cal_dict['model_poly_degree'][setup_name])

calc_waves=np.linspace(min_wave, max_wave,1000)
#sens_curve_points = np.polyval(obs_curve, stand_waves)/np.polyval(model_curve, stand_waves)
#sens_curve_points = np.polyval(obs_curve, obs_waves)/np.polyval(model_curve, obs_waves)
#sens_curve_points= np.polyval(obs_curve,calc_waves)/np.polyval(model_curve, calc_waves)
#sens_curve_fit= np.polyfit(stand_waves, sens_curve_points,5)
#sens_curve_fit= np.polyfit(obs_waves, sens_curve_points,poly_degree)
#sens_curve_fit= np.polyfit(calc_waves, sens_curve_points,poly_degree+2)
#sens_curve_fit= np.polyfit(calc_waves, sens_curve_points,cp.flux_cal_dict['obs_poly_degree'][setup_name]+division_extra_deg)

#sens_curve_fit=np.polyfit(obs_waves, obs_flux/stand_flux,poly_degree) #20190618


plt.plot(obs_waves1, obs_flux1, label= 'observed', marker= 'o', linestyle='none')
plt.plot(obs_waves, obs_flux, label= 'observed used', marker= 'o', linestyle='none')
plt.plot(obs_waves1, np.polyval(obs_curve, obs_waves1), label = 'curve')

#obs_interpolator=scinterp.CubicSpline(obs_waves,obs_flux)
#obs_interp_flux=obs_interpolator(obs_waves1)
#plt.plot(obs_waves1, obs_interp_flux, label='interpolated obs')

plt.legend(loc='best')
plt.show()

plt.plot(obs_waves1, interp_model_flux, label= 'model', marker= 'o', linestyle = 'none')
plt.plot(stand_waves, stand_flux, label= 'model used', marker= 'o', linestyle = 'none')
plt.plot(obs_waves1, np.polyval(model_curve, obs_waves1), label = 'curve')
plt.legend(loc='best')
plt.show()


plt.plot(obs_waves1, obs_flux1/interp_model_flux, label='all obs/model', marker='o', linestyle='none')

if sens_fit_method == 'poly/poly':
    sens_curve_points= np.polyval(obs_curve,calc_waves)/np.polyval(model_curve, calc_waves)
    sens_curve_fit= np.polyfit(calc_waves, sens_curve_points,poly_degree+2)
    plt.plot(calc_waves, sens_curve_points, label = 'used obs curve/model curve', marker = 'o', linestyle = 'none')
elif sens_fit_method=='empirical':
    sens_curve_points=obs_flux/stand_flux
    sens_curve_fit=np.polyfit(obs_waves, obs_flux/stand_flux,poly_degree) #20190618
    plt.plot(obs_waves, sens_curve_points, label = 'used obs curve/model curve', marker = 'o', linestyle = 'none')
elif sens_fit_method=='poly/interp':
    sens_curve_points=np.polyval(obs_curve,calc_waves)/other_interpolator(calc_waves)
    sens_curve_fit=np.polyfit(obs_waves, obs_flux/stand_flux,poly_degree) #20190618
    plt.plot(calc_waves, sens_curve_points, label = 'used obs curve/model curve', marker = 'o', linestyle = 'none')
else:
    print('no valid sens_fit_method selected, currently selected:')
    print(sens_fit_method)

poly_curve = np.polyval(sens_curve_fit,obs_waves1)

#plt.plot(stand_waves, sens_curve_points, label = 'data points', marker = 'o', linestyle = 'none')
#plt.plot(obs_waves1, obs_flux1/interp_model_flux, label='all obs/model', marker='o', linestyle='none')
#plt.plot(obs_waves, sens_curve_points, label = 'used obs/model', marker = 'o', linestyle = 'none')
#plt.plot(calc_waves, sens_curve_points, label = 'used obs curve/model curve', marker = 'o', linestyle = 'none')
#plt.plot(calc_waves, sens_curve_points, label = 'data points', marker = 'o', linestyle = 'none')
#plt.plot(obs_waves1, poly_curve, label= 'polynomial fit')
plt.plot(obs_waves1, poly_curve, label= 'polynomial fit')
plt.legend(loc='best')
plt.show()

#fcal_obs = obs_flux1/poly_curve
fcal_obs = obs_flux1/poly_curve


#plt.plot(obs_waves1, fcal_obs, label ='flux calibrated observation', marker= 'o', linestyle = 'none')
#plt.plot(stand_waves1, stand_flux1, label = 'model', marker = 'o', linestyle = 'none')
#plt.plot(obs_waves1, fcal_obs, label ='flux calibrated observation')
#plt.plot(stand_waves1, stand_flux1, label = 'model')
plt.plot(obs_waves1, fcal_obs, label ='flux calibrated observation')
#plt.plot(stand_waves1, stand_flux1, label = 'model')
plt.plot(obs_waves1, interp_model_flux, label='model')
plt.xlabel('wavelength ($\AA$)')
plt.ylabel('Flux (ergs/cm/cm/s/A 1e-16)')
#plt.legend()
#plt.show()
spt.show_plot()

def get_residuals(model_flux=interp_model_flux, plot_all = False):
    residuals= fcal_obs/model_flux
    if plot_all:
        plt.plot(obs_waves1, residuals)
        plt.xlabel('wavelength ($\AA$)')
        #plt.ylabel('Flux (ergs/cm/cm/s/A 1e-16)')
        plt.ylabel('F_obscal/ F_model')
        plt.title('residuals')
        plt.axhline(y=1, color='k')
        plt.xlim(np.nanmin(obs_waves1), np.nanmax(obs_waves1))
        #plt.show()
        spt.show_plot()
    return residuals

def limit_to_telluric(obs_waves1,residuals, plot_all= True):
    """
    Take the calculated residuals and limit them to the region of telluric absorption and set the value to be 1 everywhere else, so the remnants of other regions aren't messed up too.
    
    """
    hold_telluric_lines = cp.telluric_lines
    hold_standard_lines = standard_info['balmer_masks']
    
    hold_telluric_lines= spt.check_overlaps(hold_telluric_lines, hold_standard_lines, np.nanmax(obs_waves1), np.nanmin(obs_waves1))
    if setup_dict['setupname']=='Gemini':
        print('Gemini observation so need segment gaps')
        hold_gaps=cp.gem_gaps
        #hold_telluric_lines=spt.check_overlaps(hold_telluric_lines, hold_gaps, np.nanmax(obs_waves1), np.nanmin(obs_waves1))
    else:
        pass
    io_telluric_lines= spt.make_inside_out(hold_telluric_lines, np.nanmin(obs_waves1), np.nanmax(obs_waves1))
    #inverted masks of telluric lines, so that the things that aren't telluric lines are masked
    output_factors = np.copy(residuals)
    factor_spec= np.vstack([obs_waves1, output_factors])
    for mask in io_telluric_lines:
        factor_spec= spt.replace_range(factor_spec, mask, method='ones')
    
    if plot_all:
        plt.title('telluric transmission factors')
        plt.xlabel('Wavelength (Angstroms)')
        plt.plot(obs_waves1, residuals, label='Residual factors before')
        plt.plot(factor_spec[0], factor_spec[1], label='Telluric factors')
        plot_telluric_spectrum([3700,9000], smooth=True, pix_width=30)
        plt.xlim(np.nanmin(factor_spec[0]), np.nanmax(factor_spec[0]))
        spt.show_plot()
    return factor_spec


if cp.flux_cal_dict['sens_fit_method'][setup_name]=='poly/interp':
    print('using residuals of interpolation')
    residuals= get_residuals(model_flux=other_interp_flux, plot_all=True)
    telluric_factor_spec= limit_to_telluric(obs_waves1, residuals)
    plt.plot(obs_waves1, other_interp_flux, label='model')
else:
    residuals= get_residuals(plot_all=True)
    telluric_factor_spec= limit_to_telluric(obs_waves1, residuals)
    plt.plot(obs_waves1, interp_model_flux, label='model')


plt.title('telluric corrected spectrum')
#plt.plot(obs_waves1, fcal_obs/residuals, label='flux-calibrated observation/residuals')
plt.plot(obs_waves1, fcal_obs/telluric_factor_spec[1], label='flux-calibrated observation/residuals')
plt.xlabel('wavelength ($\AA$)')
if use_fnu:
    plt.ylabel('Flux(10**-28 erg/s/cm^2/Hz)')
else:
    plt.ylabel('Flux (ergs/cm/cm/s/A 1e-16)')
#plt.legend()
#plt.show()
spt.show_plot()

#np.savetxt(output_filename, sens_curve_fit, header = 'Airmass: ' +str(airmass) + '\tMJD: ' +str(obs_time))
#np.savetxt(standard_info['sens_filename'], sens_curve_fit, header = 'Airmass: ' +str(airmass) + '\tMJD: ' +str(obs_time))

#np.savetxt('residuals_' + standard_info['sens_filename'], residuals, header='Airmass: ' +str(airmass) + '\tMJD: ' +str(obs_time))

output_header= get_output_header(header)


#np.savetxt(output_filename, sens_curve_fit, header = 'Airmass: ' +str(airmass) + '\tMJD: ' +str(obs_time))

#np.savetxt('residuals_' +output_filename, residuals, header='Airmass: ' +str(airmass) + '\tMJD: ' +str(obs_time))

#np.savetxt('telluric_thru_'+output_filename, telluric_factor_spec.T, header='Airmass: ' +str(airmass) + '\tMJD: ' +str(obs_time))


np.savetxt(output_filename, sens_curve_fit, header =output_header)

np.savetxt('residuals_' +output_filename, residuals, header=output_header)

np.savetxt('telluric_thru_'+output_filename, telluric_factor_spec.T, header=output_header)

