"""
Created by Ben Kaiser (UNC-Chapel Hill) (date not known of original creation.)

This is the file that actually applies the sensitivity curve to the target spectra. The sensitivity curves need to
already be produced by flux_calibration.py. The barycentric velocity correction actually needs to occur after the
flux-calibration step as the measured wavelengths are the ones that correspond to the given transmissivity of 
the atmosphere and instrument.

You have to manually (-ish) generate the 'listWCTB' file that includes the target spectra in the first column 
(whose prefixes are wctb.*) and then in the column beside each one, you should insert the flux standard 
spectrum that was produced by actual_reduction.py for that target.

Step 6  in Reduction

PART OF STEP-BY-STEP REDUCTION.
"""


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
from astropy import constants as const

import spec_plot_tools as spt

speclistname = "listWCTB"

color_list = ['k', 'r', 'g', 'b', 'm', 'cyan', 'purple']

speclist = np.genfromtxt(speclistname, dtype = 'str')
speclist = speclist.T
target_list = speclist[0]
sens_curve_list = speclist[1]

norm_range=[7470, 7530]
tell_range= [7430, 7800]
bad_noise_sub = 100
do_tell_corr= True
do_rv_barycorr=False
#do_ext_corr= True
default_ext_corr=True
plot_across_night=False
do_tell_waves=False

if do_tell_waves:
    from plot_spec import plot_telluric_spectrum, plot_spectrum

else:
    pass

#parkes_location = coords.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
#cerro_pachon_location = coords.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)

#def barycentric_vel_corr(header, wavelengths):
    #input_year = header['OPENDATE'] #gps-synched date
    #input_hours = header['OPENTIME'] #gps-synched time
    #exp_time= header['EXPTIME']*u.s
    #input_times = input_year+'T'+input_hours #formatting correctly
    #obs_time = Time(input_times, format = 'isot', scale = 'utc')
    #obs_time= obs_time+exp_time/2.
    #ra = header['RA']
    #dec = header['DEC']
    #radec = coords.SkyCoord(ra, dec, frame = 'icrs', unit= (u.hourangle, u.deg))
    #bary_corr = radec.radial_velocity_correction(obstime= obs_time, location = cerro_pachon_location)
    #bary_corr = bary_corr.to(u.km/u.s)
    #lambda_rest = (wavelengths*(u.Angstrom))*const.c.to(u.km/u.s)/(-1*bary_corr+const.c.to(u.km/u.s))
    #lambda_rest = lambda_rest.value
    #return lambda_rest

#def bad_noise_vals(noise_spec):
    ##bad_inds = np.where(noise_spec < 0)
    ##noise_spec[bad_inds]= bad_noise_sub
    #noise_spec = np.abs(noise_spec)
    #bad_inds = np.where(noise_spec > bad_noise_sub)
    #noise_spec[bad_inds] = bad_noise_sub
    #return noise_spec

count = 0
summed_flux = []
summed_counts =[]
times = []
airmasses= []
model_wavelength=7607 #that's the minimum of the model band
for target_file, sens_curve_file in zip(target_list, sens_curve_list):
    sens_curve_coeffs = np.genfromtxt(sens_curve_file)
    i= fits.open(target_file)
    header = fits.getheader(target_file)
    
    if 'goodman' in header['INSTRUME'].lower():
        instrument='goodman'
        do_ext_corr=default_ext_corr
    elif 'gmos-n' in header['INSTRUME'].lower():
        instrument='gmos-n'
        do_ext_corr=False
    else:
        print("No instrument recognized for deciding extinction correction.")
        print('header["INSTRUME"]', header['INSTRUME'])
    exptime = header['EXPTIME']
    wavelengths= i[0].data
    #counts = i[1].data/np.float_(exptime) #need counts/second
    #bkg_counts = i[2].data/np.float_(exptime) #need counts/second
    #noise_spec = i[3].data #don't need to divide this by the exposure time since it's normalized already in proportion to whatever units we use.
    counts = i[1].data
    bkg_counts = i[2].data
    noise_spec = i[3].data #don't need to divide this by the exposure time since it's normalized already in proportion to whatever units we use.
    dlambda=i[4].data
    if do_tell_waves:
        print('\n\n\nCorrect the wavelengths based on the telluric features.\nor change do_tell_waves=False\n\n\n')
        tell_name= glob('telluric_thru*'+ sens_curve_file)[0]
        print('tell_name', tell_name)
        tell_corr_spec= np.genfromtxt(tell_name, skip_header=1).T
        #model_wavelength=7607 #that's the minimum of the model band
        plot_spectrum([wavelengths, counts], target_file, header, smooth=True, norm=True, pix_width=5, kernel_type='box', norm_range=norm_range)
        #plot_telluric_spectrum([3700,9000], smooth=True, pix_width=30)
        plt.plot(tell_corr_spec[0], tell_corr_spec[1], label='Telluric residuals')
        plt.axvline(x=model_wavelength, color='k', linestyle='--')
        plt.xlim(tell_range)
        plt.ylim(0,1.5)
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.legend()
        plt.show()
        if count==0:
            model_wavelength = float(raw_input("Model spec wavelength>>>"))
        obs_wavelength= float(raw_input("Observed spec wavelength ("+str(model_wavelength)+ " for no change)>>>"))
        wave_offset=model_wavelength-obs_wavelength
        wavelengths=wavelengths+wave_offset
        header.append(card=('tell_off', wave_offset, 'wavelength offset by tellurics'))
    else:
        header.append(card=('tell_off', False, 'wavelength offset by tellurics'))
        pass
    #noise_spec = bad_noise_vals(noise_spec) #remove negative noise values and exceedingly high ones
    sens_curve = np.polyval(sens_curve_coeffs,wavelengths)
    #plt.plot(sens_curve)
    #plt.show()
    obs_spec= np.vstack([wavelengths, counts])
    obs_spec= np.copy(spt.counts_to_flambda(obs_spec, dlambda))
    print('Observed spectrum in units of erg/s/cm^2/angstrom')
    if do_ext_corr:
        print("Doing atmospheric extinction correction.")
        #obs_spec= np.vstack([wavelengths, counts])
        obs_spec= np.copy(spt.correct_extinction(obs_spec, header, plot_all=False))
        #wavelengths=np.copy(obs_spec[0])
        #counts=np.copy(obs_spec[1])
        header.append(card=('ext_corr', True, 'atmospheric extinction correction'))
    else:
        header.append(card=('ext_corr', False, 'atmospheric extinction correction'))
    #flux = counts/sens_curve
    flux= np.copy(obs_spec[1])
    flux=flux/sens_curve# I had commented out the only flux calibration part of the flux calibration...
    if do_tell_corr:
        tell_name= glob('telluric_thru*'+ sens_curve_file)[0]
        print('tell_name', tell_name)
        tell_corr_spec= np.genfromtxt(tell_name, skip_header=1).T
        interp_tell_corr= np.interp(wavelengths, tell_corr_spec[0], tell_corr_spec[1])
        #plt.plot(tell_corr_spec[0], tell_corr_spec[1], label=tell_name)
        #plt.plot(wavelengths, interp_tell_corr, label='interp_'+tell_name)
        #plt.plot(wavelengths, flux/np.nanmean(flux), label='normalized '+target_file)
        #plt.xlim(np.nanmin(wavelengths), np.nanmax(wavelengths))
        #spt.show_plot()
        #plot_spectrum([wavelengths, flux], target_file, header, smooth=True, norm=True, pix_width=5, kernel_type='box', norm_range=norm_range)
        ##plot_telluric_spectrum([3700,9000], smooth=True, pix_width=30)
        #plt.plot(wavelengths, interp_tell_corr, label='Interp Telluric residuals')
        #plt.plot(tell_corr_spec[0], tell_corr_spec[1], label='Telluric residuals')
        #plt.axvline(x=model_wavelength, color='k', linestyle='--')
        #plt.xlim(tell_range)
        #plt.ylim(0,1.5)
        #plt.xlabel(r'Wavelength ($\AA$)')
        #plt.legend()
        #plt.show()
        flux= flux/interp_tell_corr
        header.append(card=('tellcorr', True, 'Telluric Corrections performed'))
    else:
        header.append(card=('tellcorr', False, 'Telluric Corrections performed'))
        pass
    total_flux = np.sum(flux)
    times.append([header['BMJD_TDB']])
    airmasses.append(header['AIRMASS'])
    summed_flux.append(total_flux)
    summed_counts.append(np.sum(counts))
    header.append(card = ('Senscurv', sens_curve_file, 'file used for flux calibration'))
    #header.append(card = ('Units', 'ergs/cm/cm/s/A 10**-16', 'Units for flux'))
    header['units']= 'ergs/cm/cm/s/A 10**-16'
    header.append(card = ('Wavlngth', 0, 'Angstroms extension for wavelengths'))
    header.append(card = ('Flux', 1, 'in flux units extension for target flux values'))
    header.append(card = ('Bkg', 2, 'in flux units extension for bkg flux values'))
    header.append(card = ('Noise', 3, 'unitless. Normalized'))
    header.append(card=('dlambda',4, 'delta wavelength value of each bin'))
    if do_rv_barycorr:
        header.append(card = ('barycorr', True, 'wavelengths corrected to barycenter'))
        wavelengths = spt.barycentric_vel_corr(header, wavelengths)
    else:
        header.append(card = ('barycorr', False, 'wavelengths corrected to barycenter'))
    
    target_file = 'f'+target_file
    bkg_flux = bkg_counts/sens_curve
    
    
    
    ##########################3
    #convert fluxes to 1e-16 values... meaning multiply by 1e16
    flux=flux*1e16
    bkg_flux=bkg_flux*1e16
    
    
    
    
    ####################
    hdu=fits.PrimaryHDU(wavelengths, header = header)
    hdu1= fits.ImageHDU(flux)
    hdu2 = fits.ImageHDU(bkg_flux)
    hdu3 = fits.ImageHDU(noise_spec)
    hdu4=fits.ImageHDU(dlambda)
    hdulist = fits.HDUList([hdu, hdu1, hdu2, hdu3,hdu4])
    hdulist.writeto(target_file, overwrite = True)
    if count%4 == 0:
        color = color_list[count%len(color_list)]
        if plot_across_night:
            plt.plot(wavelengths, flux, label = header['OPENTIME'], color= color)
        else:
            pass
        #plt.plot(wavelengths, noise_spec*flux, label = header['OPENTIME']+ ' Noise', color = color)
        #plt.plot(wavelengths, noise_spec, label = header['OPENTIME']+ ' Noise no times', color = color)
    count += 1
#plt.ylim([-10, 10])
plt.legend()
plt.show()
airmasses = np.array(airmasses)
scale_airmass = airmasses * np.mean(summed_flux)/np.mean(airmasses)

plt.title('overall brightness change... allegedly')
plt.plot(times, summed_flux)
plt.plot(times, scale_airmass, label = 'scaled airmass')
plt.legend()
plt.show()

scale_airmass = airmasses * np.mean(summed_counts)/np.mean(airmasses)

plt.title('overall counts change')
plt.plot(times, summed_counts)
plt.plot(times, scale_airmass, label = 'scaled airmass')
plt.legend()
plt.show()


plt.title('airmasses over time')
plt.scatter(times, airmasses)
plt.show()
    
    
    
    
