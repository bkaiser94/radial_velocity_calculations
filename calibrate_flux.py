"""
This is the file that actually applies the sensitivity curve to the target spectra. The sensitivity curves need to
already be produced by flux_calibration.py. The barycentric velocity correctionactually needs to occur after the
flux-calibration step as the measured wavelengths are the ones that correspond to the given transmissivity of the
atmosphere and instrument.

You have to manually (-ish) generate the 'listWCTB' file that includes the target spectra in the first column (whose 
prefixes are wctb.*) and then in the column beside each one, you should insert the flux standard spectrum that
was produced by actual_reduction.py for that target.

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


speclistname = "listWCTB"

color_list = ['k', 'r', 'g', 'b', 'm', 'cyan', 'purple']

speclist = np.genfromtxt(speclistname, dtype = 'str')
speclist = speclist.T
target_list = speclist[0]
sens_curve_list = speclist[1]

bad_noise_sub = 100


parkes_location = coords.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
cerro_pachon_location = coords.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)

def barycentric_vel_corr(header, wavelengths):
    ra = header['RA']
    dec = header['DEC']
    radec = coords.SkyCoord(ra, dec, frame = 'icrs', unit= (u.hourangle, u.deg))
    bary_corr = radec.radial_velocity_correction(obstime= Time(header['DATE-OBS'], format = 'isot', scale= 'utc'), location = cerro_pachon_location)
    bary_corr = bary_corr.to(u.km/u.s)
    lambda_rest = (wavelengths*(u.Angstrom))*const.c.to(u.km/u.s)/(-1*bary_corr+const.c.to(u.km/u.s))
    lambda_rest = lambda_rest.value
    return lambda_rest

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
for target_file, sens_curve_file in zip(target_list, sens_curve_list):
    sens_curve_coeffs = np.genfromtxt(sens_curve_file)
    i= fits.open(target_file)
    header = fits.getheader(target_file)
    exptime = header['EXPTIME']
    wavelengths= i[0].data
    counts = i[1].data/np.float_(exptime) #need counts/second
    bkg_counts = i[2].data/np.float_(exptime) #need counts/second
    noise_spec = i[3].data #don't need to divide this by the exposure time since it's normalized already in proportion to whatever units we use.
    #noise_spec = bad_noise_vals(noise_spec) #remove negative noise values and exceedingly high ones
    sens_curve = np.polyval(sens_curve_coeffs,wavelengths)
    flux = counts/sens_curve
    total_flux = np.sum(flux)
    times.append([header['BMJD_TDB']])
    airmasses.append(header['AIRMASS'])
    summed_flux.append(total_flux)
    summed_counts.append(np.sum(counts))
    header.append(card = ('Senscurv', sens_curve_file, 'file used for flux calibration'))
    header.append(card = ('Units', 'ergs/cm/cm/s/A 10**-16', 'Units for flux'))
    header.append(card = ('Wavlngth', 0, 'Angstroms extension for wavelengths'))
    header.append(card = ('Flux', 1, 'in flux units extension for target flux values'))
    header.append(card = ('Bkg', 2, 'in flux units extension for bkg flux values'))
    header.append(card = ('Noise', 3, 'unitless. Normalized'))
    header.append(card = ('barycorr', True, 'wavelengths corrected to barycenter'))
    target_file = 'f'+target_file
    wavelengths = barycentric_vel_corr(header, wavelengths)
    bkg_flux = bkg_counts/sens_curve
    hdu=fits.PrimaryHDU(wavelengths, header = header)
    hdu1= fits.ImageHDU(flux)
    hdu2 = fits.ImageHDU(bkg_flux)
    hdu3 = fits.ImageHDU(noise_spec)
    hdulist = fits.HDUList([hdu, hdu1, hdu2, hdu3])
    hdulist.writeto(target_file, overwrite = True)
    if count%4 == 0:
        color = color_list[count%len(color_list)]
        plt.plot(wavelengths, flux, label = header['OPENTIME'], color= color)
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
plt.plot(times, airmasses)
plt.show()
    
    
    
    
