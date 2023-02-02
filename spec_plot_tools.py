"""
Created by Ben Kaiser (UNC-Chapel Hill) (date not known of original creation.)


"""
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy import constants as const
from astropy.time import Time
from astropy.table import Table
import scipy.interpolate as scinterp
from astropy import coordinates as coords
import time

import ref_index


import balmer_line_ranges as blr
import cal_params as cp


percentile= 50

soar_area= np.pi*(cp.soar_diameter/2.)**2 #area of SOAR light-gathering in meters (approximately)
parkes_location = coords.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
cerro_pachon_location = coords.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)


def naninfmax(input_array):
    array=np.array(input_array)
    #print(array.shape)
    #print(array.T.shape)
    array=array.T[0,0]
    mask_array=np.copy(array[~np.isnan(array)])
    mask_array= mask_array[~np.isinf(mask_array)]
    return np.max(mask_array)

def naninfmin(input_array):
    array=np.array(input_array)
    array=array.T[0,0]
    mask_array=np.copy(array[~np.isnan(array)])
    mask_array= mask_array[~np.isinf(mask_array)]
    return np.min(mask_array)

def air_to_vac(wavelengths):
    
    """
    convert wavelengths in air to vacuum using ref_index
    
    INPUT: wavelengths - in angstroms
    
    OUTPUT: wavelengths - in angstroms
    
    """
    print('wavelengths.min: ', wavelengths.min())
    new_wavelengths= wavelengths /10.
    new_wavelengths = ref_index.air2vac(new_wavelengths)*10. #back to angstroms
    print('new_wavelengths.min():', new_wavelengths.min())
    return new_wavelengths

def vac_to_air(wavelengths):
    """
    convert wavelengths in air to vacuum using ref_index
    
    INPUT: wavelengths - in angstroms
    
    OUTPUT: wavelengths - in angstroms
    
    """
    print('wavelengths.min: ', wavelengths.min())
    new_wavelengths= wavelengths /10.
    new_wavelengths = ref_index.vac2air(new_wavelengths)*10. #back to angstroms
    print('new_wavelengths.min():', new_wavelengths.min())
    return new_wavelengths

def make_inside_out(input_list, min_val, max_val):
    """
    Takes a list of lists (usually something like a mask_list from the other scripts) and then makes the selected regions become the outer boundaries.
    
    I.E. make_inside_out([[10,12],[14,18]], 8, 30)
    returns
    [[8,10],[12,14],[18,30]]
    which can then be used as a mask_list to do other stuff
    
    from balmer_line_ranges.py, but I figured it should just be in here actually
    """
    new_list=[[]]
    for index in range(0,len(input_list)+1):
        if index==0:
            new_list.append([min_val,input_list[index][0]])
        elif index==len(input_list):
            new_list.append([input_list[index-1][1],max_val])
        else:
            new_list.append([input_list[index-1][1],input_list[index][0]])
        #print new_list
    new_list= new_list[1:]
    #print new_list
    return new_list
        
def check_overlaps(mask_list1, mask_list2, wave_max, wave_min):
    """
    Check mask_list1 for any ranges that overlap from mask_list2, and return only those masks from mask_list1
    that do not overlap with any of the masks from mask_list2
    
    Also, checks if it overlaps with the bounds of some wavelength range so you aren't running up against the 
    edges.
    
    Using the lessons learned from rebin_spec()
    """
    output_mask_list=[]
    for mask1 in mask_list1:
        overlap_checks= []
        for mask2 in mask_list2:
            dif = min(mask1[1],mask2[1])-max(mask1[0],mask2[0])
            if dif <= 0:
                overlap_checks.append(0)
            else:
                overlap_checks.append(1)
        if ((np.sum(overlap_checks)<0.1) and (mask1[0] > wave_min) and (mask1[1] < wave_max)):
            output_mask_list.append(mask1)
        else:
            print(mask1, ' mask overlapped with another mask')
    return output_mask_list

def get_doppler_shifted(wavelengths, radial_velocity):
    #print "doppler shifting by ", radial_velocity
    lambda_obs = wavelengths * (radial_velocity*u.km/u.s + const.c.to(u.km/u.s)) / const.c.to(u.km/u.s)
    return lambda_obs.value

def dopp_shift_list(input_list, radial_velocity):
    dopp_list = []
    for waves in input_list:
        shift_waves = get_doppler_shifted(waves, radial_velocity)
        dopp_list.append(shift_waves)
    return dopp_list

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

def replace_range(input_spec, bound_list, method='ones'):
    """
    Different method of masking the desired region without literally removing the datapoints; it should reassign
    the values corresponding to certain wavelengths
    
    """
    wave_array= input_spec[0]
    other_array= input_spec[1]
    lower_bound = bound_list[0]
    upper_bound= bound_list[1]
    inbound_mask= np.where((wave_array>=bound_list[0]) & (wave_array<=bound_list[1]))
    if method=='ones':
        input_spec[1][inbound_mask]= 1.
    else:
        pass
    return input_spec


def get_pixel_scale(header):
    """
    input header
    
    return the arseconds per pixel of the image whose header you have
    
    """
    binning= header['CCDSUM'].split(' ')
    xbinning= int(binning[0])
    ybinning= int(binning[1])
    if xbinning != ybinning:
        print('\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print('x-binning and y-binning are not equal!\nSo... be careful I guess?')
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n')
    xscale= xbinning*cp.goodman_unbinned_pixscale
    yscale= ybinning*cp.goodman_unbinned_pixscale
    return xscale, yscale

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


def get_med_val(input_spec, wave_range):
    sub_spec = trim_spec(input_spec, wave_range[0], wave_range[1])
    length = sub_spec[0].shape[0]
    #if length%2 == 0:
        ##even number
        #sub_spec = sub_spec[:, :-1] #trim off the last point
    
    #med_val = np.nanmedian(sub_spec, axis =1)
    #med_flux = np.nanmedian(sub_spec[1])
    med_flux = np.percentile(sub_spec[1], percentile)
    med_index = np.where(sub_spec[1] == med_flux)[0]
    try:
        med_wave = sub_spec[0, med_index][0]
    except IndexError:
        min_index = np.argmin(np.abs(med_flux-sub_spec[1]))
        med_wave = sub_spec[0][min_index]
    med_val =[med_wave, med_flux]
    #print med_val
    return med_val

def make_continuum(input_spec, continuum_list= []):
    waves= []
    flux = []
    for ranges in continuum_list:
        #print ranges
        try:
            new_vals = get_med_val(input_spec, ranges)
            waves.append(new_vals[0])
            flux.append(new_vals[1])
        except IndexError as error:
            #print "Continuum construction point outside spectrum, so ignoring: ", ranges
            #print error
            pass
    wave_array = np.array(waves)
    flux_array = np.array(flux)
    continuum_spec = np.vstack([wave_array, flux_array])
    #plt.plot(input_spec[0], input_spec[1], label = 'input_spec')
    #plt.plot(continuum_spec[0], continuum_spec[1], linestyle = 'none', marker = 'o', label = 'continuum')
    #plt.legend()
    #plt.show()
    return continuum_spec

def get_norm_polynomial(input_spec, continuum_list = [], poly_degree = 3, plot_all = False, radial_velocity=0):
    #continuum_spec = make_continuum(input_spec, continuum_list = continuum_list)
    masks= dopp_shift_list(blr.balmer_fit_ranges, radial_velocity)
    #masks= dopp_shift_list(blr.io_balmer_norm_ranges, radial_velocity)
    #masks= dopp_shift_list(blr.io_continuum_list, radial_velocity)
    #continuum_spec= clean_spectrum(input_spec, np.nanmin(input_spec[0]), np.nanmax(input_spec[0]),masks)
    continuum_spec= input_spec
    #continuum_spec= clean_spectrum(input_spec, np.nanmin(input_spec[0]), np.nanmin(input_spec[0]).max, blr.balmer_norm_masks )
    #continuum_spec= clean_spectrum(input_spec, np.nanmin(input_spec[0]), np.nanmax(input_spec[0]), blr.balmer_fit_ranges )
    #print continuum_spec.shape
    poly_coeffs= np.polyfit(continuum_spec[0], continuum_spec[1], poly_degree)
    if plot_all:
        plt.plot(input_spec[0], input_spec[1], label = 'input_spec')
        plt.plot(continuum_spec[0], continuum_spec[1], linestyle = 'none', marker = 'o', label = 'continuum', color = 'r')
        #plt.plot(continuum_spec[0], continuum_spec[1], label = 'continuum', color = 'r')
        plt.plot(input_spec[0], np.polyval(poly_coeffs, input_spec[0]), label = 'fit')
        plt.title(continuum_list[0])
        plt.legend()
        plt.show()
    else:
        pass
    return poly_coeffs

def poly_norm_spec(input_spec, continuum_list = [], poly_degree = 3, plot_all  = False, radial_velocity=0):
    poly_coeffs = get_norm_polynomial(input_spec, continuum_list = continuum_list, poly_degree = poly_degree, plot_all = plot_all, radial_velocity=radial_velocity)
    poly_vals = np.polyval(poly_coeffs, input_spec[0])
    input_spec[1]= np.float_(input_spec[1])/poly_vals
    return input_spec



def retrieve_spec(filename, scale_noise= True):
    """
    Input: filename for the target spectrum you want to get
    
    Output: Spectrum made of a 2xN numpy array, header of the fits file you loaded it from
    """
    #print filename
    i=fits.open(filename)
    header = fits.getheader(filename)
    file_waves= np.copy(i[0].data)
    file_flux = np.copy(i[1].data)
    #set flux=0 as a slightly non-zero value to protect against NaNs
    zero_fluxes= np.where(np.abs(file_flux)<1e-10)
    #print('zero flux indices:',zero_fluxes)
    #plt.plot(file_waves, file_flux, label='pre fixes', marker='o')
    try:
        #if zero_fluxes[0].shape[0] > 0:
        file_flux[zero_fluxes]=1e-10
        #else:
            #pass
    except IndexError:
        pass
    nan_fluxes= np.where(file_flux==np.nan)
    #print('nan flux indices:',nan_fluxes)
    try:
        #if nan_fluxes[0].shape[0]> 0:
        file_flux[nan_fluxes]=1e-10
        #else:
            #pass
    except IndexError:
        pass
    #plt.plot(file_waves, file_flux, label='post fixes')
    #plt.title('spec')
    #plt.legend()
    #plt.show()
    file_noise = np.copy(i[3].data)
    #plt.plot(file_waves, file_noise, label='pre fixes', marker='o')
    nan_sigma= np.where(file_noise==np.nan)
    #print('nan sigma', nan_sigma)
    try:
        #if nan_sigma[0].shape[0]> 0:
        file_noise[nan_sigma]=1e10
        #else:
            #pass
    except IndexError:
        pass
    
    file_spec = np.vstack([file_waves, file_flux])
    file_noise_spec = np.vstack([file_waves, file_noise])
    if scale_noise:
        file_noise_spec[1] = file_spec[1]*file_noise_spec[1]
    else:
        pass
    return file_spec, header, file_noise_spec


def retrieve_sdss_spec(filename,scale_noise=True, wave_medium= 'air'):
    """
    Input: filename for a spectrum from SDSS (so it has the SDSS headers and fits format)
    
    Output:    Output: Spectrum made of a 2xN numpy array, header of the fits file you loaded it from, and noise spectrum. With wavelengths converted to AIR
    
    I.E. the same output format as retrieve_spec(), so that this is effectively a way of bringing all of the spectra together in the same format that I'm already tooled to use.
    
    """
    spec_hdu= fits.open(filename)
    spec_array=spec_hdu[1].data
    waves= 10.**np.copy(spec_array['loglam'])
    flux= np.copy(spec_array['flux'])
    flux= flux/10. #convert from 10**-17 to 10**-16
    
    if wave_medium=='air':
        waves= vac_to_air(waves) #conversion to air wavelengths
    elif wave_medium=='vac':
        pass
    
    try:
        noise= np.copy(spec_array['PropErr'])
    except KeyError as error:
        print(error)
        try:
            print('trying "ivar"')
            noise=np.copy(spec_array['ivar'])
            noise=np.sqrt(1./noise) #square root of the 1/variance value provided by SDSS
        except KeyError as error:
            print(error)
            print('setting noise=1')
            noise=np.ones(waves.shape[0])
        #print('setting noise=1')
        #noise=np.ones(waves.shape[0])
    file_spec=np.vstack([waves, flux])
    if scale_noise:
        file_noise_spec= np.vstack([waves, noise])
    else:
        file_noise_spec=np.vstack([waves, noise/flux])
    header= spec_hdu[1].header
    return file_spec, header, file_noise_spec

def retrieve_model_spec(filename):
    all_array=np.genfromtxt(filename).T
    return all_array

def retrieve_telluric_model(filename, wave_range):
    filename=cp.tell_dir+filename
    hdu=fits.open(filename)
    tell_array=hdu[1].data
    waves=tell_array['lam']*10000.
    transmission=tell_array['trans']
    tell_spec=np.vstack([waves,transmission])
    tell_spec=clean_spectrum(tell_spec, wave_range[0],wave_range[1], [])
    return tell_spec

def rescale_spectrum(input_spec, reference_spec, scale_range):
    input_value = get_med_val(input_spec, scale_range)[1]
    reference_value = get_med_val(reference_spec, scale_range)[1]
    scale_factor = reference_value/np.float_(input_value)
    input_spec[1] = input_spec[1]*scale_factor
    return input_spec
    
    
    
def plot_telluric():
    for region in cp.telluric_lines:
        plt.axvspan(region[0],region[1], alpha=0.05, color='k')
    return
    

def retrieve_nist_list(nist_file):
    nist_table=Table.read(nist_file, format='ascii.csv')
    #print(ne_table['intens'])
    #ne_table=remove_stars(ne_table)
    #print(ne_table['intens'])
    #try:
        #good_inds= np.where(ne_table['use']>0)
        #ne_table=ne_table[good_inds]
    #except KeyError:
        #pass
    #for count, row in enumerate(ne_table['obs_wl_air(A)']):
        #print(count, row)
    for count, row in enumerate(nist_table):
        #print(count, row['intens'])
        try:
            row['intens']=int(row['intens'])
            nist_table[count]['intens']=float(row['intens'])
            #print('converted to int!', row['intens'])
        except ValueError:
            print('ValueError ^^^^')
    #ne_table['intens']=ne_table['intens'].astype(float) #changing this column that gets read as strings for whatever reason
    #nist_table.pprint()
    return nist_table

def plot_line_markers(nist_file, wavelength_key='obs_wl_vac(A)', convert_to_air=False, label_pos=1.5,show_label=True):
    nist_table= retrieve_nist_list(nist_file)
    for row in nist_table:
        #print(name+name2, type(name))
        try:
            if row['show_line']:
                show_line=True
            else:
                show_line= False
        except KeyError:
            show_line=False
        if show_line:
            try:
                air_name=row['element']+' '+str(row['sp_num'])+'-'+str(row[wavelength_key])[:4]
            except KeyError as error:
                print('KeyError:', error)
                air_name='? '+str(row[wavelength_key])[:4]
            if convert_to_air:
                plt.axvline(x=vac_to_air(row[wavelength_key]), linestyle='--', color=cp.line_color_dict[row['element']])
            else:
                plt.axvline(x=row[wavelength_key], linestyle='--', color=cp.line_color_dict[row['element']])
            #plt.text(row['obs_wl_air(A)'], np.nanmax(counts), air_name, color='g', rotation=90)
            #plt.text(row[wavelength_key], 1. , air_name, color=cp.line_color_dict[row['element']], rotation=90, transform=ax.transAxes)
            #plt.text(row[wavelength_key], 1. , air_name, color=cp.line_color_dict[row['element']], rotation=90)
            if show_label:
                plt.text(row[wavelength_key], label_pos, air_name, color=cp.line_color_dict[row['element']], rotation=270)
            else:
                pass
        else:
            pass
    return

def show_plot(show_telluric=True, show_legend=True, line_id='', convert_to_air=False, label_pos=1.5,actually_show=True,show_label=True):
    if show_legend:
        plt.legend(loc='best')
    else:
        pass
    if show_telluric:
        plot_telluric()
    else:
        pass
    if line_id !='':
        line_list_file= cp.line_list_dir+cp.line_id_dict[line_id]
        print(line_list_file)
        plot_line_markers(line_list_file, convert_to_air=convert_to_air, label_pos=label_pos,show_label=show_label)
    else:
        pass
    if actually_show:
        plt.show()
    else:
        pass
    return

def get_photon_energy(wavelengths):
    """
    input wavelengths that need energy values (ideally in angstroms) or else needs to be an astropy unitted 
    object
    
    returns energy for each photon at each wavelength value in erg, but not as an astropy quantity object
    """
    try:
        print(wavelengths.value)
    except AttributeError:
        print('Input wavelengths do not have units for photon energy calculation.\nAssuming units of angstroms.')
        wavelengths= wavelengths*u.angstrom
    energy=(const.c*const.h/ wavelengths).cgs #energy in ergs
    return energy.value

def counts_to_flambda(input_spec, dlambda):
    """
    Returns a 'spec' type array after being given an input spec and the associated width of the wavelength bins
    
    output units: erg/s/cm^2/A, assuming the input spectrum is in units of counts/s (photons/s), which it should be...
    """
    photon_energies= get_photon_energy(input_spec[0]) #erg (probably erg/s)
    #dlambda= dlambda #delta wavelengths in angstroms
    soar_area_cm2= (soar_area*(u.meter**2)).to(u.cm**2)
    soar_area_cm2= soar_area_cm2.value #making it not be an astropy quantity
    
    flambda= input_spec[1]*photon_energies/dlambda/soar_area_cm2
    print(flambda)
    output_spec= np.vstack([input_spec[0], flambda])
    return output_spec

def counts_to_fnu(input_spec, dlambda):
    
    
    return
    
    
def flambda_to_fnu(input_spec, dlambda=0.):
    """
    returns a 'spec' type array after being given an input spec and the associated width of the wavelength bins
    
    dlambda defaults to zero, which is essentially the infinitesimal version of this, which is safe for actual 
    spectroscopy....I'm pretty sure.
    
    output_spec is in units of 10^-28 erg/s/cm^2/Hz assuming the input was in units of 10^-16 erg/s/cm^2/A
    """
    flambda= np.copy(input_spec[1])
    flambda= flambda*u.erg/(u.cm**2)/u.s/u.angstrom
    waves= np.copy(input_spec[0])
    waves= waves*u.angstrom
    delta_lambda= np.copy(dlambda)*u.angstrom
    fnu= flambda*(waves**2-0.25*delta_lambda**2)/const.c
    print('fnu', fnu)
    fnu=fnu.to(u.erg/u.cm/u.cm/u.s/u.hertz)
    print('fnu', fnu)
    fnu=fnu*1e12 #reducing the decimal stuff and making the units match the description
    output_spec= np.vstack([waves.value, fnu.value])
    return output_spec


def correct_extinction(input_spec, header, plot_all=False):
    """
    Read in the extinction curves and correct a spectrum for atmospheric extinction
    
    The extinction curve is in mags/airmass
    
    """
    extinction_filename=cp.ref_dir+'extinction/Stritzinger_2005_extinction_curve.csv'
    extinction_table= Table.read(extinction_filename)
    extinction_interp= scinterp.CubicSpline(extinction_table['lambda'], extinction_table['extinction'])
    extinction_vals= extinction_interp(input_spec[0])
    corr_flux= input_spec[1]*10.**(0.4*extinction_vals*header['airmass'])
    if plot_all:
        plt.plot(extinction_table['lambda'], extinction_table['extinction'], label='Original')
        plt.plot(input_spec[0], extinction_vals, label='Interpolated')
        plt.title('Extinction curve comparison')
        plt.ylabel('Extinction (mags/airmass)')
        plt.xlabel('Wavelength (angstroms)')
        plt.legend()
        plt.show()
        
        plt.plot(extinction_table['lambda'], 10.**(0.4*extinction_table['extinction']*header['airmass']), label='Original')
        plt.plot(input_spec[0], 10.**(0.4*extinction_vals*header['airmass']), label='Interpolated')
        plt.title('Extinction curve comparison')
        plt.ylabel('Extinction Coefficients')
        plt.xlabel('Wavelength (angstroms)')
        plt.legend()
        plt.show()
        
        plt.plot(input_spec[0], input_spec[1], label='Original spectrum')
        plt.plot(input_spec[0], corr_flux, label='Extinction Corrected')
        plt.xlabel('Wavelength (angstroms)')
        #plt.ylabel(header['units'])
        plt.ylabel('Flux (erg/s/cm^2/Angstrom)')
        plt.title("Spectrum correction")
        plt.legend()
        plt.show()
        
    else:
        pass
    input_spec[1]= corr_flux
    return input_spec
    


def barycentric_vel_corr(header, wavelengths):
    """
    Correct wavelengths to the barycenter of the solar system
    """
    input_year = header['OPENDATE'] #gps-synched date
    input_hours = header['OPENTIME'] #gps-synched time
    exp_time= header['EXPTIME']*u.s
    input_times = input_year+'T'+input_hours #formatting correctly
    obs_time = Time(input_times, format = 'isot', scale = 'utc')
    obs_time= obs_time+exp_time/2.
    ra = header['RA']
    dec = header['DEC']
    radec = coords.SkyCoord(ra, dec, frame = 'icrs', unit= (u.hourangle, u.deg))
    bary_corr = radec.radial_velocity_correction(obstime= obs_time, location = cerro_pachon_location)
    bary_corr = bary_corr.to(u.km/u.s)
    lambda_rest = (wavelengths*(u.Angstrom))*const.c.to(u.km/u.s)/(-1*bary_corr+const.c.to(u.km/u.s))
    lambda_rest = lambda_rest.value
    return lambda_rest

def barycentric_vel_uncorr(header, wavelengths, sys_vel= 0.0):
    """
    Shift wavelengths from barycenter back to the way they would be seen from Earth
    """
    input_year = header['OPENDATE'] #gps-synched date
    input_hours = header['OPENTIME'] #gps-synched time
    exp_time= header['EXPTIME']*u.s
    input_times = input_year+'T'+input_hours #formatting correctly
    obs_time = Time(input_times, format = 'isot', scale = 'utc')
    obs_time= obs_time+exp_time/2.
    ra = header['RA']
    dec = header['DEC']
    radec = coords.SkyCoord(ra, dec, frame = 'icrs', unit= (u.hourangle, u.deg))
    bary_corr = radec.radial_velocity_correction(obstime= obs_time, location = cerro_pachon_location)
    bary_corr = -1. * bary_corr.to(u.km/u.s) #you need the negative of the correction to get the barycentric velocity value
    bary_corr=bary_corr+sys_vel*u.km/u.s
    print('barycentric velocity:', bary_corr)
    lambda_obs= wavelengths + bary_corr*wavelengths/const.c.to(u.km/u.s)
    return lambda_obs
    
    
    
    


def rebin_image(im_array, rebin_axis=1, rebin_num= 10, plot_all= False):
    """
    rebin a Goodman image (technically any 2-d array, but I'm assuming it's Goodman)
    
    Input: 2-d image array
    kwarg: rebin_axis (1 is dispersion axis, x)
                rebin_num: number of pixels to bin together
    
    outputs: 
    output_im: 2-d image array that is rebinned
    output_indices: indices of the output array because they'll be missing pixels

    
    """
    copy_im = np.copy(im_array) #just in case I do something that ends up messing with the input; I'm paranoid
    indices= np.indices(copy_im.shape) #going to need to keep the location of the pixels that we're messing with
    binned_im=[]
    binned_indices=[] #initialize list to append to, yes I'm doing a for-loop because I'm inefficient.
    ax_len= copy_im.shape[rebin_axis]
    low_edges= np.arange(0,ax_len, rebin_num)
    if plot_all:
        plt.title('unbinned image')
        plt.imshow(np.log10(copy_im), cmap='hot')
        plt.show()
    else:
        pass
    for bin_edge in low_edges:
        if rebin_axis==1:
            sub_im= np.copy(copy_im[:,bin_edge:bin_edge+rebin_num])
            sub_indices= np.copy(indices[:, :, bin_edge:bin_edge+rebin_num])
        elif rebin_axis==0:
            pass
        sub_line= np.nanmean(sub_im, axis=rebin_axis)
        sub_line_inds= np.nanmean(sub_indices, axis=rebin_axis+1)
        #print('sub_line_inds', sub_line_inds)
        #print('sub_line.shape', sub_line.shape)
        binned_im.append(sub_line)
        binned_indices.append(sub_line_inds[rebin_axis])
        #plt.title('rebinned element')
        #plt.imshow(np.log10(sub_line), cmap='hot')
        #plt.show()
    
    binned_im_array= np.array(binned_im)
    if rebin_axis==1:
        binned_im_array=binned_im_array.T
    else:
        pass
    print('binned_im_array.shape', binned_im_array.shape)
    binned_indices_array= np.array(binned_indices)
    print('binned_indices_array.shape', binned_indices_array.shape)
    #plt.imshow(np.log10(binned_im_array), cmap='hot')
    if plot_all:
        plt.title('rebinned' + str(rebin_num) +' pixels on axis' + str(rebin_axis))
        plt.imshow(binned_im_array, cmap='hot')
        plt.show()
    else:
        pass
    return binned_im_array, binned_indices_array


def discrete_int(input_val):
    """
    inputs:
        input_val - some float value (created to handle a single point from a polynomial output, but it needs to be 
            an input
    outputs:
        output_val - a numpy integer dtype that was first rounded to the nearest whole number
    """
    output_val= np.int_(np.around(input_val, decimals=0))
    return output_val
    
    
def get_edges(waves, dlambda):
    return waves-0.5*dlambda, waves+0.5*dlambda


def rebin_spec(input_filename, desired_waves, desired_dlambda):
    """
    Take a spectrum (with all of its extensions) and rebin it to some new wavelength bins
    
    All values need to be in units of f_lambda ... so the noise2 calculation is almost certainly wrong at the moment
    
    Created by Ben Kaiser (UNC-Chapel Hill)
    """
    target_spec, header, target_noise = retrieve_spec(input_filename)
    noise2= target_noise[1]**2
    hdu=fits.open(input_filename)
    sky=np.copy(hdu[2].data)
    dlambda_spec= np.copy(hdu[4].data)
    target_flux = target_spec[1]
    target_low_edges, target_high_edges= get_edges(target_spec[0],dlambda_spec)
    desired_low_edges, desired_high_edges= get_edges(desired_waves, desired_dlambda)
    rebin_flux_list= []
    rebin_sky_list=[]
    rebin_noise2_list=[]
    for des_low, des_high in zip(desired_low_edges, desired_high_edges):
        stretch_high = np.ones(target_high_edges.shape)*des_high
        stretch_low= np.ones(target_low_edges.shape)*des_low
        upper= np.nanmin([stretch_high, target_high_edges], axis=0)
        lower= np.nanmax([stretch_low, target_low_edges], axis=0)
        dif = upper-lower
        negatives = np.where(dif<=0)
        dif[negatives]=0.
        #rebin_factors = dif /dlambda_spec
        rebin_factors= dif/(des_high-des_low)
        
        rebin_flux_list.append(np.sum(target_flux*rebin_factors))
        rebin_sky_list.append(np.sum(sky* rebin_factors))
        rebin_noise2_list.append(np.sum(noise2*rebin_factors**2))
    
    ##plt.plot(desired_waves, rebin_flux_list, label='rebinned spectrum')
    ##plt.plot(target_spec[0], target_spec[1], label='original spectrum')
    ##plt.xlabel('Wavelength (Angstroms)')
    ##plt.ylabel(header['UNITS'])
    ##plt.legend()
    ##plt.show()
    ####### This is an attempt at masking out the values from the merging of spectra that cover different setups
    #rebin_flux_array=np.array(rebin_flux_list)
    #rebin_sky_array= np.array(rebin_sky_list)
    #rebin_noise2_array=np.array(rebin_noise2_list)
    #nan_wants= np.where(rebin_flux_array==0.)
    #print('nan_wants:',nan_wants)
    #rebin_flux_array[nan_wants]==np.nan
    #rebin_sky_array[nan_wants]==np.nan
    #rebin_noise2_array[nan_wants]=np.nan
    #return [desired_waves, rebin_flux_array, rebin_sky_array, rebin_noise2_array, desired_dlambda]
    ########^^^ That isn't ready yet #######
    return [desired_waves, np.array(rebin_flux_list), np.array(rebin_sky_list), np.array(rebin_noise2_list), desired_dlambda]
    
    
    
def rebin_generic_spec(input_spec, input_dlambda, desired_waves, desired_dlambda):
    """
    Takes any spectrum as an input and rebins it; this version is less efficient than rebin_spec() if you're doing it 
    on the reduced Goodman spectra because it wont' simultaneously rebin the sky and noise and stuff, which 
    the other one does. This function does, however, actually just work on whatever you give it.
    
    Created by Ben Kaiser (UNC-Chapel Hill)
    
    """
    target_low_edges, target_high_edges= get_edges(input_spec[0],input_dlambda)
    desired_low_edges, desired_high_edges= get_edges(desired_waves, desired_dlambda)
    rebin_flux_list= []
    for des_low, des_high in zip(desired_low_edges, desired_high_edges):
        stretch_high = np.ones(target_high_edges.shape)*des_high
        stretch_low= np.ones(target_low_edges.shape)*des_low
        upper= np.nanmin([stretch_high, target_high_edges], axis=0)
        lower= np.nanmax([stretch_low, target_low_edges], axis=0)
        dif = upper-lower
        negatives = np.where(dif<=0)
        dif[negatives]=0.
        #rebin_factors = dif /dlambda_spec
        rebin_factors= dif/(des_high-des_low)
        
        rebin_flux_list.append(np.sum(input_spec[1]*rebin_factors))
    output_spec= np.vstack([desired_waves, rebin_flux_list])
    return output_spec



def get_slit_width(header):
    if '3.2' in header['slit']:
        slit_width = 3.2
    elif '1.0' in header['slit']:
        slit_width= 1.0
    
    return slit_width


def get_generic_ew(input_spec, wave_range, noise=1e-10, dlams=[], cont_method= 'avg', cont_width=20, plot_all=False, noise_method='prop'):
    """
    This one should get the equivalent width on just an arbitrary spectrum instead of starting from a file
    
    The noise kwarg has what I assume is a relatively small value default entered (but not zero) so that you don't end up with Nan's but so that there's a built-in noise handling method.
    
    Assumes that input spectrum is in units of f_lambda
    
    """
    try:
        print(noise.shape)
        noise_spec=noise
    except AttributeError:
        if noise_method != 'rms':
            print("Probably should be using noise_method='rms' since you didn't provide a noise value")
        else:
            pass
        noise_spec=np.ones(input_spec.shape)
        noise_spec[1]= noise*noise_spec[1]
        noise_spec[0]=np.copy(input_spec[0])
    print('np.nanmean(noise_spec[1])', np.nanmean(noise_spec[1]))
    
    if dlams ==[]:
        #creating delta lambda values in case the spectrum doesn't actually have them by just taking the difference of centers.
        print("No dlams provided as a kwarg, so we'll make them up.")
        dlams= np.copy(input_spec[0][1:]-input_spec[0][:-1])
        dlams=np.append(dlams, dlams[-1])
        dlams_spec= np.vstack([np.copy(input_spec[0]),dlams])
    else:
        dlams_spec=dlams
    cont_spec= clean_spectrum(input_spec, wave_range[0]-cont_width, wave_range[1]+cont_width, [wave_range]) #make a spectrum of only the 'continuum' regions
    cont_noise= clean_spectrum(noise_spec, wave_range[0]-cont_width, wave_range[1]+cont_width, [wave_range])
    abs_spec= clean_spectrum(input_spec, wave_range[0], wave_range[1],[])
    abs_noise=clean_spectrum(noise_spec, wave_range[0], wave_range[1],[])
    abs_dlams= clean_spectrum(dlams_spec, wave_range[0], wave_range[1], [])
    
    
    #print('sum abs_spec/merged noise:', np.sum(abs_spec[1])/np.sqrt(np.sum(abs_noise[1]**2)))
    #plt.plot(abs_noise[0],abs_noise[1],label='og')
    if cont_method=='avg':
        #print('cont_spec.shape', cont_spec.shape)
        #print('abs_spec.shape', abs_spec.shape)
        cont_val_abs= np.nanmean(cont_spec[1])
        cont_val=cont_val_abs
        #print('mean of cont_noise:', np.nanmean(cont_noise[1]))
        #print('np.nanmean(noise)/sqrt(cont_noise.shape[1]))', np.nanmean(cont_noise[1])/np.sqrt(cont_noise.shape[1]))
        if noise_method=='prop':
            cont_noise_merge= np.sqrt(np.sum(cont_noise[1]**2)/cont_noise.shape[1]**2)
        elif noise_method=='rms':
            #cont_noise_single=np.std(cont_spec[1])*np.sqrt(cont_spec.shape[1]) #you have to deaden the sigma value to recover the  standard deviation as the combined noise value of the inferred noise level
            cont_noise_single=np.std(cont_spec[1]) #Maybe one shouldn't 
            #print('np.std(cont_spec[1])', np.std(cont_spec[1]))
            cont_noise_merge=np.copy(np.sqrt((cont_noise.shape[1]*cont_noise_single**2)/cont_noise.shape[1]**2))
            print('cont_noise_single', cont_noise_single, 'cont_noise_merge', cont_noise_merge)
            cont_noise_single=cont_noise_single/cont_val #rescaling the noise level to the flux now that I'm about to use it...
            abs_noise[1]=cont_noise_single*abs_spec[1] #rescaling that expected noise level per pixel to the flux levels in the different cells of interest
            print('abs_noise[1]',np.nanmean(abs_noise[1]))
            print('cont_spec.shape', cont_spec.shape)
            print('abs_spec.shape', abs_spec.shape)
    elif cont_method=='poly1':
        #print('cont_spec.shape', cont_spec.shape)
        #print('abs_spec.shape', abs_spec.shape)
        #cont_val= np.nanmean(cont_spec[1])
        poly_coeffs=np.polyfit(cont_spec[0], cont_spec[1], 1)
        cont_val_abs= np.polyval(poly_coeffs, abs_spec[0])
        cont_val=np.nanmean(cont_val_abs)
        #print('mean of cont_noise:', np.nanmean(cont_noise[1]))
        #print('np.nanmean(noise)/sqrt(cont_noise.shape[1]))', np.nanmean(cont_noise[1])/np.sqrt(cont_noise.shape[1]))
        if noise_method=='prop':
            cont_noise_merge= np.sqrt(np.sum(cont_noise[1]**2)/cont_noise.shape[1]**2)
        elif noise_method=='rms':
            #cont_noise_single=np.std(cont_spec[1])*np.sqrt(cont_spec.shape[1]) #you have to deaden the sigma value to recover the  standard deviation as the combined noise value of the inferred noise level
            
            #cont_noise_rms= np.sqrt(np.nanmean((cont_spec[1]-np.polyval(poly_coeffs, cont_spec[0])) **2 ))
            cont_noise_rms=np.sqrt(np.sum((cont_spec[1]-np.polyval(poly_coeffs, cont_spec[0]))**2)/(cont_spec.shape[1]-2)) #minus 2 because there's supposed to be a minus 1 in the first place, and the sloap is an additional free parameter in this case
            cont_noise_single=np.copy(cont_noise_rms)
            
            
            cont_noise_merge= cont_noise_rms
            #print('np.std(cont_spec[1])', np.std(cont_spec[1]))
            #cont_noise_single= cont_noise_rms*np.sqrt(np.float_(cont_spec.shape[1]))
            #cont_noise_mean=np.copy(np.sqrt((cont_noise.shape[1]*cont_noise_single**2)/cont_noise.shape[1]**2))
            print('cont_noise_single', cont_noise_single, 'cont_noise_rms', cont_noise_merge)
            #cont_noise_single=cont_noise_single/cont_val #rescaling the noise level to the flux now that I'm about to use it...
            print('cont_noise_single', cont_noise_single)
            abs_noise[1]=cont_noise_single
            #abs_noise[1]=cont_noise_single*abs_spec[1] #rescaling that expected noise level per pixel to the flux levels in the different cells of interest
            print('abs_noise[1]',np.nanmean(abs_noise[1]))
            print('cont_spec.shape', cont_spec.shape)
            print('abs_spec.shape', abs_spec.shape)
        elif noise_method=='rms_wrong':
            #cont_noise_single=np.std(cont_spec[1])*np.sqrt(cont_spec.shape[1]) #you have to deaden the sigma value to recover the  standard deviation as the combined noise value of the inferred noise level
            print("\n=======")
            print("wrong method (I think) of RMS error calculation used here, where you don't rescale the noise value by the sqrt of the number of measurements that went into it to get the single value.")
            cont_noise_rms= np.sqrt(np.nanmean((cont_spec[1]-np.polyval(poly_coeffs, cont_spec[0]) )**2 ))
            cont_noise_merge= cont_noise_rms
            #print('np.std(cont_spec[1])', np.std(cont_spec[1]))
            cont_noise_single= cont_noise_rms
            #cont_noise_mean=np.copy(np.sqrt((cont_noise.shape[1]*cont_noise_single**2)/cont_noise.shape[1]**2))
            print('cont_noise_single', cont_noise_single, 'cont_noise_rms', cont_noise_merge, cont_noise_rms*np.sqrt(np.float_(cont_spec.shape[1])))
            cont_noise_single=cont_noise_single/cont_val #rescaling the noise level to the flux now that I'm about to use it...
            print('cont_noise_single', cont_noise_single)
            abs_noise[1]=cont_noise_single*abs_spec[1] #rescaling that expected noise level per pixel to the flux levels in the different cells of interest
            print('abs_noise[1]',np.nanmean(abs_noise[1]))
            print('cont_spec.shape', cont_spec.shape)
            print('abs_spec.shape', abs_spec.shape)
        else:
            pass
        #plt.plot(abs_noise[0],abs_noise[1],label='rms')
        #plt.legend(loc='best')
        #plt.show()
        #print('combined noise:', cont_noise_mean)
        #print('cont std:', np.std(cont_spec[1]))
    cont_energy= np.sum(abs_dlams[1]*cont_val_abs)
    #cont_energy_noise=np.sqrt(np.sum((abs_dlams[1]*cont_noise_mean)**2))
    #cont_energy_noise=np.sqrt(np.sum((abs_dlams[1]*cont_noise_merge)**2))
    cont_energy_noise=np.sqrt(np.sum((abs_dlams[1]*cont_noise_merge)**2))
    #print('divisions:', cont_val_abs/cont_noise_merge, cont_energy/cont_energy_noise)
    print('cont_energy', cont_energy, '+/-', cont_energy_noise)
    #else:
        #print("no valid 'cont_method' specified")
        #pass
    
    abs_energy= np.sum(abs_spec[1]*abs_dlams[1])
    abs_energy_noise= np.sqrt(np.sum((abs_noise[1]*abs_dlams[1])**2))
    print('abs_energy/abs_energy_noise:', abs_energy/abs_energy_noise)
    print("abs_energy", abs_energy, '+/-', abs_energy_noise)
    total_noise= np.sqrt(abs_energy_noise**2+cont_energy_noise**2)
    energy_dif= cont_energy-abs_energy
    print("energy_dif", energy_dif, '+/-', total_noise)
    #mid_index= int(abs_spec.shape[1]/2.)
    #used_dlam= abs_dlams[1][mid_index]
    ew_noise= total_noise/cont_val
    ew=energy_dif/cont_val
    if plot_all:
        plt.plot(abs_spec[0], cont_val_abs*np.ones(abs_spec[0].shape), label='continuum used')
        plt.plot(cont_spec[0], cont_spec[1], label='continuum source')
        plt.plot(abs_spec[0], abs_spec[1], label='absorption')
        #plt.errorbar(abs_spec[0],abs_spec[1],yerr=abs_noise[1], label='absorption')
        plt.axvspan(np.mean(wave_range)-0.5*ew, np.mean(wave_range)+0.5*ew, alpha=0.1, color='r',label='EW')
        plt.legend(loc='best')
        plt.show()
    else:
        pass
    print('EW:', ew, '+/-', ew_noise)
    
    return ew, ew_noise
        
def get_ew(filename, wave_range, cont_method= 'avg', cont_width=20, plot_all=False, noise_method='prop'):
    """
    Take a FITS file from the Goodman reduction process (most likely a ravg_fwctb file)
    and then get the equivalent width for some portion of the spectrum.
    
    assumes input spectrum is in units of f_lambda
    """
    target_spec, header, target_noise= retrieve_spec(filename)
    print('\n=====\n')
    i=fits.open(filename)
    dlams= np.copy(i[4].data)
    dlams_spec=np.vstack([np.copy(target_spec[0]), dlams])
    cont_spec= clean_spectrum(target_spec, wave_range[0]-cont_width, wave_range[1]+cont_width, [wave_range]) #make a spectrum of only the 'continuum' regions
    cont_noise= clean_spectrum(target_noise, wave_range[0]-cont_width, wave_range[1]+cont_width, [wave_range]) #make a spectrum of only the 'continuum' regions
    #if plot_all:
        #plt.plot(target_noise[0], target_spec[1]/target_noise[1])
        #plt.ylabel('S/N')
        #plt.xlim(wave_range[0]-cont_width, wave_range[1]+cont_width)
        #plt.show()
    #else:
        #pass
    ew, ew_noise=get_generic_ew(target_spec, wave_range, cont_method=cont_method, noise=target_noise, cont_width=cont_width, plot_all=plot_all, noise_method=noise_method, dlams=dlams_spec)
    #abs_spec= clean_spectrum(target_spec, wave_range[0], wave_range[1],[])
    #abs_noise=clean_spectrum(target_noise, wave_range[0], wave_range[1],[])
    #abs_dlams= clean_spectrum(dlams_spec, wave_range[0], wave_range[1], [])
    ##print('sum abs_spec/merged noise:', np.sum(abs_spec[1])/np.sqrt(np.sum(abs_noise[1]**2)))
    
    #if cont_method=='avg':
        ##print('cont_spec.shape', cont_spec.shape)
        ##print('abs_spec.shape', abs_spec.shape)
        #cont_val= np.nanmean(cont_spec[1])
        ##print('mean of cont_noise:', np.nanmean(cont_noise[1]))
        ##print('np.nanmean(noise)/sqrt(cont_noise.shape[1]))', np.nanmean(cont_noise[1])/np.sqrt(cont_noise.shape[1]))
        ##cont_noise_mean= np.sqrt(np.sum(cont_noise[1]**2)/cont_noise.shape[1]**2)
        ##print('combined noise:', cont_noise_mean)
        ##print('cont std:', np.std(cont_spec[1]))
        ##cont_energy= np.sum(abs_dlams[1]*cont_val)
        ##cont_energy_noise=np.sum(abs_dlams[1]*cont_noise_mean)
        ##cont_energy_noise=np.sqrt(np.sum((abs_dlams[1]*cont_noise_mean)**2))
        ##print('divisions:', cont_val/cont_noise_mean, cont_energy/cont_energy_noise)
        ##print('cont_energy', cont_energy, '+/-', cont_energy_noise)
        #if noise_method=='prop':
            #cont_noise_mean= np.sqrt(np.sum(cont_noise[1]**2)/cont_noise.shape[1]**2)
        #elif noise_method=='rms':
            #cont_noise_single=np.std(cont_spec[1])*np.sqrt(cont_spec.shape[1])
            #cont_noise_mean=np.sqrt((cont_noise.shape[1]*cont_noise_single**2)/cont_noise.shape[1]**2)
            #abs_noise[1]=cont_noise_single
        #cont_energy= np.sum(abs_dlams[1]*cont_val)
        #cont_energy_noise=np.sum(abs_dlams[1]*cont_noise_mean)
        #if plot_all:
            #plt.plot(abs_spec[0], cont_val*np.ones(abs_spec[0].shape), label='continuum used')
        #else:
            #pass
    #else:
        #print("no valid 'cont_method' specified")
        #pass
    #if plot_all:
        #plt.plot(cont_spec[0], cont_spec[1], label='continuum source')
        #plt.plot(abs_spec[0], abs_spec[1], label='absorption')
        #plt.legend(loc='best')
        #plt.show()
    #else:
        #pass
    #abs_energy= np.sum(abs_spec[1]*abs_dlams[1])
    #abs_energy_noise= np.sqrt(np.sum((abs_noise[1]*abs_dlams[1])**2))
    #print('abs_energy/abs_energy_noise:', abs_energy/abs_energy_noise)
    #print("abs_energy", abs_energy, '+/-', abs_energy_noise)
    #total_noise= np.sqrt(abs_energy_noise**2+cont_energy_noise**2)
    #energy_dif= cont_energy-abs_energy
    #print("energy_dif", energy_dif, '+/-', total_noise)
    ##mid_index= int(abs_spec.shape[1]/2.)
    ##used_dlam= abs_dlams[1][mid_index]
    #ew_noise= total_noise/cont_val
    #ew=energy_dif/cont_val
    #print('EW:', ew, '+/-', ew_noise)
   
    return ew, ew_noise
    
    
def initiate_science_plot():
    """
    Create all of the conditions required for a Science figure. I'm going to be calling this a ton of times anyway, so I figured I should get it right now.
    
    """
    #import matplotlib
    #matplotlib.use('pdf')
    #import importlib
    #importlib.reload(plt)
    plt.rc('font',family='Microsoft Sans Serif',size=9)


    
    return

def start_ApJ_fig(width_cols=1, constrained_layout=True, width_height=[1.,1.]):
    """
    width_cols: the width of the figure in column widths for 2-column ApJ layout. Default is the width of a single column in the 2-column format.. Must be an integer! (and only 1 or 2)
    
    
    width_height: is a list of the width and height dimensions of the desired figure that you already have. It is only the ratio of these two values that matters as they are rescaled by the desired width in units of columns as specified by width_cols
    """
    
    single_col_width=3.3522420091324205
    double_col_width=7.100005949910059
    if width_cols==1:
        new_width=single_col_width
    elif width_cols==2:
        new_width=double_col_width
    else:
        print('No valid width_cols value specified:', width_cols)
    height_over_width=width_height[1]/width_height[0]
    new_height=new_width*height_over_width
    fig= plt.figure(figsize=(new_width, new_height), constrained_layout=constrained_layout)
    
    return fig


def make_spec_hump(cent_wave,wave_sigma,amp):
    
    return

def clean_color_string(input_table, color_header='plot_color'):
    """
    Take in the hexadecimal color stored in an astropy table and fix it to not have extra 
    quotation marks in it that Excel loves to add because screw anyone trying to use 
    something not in Excel... apparently.
    
    """
    index_range=range(0,len(input_table))
    for index in index_range:
        color_string= input_table[color_header][index]
        stripped_string=color_string.replace('"', '')
        input_table[color_header][index]=stripped_string
    return input_table

def calculate_atmospheric_refraction(header, offset=0.):
    """
    use Wikipedia-provided simplified refraction calculation from like 1890 or whatever.
    There's an astropy method that I can dig into later if it seems worth it.
    
    This assumes a Goodman header and that all of the necessary fields are populated.
    
    offset in arcminutes up or down.
    """
    from astropy.units import cds
    cds.enable()
    pressure=(header['ENVPRE']*100.*u.pascal).to(cds.mmHg).value
    temperature=header['ENVTEM']
    cotangent=1./np.tan((header['MOUNT_EL']+(offset/60.))/180.*np.pi)
    
    
    return (21.5*pressure)/(273+temperature) * cotangent



def time_string():
    start = time.time()
    print(start)
    time_out=str(start).split('.')[0]
    return time_out






