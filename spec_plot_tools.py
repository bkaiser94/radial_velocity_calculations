"""
Created by Ben Kaiser (UNC-Chapel Hill) (date not known of original creation.)


"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy import constants as const
from astropy.time import Time
from astropy.table import Table
import scipy.interpolate as scinterp
from astropy import coordinates as coords

import balmer_line_ranges as blr
import cal_params as cp

percentile= 50

soar_area= np.pi*(cp.soar_diameter/2.)**2 #area of SOAR light-gathering in meters (approximately)
parkes_location = coords.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
cerro_pachon_location = coords.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)



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
    low_mask = np.where(wave_array < lower_bound)
    high_mask= np.where(wave_array > upper_bound)
    low_waves= wave_array[low_mask]
    high_waves= wave_array[high_mask]
    low_other = other_array[low_mask]
    high_other= other_array[high_mask]
    #merge_waves= np.append(low_waves, high_waves)
    #merge_other= np.append(low_other, high_other)
    return


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
    return clean_spec


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


def retrieve_sdss_spec(filename,scale_noise=True):
    """
    Input: filename for a spectrum from SDSS (so it has the SDSS headers and fits format)
    
    Output:    Output: Spectrum made of a 2xN numpy array, header of the fits file you loaded it from, and noise spectrum.
    
    I.E. the same output format as retrieve_spec(), so that this is effectively a way of bringing all of the spectra together in the same format that I'm already tooled to use.
    
    """
    spec_hdu= fits.open(filename)
    spec_array=spec_hdu[1].data
    waves= 10.**np.copy(spec_array['loglam'])
    flux= np.copy(spec_array['flux'])
    flux= flux/10. #convert from 10**-17 to 10**-16
    try:
        noise= np.copy(spec_array['PropErr'])
    except KeyError as error:
        print(error)
        print('setting noise=1')
        noise=np.ones(waves.shape[0])
    file_spec=np.vstack([waves, flux])
    if scale_noise:
        file_noise_spec= np.vstack([waves, noise])
    else:
        file_noise_spec=np.vstack([waves, noise/flux])
    header= spec_hdu[1].header
    return file_spec, header, file_noise_spec

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
        print(count, row['intens'])
        try:
            row['intens']=int(row['intens'])
            nist_table[count]['intens']=float(row['intens'])
            print('converted to int!', row['intens'])
        except ValueError:
            print('ValueError ^^^^')
    #ne_table['intens']=ne_table['intens'].astype(float) #changing this column that gets read as strings for whatever reason
    nist_table.pprint()
    return nist_table

def show_plot(show_telluric=True, show_legend=True):
    if show_legend:
        plt.legend(loc='best')
    else:
        pass
    if show_telluric:
        plot_telluric()
    else:
        pass
    plt.show()
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

    
    
