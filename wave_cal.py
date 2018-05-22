"""
This script should extract the 1-d spectra and subtract the background region. It should also then wavelength calibrate the target spectra and output the wavelengths, subtracted spectrum, raw spectrum, and background spectrum as a fits file. This file should also do the barycentric velocity correction for the wavelength values. I guess it might as well also add a header for the BMJD_TDB time in the ouput fits files for the target.

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


speclistname = 'listCTB'
linefilename = 'JJ_FeAr_lines.txt'


parkes_location = coords.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
cerro_pachon_location = coords.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)


def to_barycenter(header):
    input_times = header['DATE-OBS'] #not gps-synched times
    obs_time = Time(input_times, format = 'isot', scale = 'utc',location = cerro_pachon_location)
    ra = header['RA']
    dec = header['DEC']
    target_coord = coords.SkyCoord(ra, dec, frame = 'icrs', unit= (u.hourangle, u.deg), )
    bary_corr =obs_time.tdb.light_travel_time(target_coord)
    bmjd_tdb_val = (obs_time.tdb+ bary_corr.tdb).mjd
    header.append(card = ('BMJD_TDB', bmjd_tdb_val, "value from DATE-OBS header"))
    return header

####
trace_band_mid= 105   #y-pixel that's about the center of the trace
trace_band_width = 20 #pixel width to determine the center of the trace
core_sides=  5
poly_degree = 3 #polynomial degree of the fit to the trace
bkg_width= core_sides
bkg_shift= 25
lamp_sigma_guess= 2
line_search_width = 3
lamp_p0 = [1000, 500,  lamp_sigma_guess, 0]

#####

fear_array= np.genfromtxt(linefilename, names = True)
line_x_checks = np.copy(fear_array['Pixel']) +90
print "line_x_checks should have just been created"
print line_x_checks
lamp_lines = np.copy(fear_array['User'])
line_sides = np.ones(line_x_checks.shape[0])*line_search_width


####

speclist = np.genfromtxt(speclistname, dtype = 'str')

#######


def gaussian_curve(x, a, x0, sigma,b):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+b

def fit_gaussian_curve(x_pixels, light_values, p0_list, search_width, plot_all = False):
    cut_region = np.where(x_pixels> (p0_list[1]-search_width ))
    #print '========'
    #print p0_list
    #print  "lower bound:", p0_list[1]-search_width
    #print "upper bound: ", p0_list[1]+search_width
    high_x_pixels= np.copy(x_pixels[cut_region])
    high_light_values= np.copy(light_values[cut_region])
    upper_cut = np.where(high_x_pixels < (p0_list[1]+search_width))
    cut_x_pixels = high_x_pixels[upper_cut]
    #print np.min(cut_x_pixels), np.max(cut_x_pixels), p0_list[1]
    cut_light_values= high_light_values[upper_cut]
    popt, pcov = sciop.curve_fit(gaussian_curve, cut_x_pixels, cut_light_values, p0= p0_list)
    #print "[amplitude, x0, sigma, b]"
    #print popt
    if plot_all:
        plt.plot(cut_x_pixels, cut_light_values, label = "data")
        plt.plot(cut_x_pixels, gaussian_curve(cut_x_pixels,popt[0],popt[1],popt[2],popt[3]),label ='fit')
        plt.legend()
        plt.show()
    else:
        pass
    #try:
        #cut_region = np.where(x_pixels> (popt[1]-search_width ))
        #print popt
        #print  "lower bound:", popt[1]-search_width
        #print "upper bound: ", popt[1]+search_width
        #high_x_pixels= np.copy(x_pixels[cut_region])
        #high_light_values= np.copy(light_values[cut_region])
        #upper_cut = np.where(high_x_pixels < (popt[1]+search_width))
        #cut_x_pixels = high_x_pixels[upper_cut]
        #print np.min(cut_x_pixels), np.max(cut_x_pixels), popt[1]
        #cut_light_values= high_light_values[upper_cut]
        #popt, pcov = sciop.curve_fit(gaussian_curve, cut_x_pixels, cut_light_values, p0= popt)
    #except ValueError as error:
        #print "probably chose a center outside the bounds of the image"
        #print error
        
    #print popt
    #print '========'
    #plt.plot(cut_x_pixels, cut_light_values)
    #plt.plot(cut_x_pixels, gaussian_curve(cut_x_pixels, popt[0], popt[1], popt[2], popt[3]))
    #plt.show()
    return popt, pcov




######
def get_trace_waves(target_med, lamp_im):
    target_band=target_med[trace_band_mid-trace_band_width/2:trace_band_mid+trace_band_width/2,:]
    band_inds= np.indices(target_band.shape)
    x_positions= band_inds[1,1]
    print 'xpositionsshape', x_positions.shape
    y_positions= np.argmax(target_band,axis=0)+(trace_band_mid-trace_band_width/2)
    print 'yshape', y_positions.shape
    polynomial_fit= np.polyfit(x_positions,y_positions,poly_degree)
    print polynomial_fit
    print polynomial_fit.shape
    poly_curve_y = np.polyval(polynomial_fit, x_positions)

    plt.imshow(np.log(img_data),cmap = 'hot', interpolation = 'none')
    plt.plot(x_positions, poly_curve_y, color = 'blue', label  = 'polynomial fit')
    plt.plot(x_positions,y_positions, color = 'black', label = 'max values', linestyle = 'none', marker = '*')
    plt.legend()
    plt.show()
        
    target_light= np.array([])
    bkg_light= np.array([])
    lamp_light= np.array([])
    print target_light.shape
    for x_pos in x_positions:
        xsum= np.sum(target_med[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides),x_pos])
        target_light= np.append(target_light,[xsum])
        bkg_sum= np.sum(target_med[np.int_(poly_curve_y[x_pos]+bkg_shift-core_sides):np.int_(poly_curve_y[x_pos]+bkg_shift+core_sides),x_pos])
        bkg_light= np.append(bkg_light,[bkg_sum])
        lamp_sum= np.sum(lamp_im[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides),x_pos])
        lamp_light= np.append(lamp_light,[lamp_sum])
    plt.plot(x_positions,target_light,'-')
    plt.xlabel('x (pixel)')
    plt.ylabel('Counts')
    plt.title('Target Spectrum')
    plt.show()
    
    target_light= target_light-bkg_light
    for x_spot in line_x_checks:
        plt.axvline( x= x_spot, color = 'r')
    #for x_spot in np.array(WaveList_Fe_930_12_24[0])/2.:
        #plt.axvline( x= x_spot, color = 'r')
    plt.plot(x_positions,lamp_light,'-')
    plt.xlabel('x (pixel)')
    plt.ylabel('Counts')
    plt.title('Lamp Spectrum (record corresponding dotted line and emission pixels)')
    #plt.yscale('log')
    plt.show()
    offset = 0

    #dotted_pixel = float(raw_input("dotted line pixel>>>"))
    #emission_pixel= float(raw_input("emission line pixel>>>"))
    #offset = emission_pixel-dotted_pixel
    print "skipping offsetting. Change lines 148 and 150 if you want otherwise."
   
    line_x_checks2 = np.copy(line_x_checks+offset)

    #for x_spot in line_x_checks2:
        #plt.axvline( x= x_spot, color = 'r')
    #for x_spot in np.array(WaveList_Fe_930_12_24[0])/2.:
        #plt.axvline( x= x_spot, color = 'r')
    #plt.plot(x_positions,lamp_light,'-')
    #plt.xlabel('x (pixel)')
    #plt.ylabel('Counts')
    #plt.title('Lamp Spectrum (offset applied)')
    ##plt.yscale('log')
    #plt.show()
    peaks_found=[]
    wave_peaks_found = []
    for lamp_line_guess,lamp_line_wave in zip( line_x_checks,lamp_lines):
        try:
            lamp_params, lamp_cov = fit_gaussian_curve(x_positions, lamp_light, [lamp_p0[0], lamp_line_guess, lamp_p0[2], lamp_p0[3]], line_search_width)
            if ((np.abs(lamp_params[0]) > 1.) and (np.abs(lamp_params[2])< 20) and (lamp_params[0] > 0) and (np.abs(lamp_line_guess-lamp_params[1]) < line_search_width) and  (np.abs(lamp_params[2])> 1)):
                peaks_found.append(lamp_params[1])
                wave_peaks_found.append(lamp_line_wave)
                #plt.plot(x_positions, lamp_light, label = 'lamp data', color = 'blue')
                #plt.plot(x_positions, gaussian_curve(x_positions, lamp_params[0], lamp_params[1], lamp_params[2], lamp_params[3]), color = 'r', label = 'Gaussian Fit')
                #plt.title("guess: " + str(lamp_line_guess) + ' fit:' + str(lamp_params[1]))
                #for x_spot in line_x_checks:
                    #plt.axvline( x= x_spot, color = 'k',linestyle = '--')
                #plt.axvline(x = lamp_line_guess, color = 'r', linestyle= '--')
                #plt.xlabel('Pixel')
                #plt.ylabel('Counts')
                #plt.legend()
                #plt.show()
            else:
                print "Gaussian too flat, flipped, or narrow (or not within the actual fitting region...):", lamp_params
        except RuntimeError as error:
            print error
    peaks_found = np.array(peaks_found)
    wave_peaks_found = np.array(wave_peaks_found)
    #print "line_x_checks:"
    #print line_x_checks
    #print "peaks found"
    #print peaks_found
    #print "wave_peaks_found"
    #print wave_peaks_found
    #for line,peak,wave in zip(line_x_checks, peaks_found, wave_peaks_found):
        #print line, peak, wave
    #polynomial fitting
    #poly_coeffs_lamp= np.polyfit(centroids,lamp_lines,2)
    poly_coeffs_lamp =np.polyfit(peaks_found, wave_peaks_found, 5)
    def x_to_wavelength(x_positions):
        #poly_curve_wavelength= poly_coeffs_lamp[2]+poly_coeffs_lamp[1]*x_positions + poly_coeffs_lamp[0]*(x_positions**2)
        poly_curve_wavelength= poly_coeffs_lamp[-1]+poly_coeffs_lamp[-2]*x_positions + poly_coeffs_lamp[-3]*(x_positions**2)+poly_coeffs_lamp[-4]*(x_positions**3)+poly_coeffs_lamp[-5]*(x_positions**4)+poly_coeffs_lamp[-6]*(x_positions**5)
        return poly_curve_wavelength
    poly_curve_wavelength= x_to_wavelength(x_positions)
    plt.plot(poly_curve_wavelength,target_light,'-')
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Counts')
    plt.title('Target Spectrum')
    #plt.yscale('log')
    #plt.ylim(10,200)
    plt.show()

    #plt.plot(x_positions, poly_curve_wavelength,  label = 'wavelength solution', color ='blue')
    #plt.plot(peaks_found, wave_peaks_found, marker= '*', linestyle = 'none', label = 'fitted values', color = 'red' )
    #plt.plot(line_x_checks2, lamp_lines, label = 'input points', color = 'green', marker = '*', linestyle = 'none')
    #plt.title("wavelength to pixel position")
    #plt.legend()
    #plt.show()

    plt.axhline(y=0 ,  label = 'wavelength solution', color ='blue')
    plt.plot(peaks_found, wave_peaks_found-x_to_wavelength(peaks_found), marker= '*', linestyle = 'none', label = 'fitted values', color = 'red' )
    plt.plot(line_x_checks2, lamp_lines-x_to_wavelength(line_x_checks2), label = 'input points', color = 'green', marker = '*', linestyle = 'none')
    plt.title("wavelength to pixel position Residuals")
    plt.xlabel('Pixel')
    plt.ylabel(r'Wavelength Residual $\AA$')
    plt.legend(loc= 'best')
    plt.show()
    return [polynomial_fit, poly_coeffs_lamp]
#######3


#####

last_file_lamp = False
target_stack = []
association_index = -1
polynomial_list = [] #should eventually be [[trace_polynomial, wavelength_fit_polynomial],[ trace...,wavelength...]]
#need to determine if the given image is a lamp or a target spectrum
for counter, img in enumerate(speclist):
    filename= glob(img)[0]
    if '_fe' in filename.lower():
        print 'Lamp file detected: ', filename
        print 'Updating lamp reference image.'
        lamp_i = fits.open(filename)
        lamp_header = fits.getheader(filename)
        lamp_im= lamp_i[0].data
        if last_file_lamp:
            #since the previous file was a lamp, that would make this a new run, so we'd want to use this lamp file, right?
            print "Double lamp detected, so it must be a new run."
            
            
        else:
            lamp_i = fits.open(filename)
            lamp_header = fits.getheader(filename)
            lamp_im= lamp_i[0].data
            #last_file_lamp = True
        last_file_lamp = True #since the image has to be a lamp
        association_index+=1
    else:
        #the filename doesn't contain a lamp indicator, so it must be a target spectrum
        print "Target file detected: ", filename
        i= fits.open(filename)
        header = fits.getheader(filename)
        img_data= i[0].data
        target_stack.append(img_data)
        last_file_lamp = False
        if '_fe' in speclist[counter+1].lower():
            print "Next file is a lamp, so we're going to do the trace and wavelength calibration."
            target_med = np.nanmedian(target_stack, axis = 0)
            new_coeffs= get_trace_waves(target_med, lamp_im)
            polynomial_list.append(new_coeffs)
            print "Resetting the target_stack."
            target_stack = [] #

        else:
            print "Next file is not a lamp."
        
print polynomial_list
        #last_file_lamp= False #since this image isn't a lamp


#def barycentric_vel_corr(header, wavelengths):
    #ra = header['RA']
    #dec = header['DEC']
    #radec = coords.SkyCoord(ra, dec, frame = 'icrs', unit= (u.hourangle, u.deg))
    #bary_corr = radec.radial_velocity_correction(obstime= Time(header['DATE-OBS'], format = 'isot', scale= 'utc'), location = cerro_pachon_location)
    #bary_corr = bary_corr.to(u.km/u.s)
    #lambda_rest = (wavelengths*(u.Angstrom))*const.c.to(u.km/u.s)/(-1*bary_corr+const.c.to(u.km/u.s))
    #lambda_rest = lambda_rest.value
    #return lambda_rest

last_file_lamp = False
target_stack = []
association_index = -1
new_filelist =[]
#need to determine if the given image is a lamp or a target spectrum
for counter, img in enumerate(speclist):
    filename= glob(img)[0]
    if '_fe' in filename.lower():
        print 'Lamp file detected: ', filename
        
        if last_file_lamp:
            #since the previous file was a lamp, that would make this a new run, so we'd want to use this lamp file, right?
            print "Double lamp detected, so it must be a new run."
            
            
        else:
            association_index+=1
            print "association_index increased: ", association_index
            #last_file_lamp = True
            
        last_file_lamp = True #since the image has to be a lamp
        
    else:
        #the filename doesn't contain a lamp indicator, so it must be a target spectrum
        print "Target file detected: ", filename
        
        i= fits.open(filename)
        header = fits.getheader(filename)
        img_data= i[0].data
        filename = 'w' + filename
        new_filelist.append(filename)
        polynomials = polynomial_list[association_index]
        band_inds= np.indices(img_data.shape)
        x_positions= band_inds[1,1]
        target_light= np.array([])
        bkg_light= np.array([])
        poly_curve_y = np.polyval(polynomials[0], x_positions)
        poly_curve_wavelength= np.polyval(polynomials[1], x_positions)
        for x_pos in x_positions:
            xsum= np.sum(img_data[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides),x_pos])
            target_light= np.append(target_light,[xsum])
            bkg_sum= np.sum(img_data[np.int_(poly_curve_y[x_pos]+bkg_shift-core_sides):np.int_(poly_curve_y[x_pos]+bkg_shift+core_sides),x_pos])
            bkg_light= np.append(bkg_light,[bkg_sum])
        plt.plot(x_positions,target_light,'-')
        plt.xlabel('x (pixel)')
        plt.ylabel('Counts')
        plt.title('Target Spectrum')
        plt.show()
        target_light= target_light-bkg_light
        header= to_barycenter(header) #append the BMJD_TDB value
        #poly_curve_wavelength= barycentric_vel_corr(header, poly_curve_wavelength) #correction of Earth's orbital motion
        hdu = fits.PrimaryHDU(poly_curve_wavelength, header = header)
        hdu1= fits.ImageHDU(target_light)
        hdu2= fits.ImageHDU(bkg_light)
        hdulist= fits.HDUList([hdu, hdu1, hdu2])
        hdulist.writeto(filename, overwrite= True)
        #target_stack.append(img_data)
        last_file_lamp = False
        if '_fe' in speclist[counter+1].lower():
            print "Next file is a lamp"
            pass
            #print "Next file is a lamp, so we're going to do the trace and wavelength calibration."
            

        else:
            print "Next file is not a lamp."
            
            
            
        
