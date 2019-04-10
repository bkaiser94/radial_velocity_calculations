"""
This script should extract the 1-d spectra and subtract the background region. It should also then wavelength
calibrate the target spectra and output the wavelengths, subtracted spectrum, raw spectrum, and background
spectrum as a fits file.  I guess it might as well also add a header for the BMJD_TDB time in the ouput fits files
for the target.

This file should create the error spectrum. I should probably tack it on after the sky spectrum. I cannot remove
the sky spectrum because it provides a ready check on if a star or something fell into the background region.

STEP 3 of Reduction

The ListCTB file is already  created by straight_reduction.py
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

zerolistname= 'listZero'

speclistname = 'listCTB'
flatlistname='listCTBflat'
masterflatfile= 'mctb.master_flat.fits'
linefilename = 'JJ_FeAr_lines.txt'
zerolist = np.genfromtxt(zerolistname, dtype ='str')
print "zerolist.shape",zerolist.shape
n_biases= zerolist.shape[0]
parkes_location = coords.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
cerro_pachon_location = coords.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)

skip_flat= True

def to_barycenter(header):
    #input_times = header['DATE-OBS'] #not gps-synched times
    input_year = header['OPENDATE'] #gps-synched date
    input_hours = header['OPENTIME'] #gps-synched time
    exp_time= header['EXPTIME']*u.s
    input_times = input_year+'T'+input_hours #formatting correctly
    obs_time = Time(input_times, format = 'isot', scale = 'utc',location = cerro_pachon_location)
    obs_time= obs_time+exp_time/2.
    ra = header['RA']
    dec = header['DEC']
    target_coord = coords.SkyCoord(ra, dec, frame = 'icrs', unit= (u.hourangle, u.deg), )
    bary_corr =obs_time.tdb.light_travel_time(target_coord)
    bmjd_tdb_val = (obs_time.tdb+ bary_corr.tdb).mjd
    header.append(card = ('BMJD_TDB', bmjd_tdb_val, "exp. midpoint value from OPENDATE and OPENTIME headers"))
    return header

####
#trace_band_mid= 85   #y-pixel that's about the center of the trace #old one as of 2018-10-31
#trace_band_mid= 95   #y-pixel that's about the center of the trace J1431
trace_band_mid=100 #y-pixel for Keaton's object 2019-03-07 2019-03-25 commented out
trace_band_width = 40 #pixel width to determine the center of the trace 2019-03-25 commented out
#trace_band_mid=95 #y-pixel for secondary of wisea0615 2019-03-07
#trace_band_mid=115 #y-pixel for actual wisea0615
#trace_band_width = 10 #pixel width to determine the center of the trace
core_sides=  5
#core_sides=  7
bkg_core_sides= core_sides #This should be changed most likely to make the value be higher to further reduce noise.
y_trace_width= core_sides*2+1 #the actual number of pixels in the vertical direction that are in the trace (or background)
poly_degree = 3 #polynomial degree of the fit to the trace
flat_poly= 7
#bkg_shift= 25 #2019-03-25 commented out
bkg_shift = 50
lamp_sigma_guess= 2
line_search_width = 3
lamp_p0 = [1000, 500,  lamp_sigma_guess, 0]
lamp_bounds = ([0,-np.inf,0,0],[30000,np.inf,20,5000 ])
seeing_range = [1200, 1220]
#seeing_p0= [1000, trace_band_width/2, lamp_sigma_guess, 0] #p0 list for the gaussian fit to the vertical
#seeing_p0= [1000, 5, lamp_sigma_guess, 0] #p0 list for the gaussian fit to the vertical
#seeing_p0= [1000, 20, lamp_sigma_guess, 0] #p0 list for the gaussian fit to the vertical
#seeing_p0= [2000, 20, lamp_sigma_guess, 0] #p0 list for the gaussian fit to the vertical
#see_fit_bounds = ([50, 0, 0.7, 0],[18000, 1000, trace_band_width, 2000]) #(lower, upper) bounds on the fit for the seeing.
seeing_p0= [2000, 20, lamp_sigma_guess, 0] #p0 list for the gaussian fit to the vertical
see_fit_bounds = ([50, 0, 0.7, 0],[1e8, 40, trace_band_width, 1e8]) #(lower, upper) bounds on the fit for the seeing.
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

def fit_gaussian_curve(x_pixels, light_values, p0_list, search_width, plot_all = False, bounds = (-np.inf, np.inf)):
    """
    Those bounds are the default for scipy.optimize.curve_fit(), so now changing them changes the bounds
    """
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
    popt, pcov = sciop.curve_fit(gaussian_curve, cut_x_pixels, cut_light_values, p0= p0_list, bounds = bounds)
    #print "[amplitude, x0, sigma, b]"
    #print popt
    if plot_all:
        print "popt", popt
        print "bounds", bounds
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


def normalize_flat(masterflatfile=masterflatfile, plot_all = False):
    """
    Normalize the master flat file by fitting a polynomial to the trace-extraction region to remove the spectral features of the lamp
    
    Assumes it's a Quartz lamp, but I would think this should work for any lamp...
    
    I'M NOT DEALING WITH FLAT UNCERTAINTIES SINCE IT'S GOING TO BE A DIVISION, WHICH MEANS ASYMMETRIC ERRORS AND THEY ARE SMALL ENOUGH TO START WITH THAT I'M NOT TRACKING THEM FOR NOW
    
    """
    i= fits.open(masterflatfile)
    header = fits.getheader(masterflatfile)
    master_flat= i[0].data
    master_flat_err= i[1].data #normalized sigma values to the original counts values. 
    readnoise = header['RDNOISE']
    roi= master_flat[trace_band_mid-trace_band_width/2:trace_band_mid+trace_band_width/2,:] #region for searching for the trace in the future. we're using it to determine the spectrum polynomial of the flat lamp
    band_inds= np.indices(roi.shape)
    x_positions= band_inds[1,1]
    summed_roi= np.mean(roi, axis=0)
    polynomial_fit = np.polyfit(x_positions,summed_roi, flat_poly)
    poly_curve= np.polyval(polynomial_fit, x_positions)
    if plot_all:
        plt.plot(x_positions, summed_roi, color ='b', label= 'data')
        plt.plot(x_positions, poly_curve, color = 'r', label= 'polynomial fit')
        plt.title('Flat polynomial fit')
        plt.legend()
        plt.show()
        plt.plot(x_positions, summed_roi/poly_curve)
        plt.title('Divided by polynomial fit')
        plt.show()
        plt.plot(x_positions, roi[0]/poly_curve)
        plt.title('single row Divided by polynomial fit')
        plt.show()
        plt.plot(x_positions, summed_roi-poly_curve)
        plt.title('Residuals of  polynomial fit')
        plt.show()
    else:
        pass
    normed_flat = master_flat/poly_curve
    print "max value in normalized flat" , np.max(normed_flat)
    print "min value in normalized flat", np.min(normed_flat)
    return normed_flat

######
def get_trace_waves(target_med, lamp_im):
    target_band=target_med[trace_band_mid-trace_band_width/2:trace_band_mid+trace_band_width/2,:]
    band_inds= np.indices(target_band.shape)
    x_positions= band_inds[1,1]
    y_pos = band_inds[0].T[0]
    print 'xpositionsshape', x_positions.shape
    y_positions= np.argmax(target_band,axis=0)+(trace_band_mid-trace_band_width/2)
    print 'yshape', y_positions.shape
    
    #seeing_band = np.sum(np.copy(target_band[:,seeing_range[0]:seeing_range[1]]),axis=1)
    #seeing_popt, seeing_pcov = fit_gaussian_curve(y_pos, seeing_band, seeing_p0, trace_band_width, plot_all=True, bounds = see_fit_bounds)
    #seeing_sigma = seeing_popt[2]
    #print seeing_popt
    #print "Seeing sigma: ", seeing_popt[2]
    polynomial_fit= np.polyfit(x_positions,y_positions,poly_degree)
    print polynomial_fit
    print polynomial_fit.shape
    poly_curve_y = np.polyval(polynomial_fit, x_positions)
    std_dev = np.std(poly_curve_y-y_positions)
    plt.imshow(np.log(img_data),cmap = 'hot', interpolation = 'none')
    plt.plot(x_positions, poly_curve_y, color = 'blue', label  = 'polynomial fit')
    plt.plot(x_positions,y_positions, color = 'black', label = 'max values', linestyle = 'none', marker = '*')
    plt.plot(x_positions, poly_curve_y+core_sides, color = 'blue', linestyle = '--')
    plt.plot(x_positions, poly_curve_y-core_sides, color = 'blue', linestyle= '--')
    plt.plot(x_positions, np.int_(poly_curve_y+bkg_shift), color = 'cyan', label = 'background')
    plt.plot(x_positions, np.int_(poly_curve_y+bkg_shift-core_sides), color = 'cyan', linestyle= '--')
    plt.plot(x_positions, np.int_(poly_curve_y+bkg_shift+core_sides), color = 'cyan', linestyle = '--')
    plt.legend()
    plt.show()
    plt.imshow(target_band[:,seeing_range[0]:seeing_range[1]], cmap='hot')
    plt.show()
    seeing_band = np.sum(np.copy(target_band[:,seeing_range[0]:seeing_range[1]]),axis=1)
    seeing_p0[1]=np.argmax(seeing_band)
    seeing_popt, seeing_pcov = fit_gaussian_curve(y_pos, seeing_band, seeing_p0, trace_band_width, plot_all=True, bounds = see_fit_bounds)
    seeing_sigma = seeing_popt[2]
    print seeing_popt
    print "Seeing sigma: ", seeing_popt[2]
    target_light= np.array([])
    bkg_light= np.array([])
    lamp_light= np.array([])
    print target_light.shape
    for x_pos in x_positions:
        xsum= np.sum(target_med[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides+1),x_pos])
        target_light= np.append(target_light,[xsum])
        bkg_sum= np.sum(target_med[np.int_(poly_curve_y[x_pos]+bkg_shift-core_sides):np.int_(poly_curve_y[x_pos]+bkg_shift+core_sides+1),x_pos])
        bkg_light= np.append(bkg_light,[bkg_sum])
        lamp_sum= np.sum(lamp_im[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides+1),x_pos])
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
    #offset = 0

    dotted_pixel = float(raw_input("dotted line pixel>>>"))
    emission_pixel= float(raw_input("emission line pixel>>>"))
    #dotted_pixel=0
    #emission_pixel=0
    offset = emission_pixel-dotted_pixel
    #print "skipping offsetting. Change lines 261 - 264 if you want otherwise."
   
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
    for lamp_line_guess,lamp_line_wave in zip( line_x_checks2,lamp_lines):
        try:
            lamp_params, lamp_cov = fit_gaussian_curve(x_positions, lamp_light, [lamp_p0[0], lamp_line_guess, lamp_p0[2], lamp_p0[3]], line_search_width, bounds= lamp_bounds)
            if ((np.abs(lamp_params[0]) > 1.) and (np.abs(lamp_params[2])< 20) and (lamp_params[0] > 0) and (np.abs(lamp_line_guess-lamp_params[1]) < line_search_width) and  (np.abs(lamp_params[2])> 1)):
                peaks_found.append(lamp_params[1])
                wave_peaks_found.append(lamp_line_wave)
                plt.plot(x_positions, lamp_light, label = 'lamp data', color = 'blue')
                plt.plot(x_positions, gaussian_curve(x_positions, lamp_params[0], lamp_params[1], lamp_params[2], lamp_params[3]), color = 'r', label = 'Gaussian Fit')
                plt.title("guess: " + str(lamp_line_guess) + ' fit:' + str(lamp_params[1]))
                for x_spot in line_x_checks2:
                    plt.axvline( x= x_spot, color = 'k',linestyle = '--')
                plt.axvline(x = lamp_line_guess, color = 'r', linestyle= '--')
                plt.xlabel('Pixel')
                plt.ylabel('Counts')
                plt.legend()
                plt.show()
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
    return [polynomial_fit, poly_coeffs_lamp], seeing_sigma
#######3


######3 Flat handling

print "skip_flat=", skip_flat
normed_flat= normalize_flat(plot_all = True)
print "\n==================\n"
print "skip_flat=", skip_flat
print "\n\n"
if skip_flat:
    normed_flat= np.ones(normed_flat.shape)
    print "skipping flat-fielding.\nYes, I know it just went to all the trouble of calculating the flat stuff."
print "\n\n==================\n"


#####

last_file_lamp = False
target_stack = []
association_index = -1
polynomial_list = [] #should eventually be [[trace_polynomial, wavelength_fit_polynomial],[ trace...,wavelength...]]
sigma_list= [] #to be appended to the headers for 
seeing_list = []
#need to determine if the given image is a lamp or a target spectrum
for counter, img in enumerate(speclist):
    filename= glob(img)[0]
    #if '_fe' in filename.lower():
    if '_fe.' in filename.lower():
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
        img_data= np.copy(i[0].data)
        img_data= img_data/normed_flat #division by the flat. The noise from the flat is not accounted for currently
        target_stack.append(img_data)
        last_file_lamp = False
        #if '_fe' in speclist[counter+1].lower():
        if '_fe.' in speclist[counter+1].lower():
            print "Next file is a lamp, so we're going to do the trace and wavelength calibration."
            print "Using last lamp file as calibration lamp"
            filename= speclist[counter+1]
            lamp_i = fits.open(filename)
            lamp_header = fits.getheader(filename)
            lamp_im= lamp_i[0].data
            target_med = np.nanmedian(target_stack, axis = 0)
            new_coeffs, seeing_sig= get_trace_waves(target_med, lamp_im)
            sigma_list.append(seeing_sig)
            seeing_list.append(2*np.sqrt(2*np.log(2))*seeing_sig) #assuming normal distribution for that
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
    #if '_fe' in filename.lower():
    if '_fe.' in filename.lower():
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
        img_data= np.copy(i[0].data)
        img_data= img_data/normed_flat #also have to divide it here... since the other place was for seeing
        filename = 'w' + filename
        new_filelist.append(filename)
        polynomials = polynomial_list[association_index]
        seeing_sig = sigma_list[association_index]
        seeing_FWHM = seeing_list[association_index]
        band_inds= np.indices(img_data.shape)
        x_positions= band_inds[1,1]
        target_light= np.array([])
        target_noise2_list= np.array([])
        bkg_light= np.array([])
        bkg_noise2_list= np.array([])
        poly_curve_y = np.polyval(polynomials[0], x_positions)
        poly_curve_wavelength= np.polyval(polynomials[1], x_positions)
        for x_pos in x_positions:
            trace_vals=img_data[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides+1),x_pos]
            xsum= np.sum(trace_vals)
            #xsum= np.sum(img_data[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides+1),x_pos]) #old way 2018-10-31
            target_light= np.append(target_light,[xsum])
            target_noise2 = np.copy(xsum+trace_vals.shape[0]*header['RDNOISE']**2+trace_vals.shape[0]*header['RDNOISE']**2/n_biases)
            target_noise2_list= np.append(target_noise2_list, [target_noise2])
            up_bkg=img_data[np.int_(poly_curve_y[x_pos]+bkg_shift-bkg_core_sides):np.int_(poly_curve_y[x_pos]+bkg_shift+bkg_core_sides+1),x_pos]
            down_bkg= img_data[np.int_(poly_curve_y[x_pos]-bkg_shift-core_sides):np.int_(poly_curve_y[x_pos]-bkg_shift+bkg_core_sides+1),x_pos]
            bkg_comb= np.append(up_bkg, down_bkg)
            #print "trace_vals.shape", trace_vals.shape
            #print "bkg_comb.shape", bkg_comb.shape
            bkg_noise2= trace_vals.shape[0]*np.copy(np.mean(bkg_comb)/bkg_comb.shape[0]+header['RDNOISE']**2/bkg_comb.shape[0]+header['RDNOISE']**2/(bkg_comb.shape[0]*n_biases))
            #bkg_sum= np.sum(img_data[np.int_(poly_curve_y[x_pos]+bkg_shift-core_sides):np.int_(poly_curve_y[x_pos]+bkg_shift+core_sides+1),x_pos])
            bkg_sum= trace_vals.shape[0]*np.copy(np.mean(bkg_comb)) #take the mean of the bkg portion of the sky
            bkg_light= np.append(bkg_light,[bkg_sum])
            bkg_noise2_list= np.append(bkg_noise2_list, [bkg_noise2]) #list of noise values for a single pixel (resulting from the mean of the sky) for a given column
            
        #plt.plot(x_positions,target_light,'-')
        #plt.xlabel('x (pixel)')
        #plt.ylabel('Counts')
        #plt.title('Target Spectrum')
        #plt.show()
        #noise_spectrum = np.copy(np.sqrt(target_light + bkg_light + y_trace_width*header['RDNOISE'])) #old way 2018-10-31
        noise_spectrum= np.copy(np.sqrt(target_noise2_list+bkg_noise2_list)) #combination of noises of the background pixels and the target trace pixels.
        print "noise_spectrum.shape", noise_spectrum.shape
        target_light= target_light-bkg_light
        noise_spectrum = noise_spectrum/target_light #normalized noise values by the target spectrum, so now unitless.
        plt.plot(x_positions,target_light,'-')
        plt.xlabel('x (pixel)')
        plt.ylabel('Counts')
        plt.title('Target Spectrum')
        plt.show()
        header= to_barycenter(header) #append the BMJD_TDB value
        header.append(card= ("pix_scal", 0.3, ' "/pixel'))
        header.append(card = ('see_sig', seeing_sig, 'Sigma of Gauss seeing fit (pixels)'))
        header.append(card = ('see_FWHM', seeing_FWHM, 'Seeing (pixels)'))
        header.append(card = ('skipflat', skip_flat, 'flatfielding skipped or not'))
        #poly_curve_wavelength= barycentric_vel_corr(header, poly_curve_wavelength) #correction of Earth's orbital motion
        hdu = fits.PrimaryHDU(poly_curve_wavelength, header = header)
        hdu1= fits.ImageHDU(target_light)
        hdu2= fits.ImageHDU(bkg_light)
        hdu3 = fits.ImageHDU(noise_spectrum)
        hdulist= fits.HDUList([hdu, hdu1, hdu2, hdu3])
        hdulist.writeto(filename, overwrite= True)
        #target_stack.append(img_data)
        last_file_lamp = False
        if '_fe' in speclist[counter+1].lower():
            print "Next file is a lamp"
            pass
            #print "Next file is a lamp, so we're going to do the trace and wavelength calibration."
            

        else:
            print "Next file is not a lamp."
            
            
            
        
