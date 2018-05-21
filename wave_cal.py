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

####
trace_band_mid= 105   #y-pixel that's about the center of the trace
trace_band_width = 20 #pixel width to determine the center of the trace
core_sides=  5
poly_degree = 3 #polynomial degree of the fit to the trace
bkg_width= core_sides
bkg_shift= 50
lamp_sigma_guess= 2
line_search_width = 3
lamp_p0 = [100, 500,  lamp_sigma_guess, 0]

#####

fear_array= np.genfromtxt(linefilename, names = True)
line_x_checks = np.copy(fear_array['Pixel']) +90
lamp_lines = np.copy(fear_array['User'])
line_sides = np.ones(line_x_checks.shape[0])*line_search_width


####

speclist = np.genfromtxt(speclistname, dtype = 'str')

#######


def gaussian_curve(x, a, x0, sigma,b):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+b

def fit_gaussian_curve(x_pixels, light_values, p0_list, search_width, plot_all = False):
    cut_region = np.where(x_pixels> (p0_list[1]-search_width ))
    print '========'
    #print p0_list
    #print  "lower bound:", p0_list[1]-search_width
    #print "upper bound: ", p0_list[1]+search_width
    high_x_pixels= np.copy(x_pixels[cut_region])
    high_light_values= np.copy(light_values[cut_region])
    upper_cut = np.where(high_x_pixels < (p0_list[1]+search_width))
    cut_x_pixels = high_x_pixels[upper_cut]
    print np.min(cut_x_pixels), np.max(cut_x_pixels), p0_list[1]
    cut_light_values= high_light_values[upper_cut]
    popt, pcov = sciop.curve_fit(gaussian_curve, cut_x_pixels, cut_light_values, p0= p0_list)
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
        
    print popt
    print '========'
    #plt.plot(cut_x_pixels, cut_light_values)
    #plt.plot(cut_x_pixels, gaussian_curve(cut_x_pixels, popt[0], popt[1], popt[2], popt[3]))
    #plt.show()
    return popt, pcov




######
def get_trace_waves(target_stack, lamp_im):
    target_band=img_data[trace_band_mid-trace_band_width/2:trace_band_mid+trace_band_width/2,:]
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
        
    dotted_pixel = float(raw_input("dotted line pixel>>>"))
    emission_pixel= float(raw_input("emission line pixel>>>"))
    offset = emission_pixel-dotted_pixel
    line_x_checks = np.copy(line_x_checks+offset)

    for x_spot in line_x_checks:
        plt.axvline( x= x_spot, color = 'r')
    #for x_spot in np.array(WaveList_Fe_930_12_24[0])/2.:
        #plt.axvline( x= x_spot, color = 'r')
    plt.plot(x_positions,lamp_light,'-')
    plt.xlabel('x (pixel)')
    plt.ylabel('Counts')
    plt.title('Lamp Spectrum (offset applied)')
    #plt.yscale('log')
    plt.show()

#######3

last_file_lamp = False
target_stack = []
#need to determine if the given image is a lamp or a target spectrum
for counter, img in enumerate(speclist):
    filename= glob(img)[0]
    if filename.lower().contains('_fe'):
        print 'Lamp file detected: ', filename
        print 'Updating lamp reference image.'
        lamp_i = fits.open(filename)
        lamp_header = fits.getheader(filename)
        lamp_im= lamp_i[0].data
        #if last_file_lamp:
            ##since the previous file was a lamp, that would make this a new run, so we'd want to use this lamp file, right?
            #print "Double lamp detected, so it must be a new run."
            
        #else:
            #lamp_i = fits.open(filename)
            #lamp_header = fits.getheader(filename)
            #lamp_im= lamp_i[0].data
            #last_file_lamp = True
        #last_file_lamp = True #since the image has to be a lamp
    else:
        #the filename doesn't contain a lamp indicator, so it must be a target spectrum
        print "Target file detected: ", filename
        i= fits.open(target_file)
        header = fits.getheader(target_file)
        img_data= i[0].data
        target_stack.append(img_data)
        
        #last_file_lamp= False #since this image isn't a lamp
        
            
            
            
        
