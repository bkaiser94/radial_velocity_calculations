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
import os
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
#import cosmics
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const
from astropy.modeling import models as asmodels
from astropy.modeling import fitting as asfitting
from astropy.table import Table, Column

import get_cal_params as gcp
import cal_params as cp
import spec_plot_tools as spt
import spectral_extraction_tools as spext



zerolistname= 'listZero'

speclistname = 'listCTB'
masterflatfile= 'mctb.master_flat.fits'
#linefilename = 'JJ_FeAr_lines.txt'
zerolist = np.genfromtxt(zerolistname, dtype ='str')
print("zerolist.shape",zerolist.shape)
n_biases= zerolist.shape[0]
parkes_location = coords.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
cerro_pachon_location = coords.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)

skip_flat= True
need_offset=True
do_save_wavesoln=True
do_save_tracetrack=True
#trace_method='maxes'
trace_method='binned_gauss'


def to_barycenter(header):
    #input_times = header['DATE-OBS'] #not gps-synched times
    try:
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
    except KeyError as error:
        print("\n\n")
        print("==========")
        print("Problem converting to BMJD_TDB time, so there won't be these times for this exposure and subsequent data from it.")
        print(error)
        header.append(card = ('BMJD_TDB','Bad', "exp. midpoint value from OPENDATE and OPENTIME headers"))
        print('===========')
        print('\n\n')
    return header

####
#trace_offset =-16.3#amount by which the calculated trace needs to be offset to end up on the dimmer desired target. Should normally be 0 unless doing a specific extraction.
trace_offset=-6
#trace_band_mid= 110   #y-pixel that's about the center of the trace #old one as of 2018-10-31
#trace_band_mid= 96   #y-pixel that's about the center of the trace J1431
#trace_band_mid=100
trace_band_mid=105 #usual extraction search center point
#trace_band_mid= 112 #y-pixel for SDSSJ1159 400M1
#trace_band_mid= 90 #y-pixel for SDSSJ1159 400M2
#trace_band_mid=129
#trace_band_mid=77 #
#trace_band_width=10
#trace_band_width = 100 #pixel width to determine the center of the trace 2019-03-25 commented out
#trace_band_width = 60#pixel width to determine the center of the trace 2019-03-25 commented out
trace_band_width=190#usual extraction search window
#trace_band_width= 14 #SDSSJ1159
#trace_band_mid=95 #y-pixel for secondary of wisea0615 2019-03-07
#trace_band_mid=115 #y-pixel for actual wisea0615
#trace_band_width = 15 #pixel width to determine the center of the trace
#sigma_multi_side= 4 #multiple of sigma value of trace gaussian that should be distance out to go for extraction window
sigma_multi_side=1.5 #multiple of sigma value of trace gaussian that should be distance out to go for extraction window
#sigma_multi_side=2
#sigma_multi_side=3
#sigma_multi_side=2.5
#sigma_multi_side= 1 #multiple of sigma value of trace gaussian that should be distance out to go for extraction window

#core_sides=  5
#core_sides=  7

#y_trace_width= core_sides*2+1 #the actual number of pixels in the vertical direction that are in the trace (or background)
poly_degree = 3 #polynomial degree of the fit to the trace
lamp_poly_degree=5 #setting prior to change on 2025-01-09. This was the default for everything up to this date.
#lamp_poly_degree=3 #new setting changed to on 2025-01-09 for HgAr 400mshort spectra because the red is doubling back.
flat_poly= 7
#bkg_shift= 25 #2019-03-25 commented out
#bkg_shift = 50 #20190412 previously in place
#bkg_shift= 30 #standard shift used
#bkg_shift=20
bkg_shift=40 #shift used frequently for spectrophotometric standard extractions
#bkg_shift=35
#bkg_shift= 47
#bkg_shift=10
#bkg_shift=48
#bkg_shift=6
#bkg_core_sides= 2*core_sides #This should be changed most likely to make the value be higher to further reduce noise.
#bkg_side_multi= 1. #
#bkg_side_multi= 1.5 #mutliple of core_sides that that  bkg_core_sides should be later
#bkg_side_multi=2. #mutliple of core_sides that that  bkg_core_sides should be later
bkg_side_multi=3. #mutliple of core_sides that that  bkg_core_sides should be later
#bkg_side_multi=4. #mutliple of core_sides that that  bkg_core_sides should be later

bkg_max_side= bkg_shift/2.-5
lamp_sigma_guess= 2
line_search_width = 3#formerly 3 20190502
#lamp_p0 = [1000, 500,  lamp_sigma_guess, 0]
lamp_p0 = [10000, 500,  lamp_sigma_guess, 0]
#lamp_bounds = ([0,-np.inf,0,0],[30000,np.inf,20,5000 ])
lamp_bounds = ([1,-np.inf,0,0],[1e6,np.inf,40,10000 ])
seeing_range = [1200, 1220]
#seeing_p0= [1000, trace_band_width/2, lamp_sigma_guess, 0] #p0 list for the gaussian fit to the vertical
#seeing_p0= [1000, 5, lamp_sigma_guess, 0] #p0 list for the gaussian fit to the vertical
#seeing_p0= [1000, 20, lamp_sigma_guess, 0] #p0 list for the gaussian fit to the vertical
#seeing_p0= [2000, 20, lamp_sigma_guess, 0] #p0 list for the gaussian fit to the vertical
#see_fit_bounds = ([50, 0, 0.7, 0],[18000, 1000, trace_band_width, 2000]) #(lower, upper) bounds on the fit for the seeing.
seeing_p0= [2000, 20, lamp_sigma_guess, 0] #p0 list for the gaussian fit to the vertical
#see_fit_bounds = ([50, 0, 0.7, 0],[1e8, trace_band_width, trace_band_width, 1e8]) #(lower, upper) bounds on the fit for the seeing.
#see_fit_bounds = ([5, 0, 0.7, 0],[1e8, trace_band_width, 7.1, 1e8]) #(lower, upper) bounds on the fit for the seeing. #way before 2021-01-25
see_fit_bounds = ([5, 0, 0.5, 0],[1e8, trace_band_width, 7.1, 1e8]) #(lower, upper) bounds on the fit for the seeing. #way after 2021-01-25
#minimum sigma value corresponds to seeing of 0.5" for 2x2 binned pixels and max is 5" seeing

box_dict= {
    'amplitude':10,
    'x_0':1000,
    'width':5,
    'slope':1
    }
######

expedited_wavecals=False
do_airglow_corr=True

default_skyline_search_width=40 #True default. This was what was used prior to 2022-05-10. I have returned to it because it served me well previously.
#default_skyline_search_width=60 #Going big. While at home. 2022-05-10


#air_off_type='lambda' #if you want the airglow offset to be applied in wavelength space, i.e. subtract a lambda value from all wavelength values
air_off_type='pixel' #if you want the offset to be applied in pixel space
#air_off_type='none' #setting for not applying the airglow correction. Realistically you should just set do_airglow_corr=False for this option


bkg_method= 'avg' #background subtraction method; 'avg'
#bkg_method= 'poly' #background subtraction method; 'avg' means the regions will be averaged together to be subtracted from each pixel in the trace region.
#bkg_method= 'opp_get_bright' #background subtraction method; gets the bright part of the frame, meaning trace_offset==0 should be true, and then the bkg_shift corresponds to ...
#bkg_method= 'opp_get_dim' #background subtraction method; gets the dim part of the frame, meaning trace_offset should be set to however far off the extraction should be , and then the bkg_shift corresponds to ...

bkg_poly= 2 #polynomial degree for background fitting.
#bkg_poly= 1 #polynomial degree for background fitting.

contaminant_offset= 0#offset to get to middle of contaminant. used only with bkg_method='opp_get_bright'. Should be equal to value used for the trace_offset when doing the 'opp_get_dim' extraction for the same frame.

if bkg_method == 'opp_get_dim':
    if trace_offset==0:
        print('trace_offset==',trace_offset)
        print('So ', bkg_method, "can't be used. Exiting script")
        sys.exit()
    else:
        pass
        #bkg_shift=int(-1.6*trace_offset)
        bkg_shift=-2*trace_offset
elif bkg_method=='opp_get_bright':
    if trace_offset!=0:
        print('trace_offset==',trace_offset)
        print('So ', bkg_method, "can't be used. Exiting script")
        sys.exit()
    else:
        bkg_shift= 2*contaminant_offset
else:
    pass





#fear_array= np.genfromtxt(linefilename, names = True)
#line_x_checks = np.copy(fear_array['Pixel']) +90
#print "line_x_checks should have just been created"
#print line_x_checks
#lamp_lines = np.copy(fear_array['User'])
#line_sides = np.ones(line_x_checks.shape[0])*line_search_width


#####

speclist = np.genfromtxt(speclistname, dtype = 'str')

arc_im_filename=speclist[1] #need to change this to be an explicit reference to whatever place actually indicates the lamp image.
header=fits.getheader(arc_im_filename)


#####
setup_dict= gcp.get_cal_params(header)
#fear_array= np.genfromtxt(setup_dict['linelistname'], names = True)
fear_array= np.genfromtxt(cp.line_list_dir+ setup_dict['linelistname'], names = True)
try:
    fear_array= fear_array[np.where(fear_array['use']>0)]
except KeyError:
    print('missing "use" column')
    pass
line_x_checks = np.copy(fear_array['Pixel']) +setup_dict['offset']
print("line_x_checks should have just been created")
print( line_x_checks)
lamp_lines = np.copy(fear_array['User'])
line_sides = np.ones(line_x_checks.shape[0])*line_search_width


####


#######

def make_wavesoln_filename(fits_filename):
    """
    make outupt filename for the wavelength solution for a given file
    
    """
    name_parts = fits_filename.split('.')
    output_filename= 'wsoln_'+name_parts[1]+'.txt'
    return output_filename

def save_wavesoln(fits_filename, wave_polynomial):
    wave_filename= make_wavesoln_filename(fits_filename)
    wave_filename= cp.wave_sol_dir+ wave_filename
    if not os.path.exists(cp.wave_sol_dir):
        os.makedirs(cp.wave_sol_dir)
    else:
        pass
    print('Saving',wave_filename,'.')
    np.savetxt(wave_filename, wave_polynomial)
    print(wave_filename,' saved.')
    return wave_filename

def save_tracetrack(fits_filename,trace_polynomial):
    name_parts = fits_filename.split('.')
    output_filename= 'tracetrack_'+name_parts[1]+'.txt'
    output_filename=cp.tracetrack_dir+output_filename
    if not os.path.exists(cp.tracetrack_dir):
        os.makedirs(cp.tracetrack_dir)
    else:
        pass
    print('Saving',output_filename,'.')
    np.savetxt(output_filename, trace_polynomial)
    print(output_filename,' saved.')
    return output_filename



def gaussian_curve(x, a, x0, sigma,b):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+b

def seeing_window(seeing_sigma):
    ext_sides= int(seeing_sigma*sigma_multi_side)+1 #setting the extraction window based on the seeing and rounding up, by adding 1
    #bkg_core_sides= np.min([int(ext_sides*bkg_side_multi), int(bkg_max_side), int(bkg_shift-ext_sides)])
    bkg_core_sides= np.min([int(ext_sides*bkg_side_multi),int(bkg_shift-ext_sides-1)])
    return ext_sides, bkg_core_sides

def fit_gaussian_curve(x_pixels, light_values, p0_list, search_width, plot_all = False, bounds = (-np.inf, np.inf), fixed_width=True):
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
    #print('p0_list it lets you use:', p0_list)
    popt, pcov = sciop.curve_fit(gaussian_curve, cut_x_pixels, cut_light_values, p0= p0_list, bounds = bounds)
    #print "[amplitude, x0, sigma, b]"
    #print popt
    if plot_all:
        print("popt", popt)
        print("bounds", bounds)
        plt.plot(cut_x_pixels, cut_light_values, label = "data")
        plt.plot(cut_x_pixels, gaussian_curve(cut_x_pixels,popt[0],popt[1],popt[2],popt[3]),label ='fit')
        #popt[1]=popt[1]+trace_offset
        plt.axvline(x=popt[1],  color='b')
        if fixed_width:
            pass
        else:
            core_sides, bkg_core_sides= seeing_window(popt[2])
        #print('\n\n\n=======\n')
        #print(filename)
        #print('bkg_core_sides right before cross section plot', bkg_core_sides)
        #print('\n=========\n\n\n')
        plt.axvline(x=popt[1]+core_sides, linestyle= '--', color='b')
        plt.axvline(x=popt[1]-core_sides,linestyle='--',  color='b')
        plt.axvline(x=popt[1]+bkg_shift,  color='cyan')
        plt.axvline(x=popt[1]+bkg_shift+bkg_core_sides, linestyle= '--', color='cyan')
        plt.axvline(x=popt[1]+bkg_shift-bkg_core_sides,linestyle='--',  color='cyan')
        plt.axvline(x=popt[1]-bkg_shift,  color='cyan')
        plt.axvline(x=popt[1]-bkg_shift+bkg_core_sides, linestyle= '--', color='cyan')
        plt.axvline(x=popt[1]-bkg_shift-bkg_core_sides,linestyle='--',  color='cyan')
        plt.title('Fit to seeing in x pixel range '+str(p0_list[1]-search_width)+'-'+str(p0_list[1]+search_width)+'\n'+'Constrained to '+str(trace_band_width)+' pixels centered on ' + str(trace_band_mid) +'\n'+'So it excludes ' + str(trace_band_mid-(trace_band_width/2)) + ' pixels on the low side\n' + 'and ' +str(200-trace_band_width/2-trace_band_mid) + ' pixels on the upper side')
        plt.legend()
        plt.show()
    else:
        pass
    return popt, pcov

def get_skyline_bounds(header, x_pixels):
    try:
        slit_width= header['SLIT'].split('"')[0] #separate the slit header to get the arcsecond value of the slit width
        slit_width = float(slit_width) #it was a string
    except ValueError as error:
        print('ValueError:', error)
        print('The Goodman headers were changed for no apparent reason on 2019-07-30, so we have to try another format for the slit width')
        slit_width= header['SLIT'].split('_')[0] 
        slit_width = float(slit_width) #it was a string
    #slit_width = float(slit_width) #it was a string
    xpix_scale, ypix_scale= spt.get_pixel_scale(header)
    slit_width_pix= slit_width/xpix_scale
    skyline_bounds={
        'amplitude':cp.slit_airline_ampbounds,
        'x_0': (np.nanmin(x_pixels), np.nanmax(x_pixels)),
        'width': (slit_width_pix*0.5, slit_width_pix*1.5)
        }
    #not doing the zeropoint bounds because I'm not sure how I'm going to do it yet...
    
    return skyline_bounds

def make_box_model(p0_dict, bounds=(-np.inf, np.inf), fixed={'amplitude':False, 'x_0:': False, "width": False}):
    box_model=asmodels.Box1D(amplitude=p0_dict['amplitude'], x_0=p0_dict['x_0'], width=p0_dict['width'], bounds= bounds, fixed=fixed)
    return box_model

def make_trap_model(p0_dict, bounds=(-np.inf, np.inf)):
    p0_dict.update({'slope': 2})
    bounds.update({'slope': (0.1, np.inf)})
    trap_model=asmodels.Trapezoid1D(amplitude=p0_dict['amplitude'], x_0=p0_dict['x_0'], width=p0_dict['width'], slope=p0_dict['slope'], bounds= bounds)
    return trap_model


def fit_slitskyline_function(x_pixels, light_values,  header, p0_dict= cp.slit_airline_p0, plot_all = False):
    skyline_bounds= get_skyline_bounds(header, x_pixels)
    #slit_model= make_box_model(p0_dict, bounds=skyline_bounds)
    #cont_box_p0= cp.cont_box_p0
    #cont_box_p0['x_0']= p0_dict['x_0']
    #cont_model= make_box_model(cont_box_p0, fixed={'amplitude': False, 'width': True, 'x_0':True}, bounds= cp.cont_box_bounds)
    #sky_model = slit_model
    #slit_model= make_trap_model(p0_dict, bounds=skyline_bounds)
    #cont_model = asmodels.Polynomial1D(degree=1)
    #box_model.bounds= {'x_0': (0, x_pixels.shape[0])} #the middle of the box plot can't be outside the range of the thing to examine; not sure if this goes away with next definition on next line
    #box_model.bounds= get_skyline_bounds(header, x_pixels )
    #fitter = asfitting.LevMarLSQFitter()
    #fitter= asfitting.SLSQPLSQFitter()
    #fitted_model = fitter(slit_model, x_pixels, light_values)
    #fitted_model = fitter(sky_model, x_pixels, light_values)
    
    p0_list_sky= [p0_dict['amplitude'], p0_dict['x_0'], p0_dict['width']/2., 1.]
    #print('p0_list:' , lamp_p0)
    print('p0_list_sky:',p0_list_sky)
    #popt, pcov = sciop.curve_fit(gaussian_curve, x_pixels, light_values,  p0=p0_list_sky) #prior to 2022-05-10
    
    temp_bounds = ([1., -1*np.inf, -1*np.inf,-1*np.inf],[np.inf, np.inf, np.inf, np.inf])#created 2022-05-10 to restrict the skyline gaussian function to being positive (concave down).
    popt, pcov = sciop.curve_fit(gaussian_curve, x_pixels, light_values,  p0=p0_list_sky, bounds=temp_bounds)

    #print('slit_model x_0:', slit_model.x_0)
    #print('fitted_model x_0:', fitted_model.x_0)
    #print('fitted_model width:' ,fitted_model.width)
    #print('fitted_model amplitude:', fitted_model.amplitude)
    #print('fitted_model fixed:', fitted_model.fixed)
    #print('slit_model fixed:', slit_model.fixed)
    #try:
        #print('fitted_model slope:' ,fitted_model.slope)
    #except AttributeError as error:
        #print('AttributeError:', error)
    #print('slit_model:',slit_model)
    #print('fitted_model:', fitted_model)
    #print('initial guesses:', p0_list_sky)
    #print('Gaussian vals:', popt)
    
    if plot_all:
        plt.plot(x_pixels, light_values, color='b')
        #plt.plot(x_pixels, fitted_model(x_pixels), color= 'r', label='model')
        plt.plot(x_pixels, gaussian_curve(x_pixels,popt[0],popt[1],popt[2],popt[3]),label ='gaussian fit')
        plt.axvline(x=popt[1],  color='g', label='gaussian center', linestyle='--')
        plt.axvline(x=p0_list_sky[1], color='k', label='guessed value', linestyle='--')
        #plt.axvline(x= fitted_model.x_0, linestyle='--', color='k', label = 'center')
        plt.xlabel('pixels')
        plt.ylabel('intensity')
        plt.legend(loc='best')
        plt.title('Sky background slit function fitting')
        plt.show()
    return popt[1]

def wavelength_to_pixel(lambda_val, in_wave_coeffs):
    """
    input wave_coeffs should already have an offset subtracted from the x-values everywhere.... you can't really 
    do that...
    """
    wave_coeffs= np.copy(in_wave_coeffs)
    wave_coeffs[-1]= wave_coeffs[-1]-lambda_val
    
    def func_to_solve(x):
        if lamp_poly_degree==5:
            return wave_coeffs[0]*x**5+ wave_coeffs[1]*x**4 +wave_coeffs[2]*x**3+ wave_coeffs[3]*x**2+wave_coeffs[4]*x + wave_coeffs[5]
        else:
            print("don't have function to solve for inversion of wavelengths for that lamp_poly_degree:", lamp_poly_degree)
            return np.polyval(wave_coeffs, x)
    #plt.plot(np.polyval(wave_coeffs, np.linspace(0,2000,2000)), label='changed wave_coeffs')
    #plt.plot(np.polyval(in_wave_coeffs, np.linspace(0,2000,2000)), label='og wave_coeffs')
    #plt.plot(func_to_solve(np.linspace(0,2000,2000)),label='func_to_solve')
    #plt.legend(loc='best')
    #plt.show()
    pixel= sciop.brentq(func_to_solve, 0,2100)
    return pixel


def find_skyline_offset(x_pixels, light_values, airline_lambda, wave_coeffs, header, search_width=default_skyline_search_width, initial_offset=0, plot_all=False):
    
    airline_guess= wavelength_to_pixel(airline_lambda, wave_coeffs)
    if plot_all:
        plt.axvline(x=airline_guess, color= 'r')
        plt.plot(x_pixels, light_values, color='b')
        plt.title('Auto-generated prediction of pixel position for air line')
        plt.show()
    else:
        pass

    cut_region = np.where(x_pixels> (airline_guess-search_width ))
    high_x_pixels= np.copy(x_pixels[cut_region])
    high_light_values= np.copy(light_values[cut_region])
    upper_cut = np.where(high_x_pixels < (airline_guess+search_width))
    cut_x_pixels = high_x_pixels[upper_cut]
    #print np.min(cut_x_pixels), np.max(cut_x_pixels), p0_list[1]
    cut_light_values= high_light_values[upper_cut]
    p0_dict= cp.slit_airline_p0
    max_ind= np.argmax(cut_light_values)
    max_point= cut_x_pixels[max_ind]
    if plot_all:
        plt.axvline(x=airline_guess, color= 'magenta', label='wave to pix guess')
    else:
        pass
    p0_dict['x_0']= max_point
    try:
        corr_pixel= fit_slitskyline_function(cut_x_pixels, cut_light_values, header, p0_dict= p0_dict, plot_all=plot_all)
        corr_wave= np.polyval(wave_coeffs, corr_pixel) #''wavelength'' of skyline on original scale
        offset_lambda= airline_lambda-corr_wave
        offset= airline_guess-corr_pixel #offset that should be added to the pixel values for the new determinations
    except RuntimeError as error:
        print(error)
        print('automatically setting this one skyline offset as 0')
        offset_lambda=0
        offset=0
    return offset, offset_lambda



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
    print("max value in normalized flat" , np.max(normed_flat))
    print("min value in normalized flat", np.min(normed_flat))
    return normed_flat

######
def get_trace_waves(target_med, lamp_im, do_wavelengths=True, poly_coeffs_lamp=[0], trace_method=trace_method,filename=''):
    try:
        target_band=target_med[trace_band_mid-trace_band_width/2:trace_band_mid+trace_band_width/2,:]
    except TypeError:
        target_band=target_med[int(trace_band_mid-trace_band_width/2):int(trace_band_mid+trace_band_width/2),:]
    band_inds= np.indices(target_band.shape)
    x_positions= band_inds[1,1]
    y_pos = band_inds[0].T[0]
    
    #20190624 moved this section up here
    plt.imshow(target_band[:,seeing_range[0]:seeing_range[1]], cmap='hot')
    plt.show()
    seeing_band = np.sum(np.copy(target_band[:,seeing_range[0]:seeing_range[1]]),axis=1)
    seeing_p0[1]=np.argmax(seeing_band)
    try:
        seeing_popt, seeing_pcov = fit_gaussian_curve(y_pos, seeing_band, seeing_p0, trace_band_width, plot_all=True, bounds = see_fit_bounds, fixed_width=False)
    except RuntimeError as error:
        print("RuntimeError:", error)
        seeing_popt=see_fit_bounds[0]
        seeing_popt[2]=1.
        print('\n\n*************\nWARNING! SEEING COULD NOT BE FIT FOR THIS OBJECT IN THE NORMAL BAND! ASSIGNING VALUE OF', seeing_popt,' FOR ALL OF seeing_popt AS THAT SHOULD BE THE BEST POSSIBLE. EXAMINE RAW IMAGES TO MAKE SURE THERE IS NOT AN ISSUE!\n**************\n\n')
        print('That corresponds to a seeing of', 2*np.sqrt(2*np.log(2))*seeing_popt[2]*0.3, '"\nand an extraction width of ', seeing_popt[2]*sigma_multi_side,'pixels\n\n')
    seeing_sigma = seeing_popt[2]
    print(seeing_popt)
    print("Seeing sigma: ", seeing_popt[2])
    core_sides, bkg_core_sides= seeing_window(seeing_sigma)
    ### end of 20190624 moved section
    
    
    if trace_method=='maxes':
        print('xpositionsshape', x_positions.shape)
        y_positions= np.argmax(target_band,axis=0)+(trace_band_mid-trace_band_width/2)
        print('yshape', y_positions.shape)
    elif trace_method=='binned_gauss':
        rebinned_im, rebinned_indices = spt.rebin_image(target_band, rebin_axis=1, rebin_num=10)
        rebinned_imT= rebinned_im.T
        rebinned_x_positions= rebinned_indices.T
        y_positions=[]
        max_fluxes=[]
        coll_seeing_sigmas= []
        rebin_counter= 0
        for pixel_column in rebinned_imT:
            #seeing_p0[1]=np.argmax(pixel_column)
            try:
                if rebin_counter%25==26:
                    subset_popt, subset_pcov = fit_gaussian_curve(y_pos, pixel_column, seeing_p0, trace_band_width, plot_all=True, bounds = see_fit_bounds, fixed_width=False)
                else:
                    subset_popt, subset_pcov = fit_gaussian_curve(y_pos, pixel_column, seeing_p0, trace_band_width, plot_all=False, bounds = see_fit_bounds, fixed_width=False)
            except RuntimeError as error:
                print('\n\n++++++++++')
                print(error)
                print("Gaussian seeing fit didn't converge for column")
                print("Setting amplitude=0, x0=trace_band_width/2., sigma=0., b=0.")
                print('++++++++++\n\n')
                subset_popt=[0., trace_band_width/2., 0., 0.]
            if (subset_popt[2]>(see_fit_bounds[0][2]+0.05)):
                #print(subset_popt[2], '>', see_fit_bounds[0][2])
                y_positions.append(subset_popt[1]+(trace_band_mid-trace_band_width/2))
                max_fluxes.append(subset_popt[0])
                coll_seeing_sigmas.append(subset_popt[2])
            else:
                print("seeing sigma was on the low boundary so we're not including it: x=", rebinned_x_positions[0][rebin_counter], subset_popt[2], 'pixels', 'seeing would be', subset_popt[2]*0.3*2*np.sqrt(2*np.log(2)),'"')
                y_positions.append(np.nan)
                max_fluxes.append(np.nan)
                coll_seeing_sigmas.append(np.nan)
            rebin_counter+=1
        y_positions=np.array(y_positions)
        x_positions=rebinned_x_positions[0]
        #print('x_positions', x_positions)
        max_fluxes= np.array(max_fluxes)
        coll_seeing_sigmas=np.array(coll_seeing_sigmas)
        
        nanmask= np.isnan(y_positions)
        y_positions= y_positions[~nanmask]
        x_positions= x_positions[~nanmask]
        max_fluxes=max_fluxes[~nanmask]
        coll_seeing_sigmas=coll_seeing_sigmas[~nanmask]
        
        plt.scatter(x_positions,2*np.sqrt(2*np.log(2))*coll_seeing_sigmas)
        plt.ylabel('Seeing (FWHM) in pixels')
        plt.xlabel('x pixel')
        plt.show()
        
        plt.scatter(max_fluxes,2*np.sqrt(2*np.log(2))*coll_seeing_sigmas)
        plt.ylabel('Seeing (FWHM) in pixels')
        plt.xlabel('amplitude of gaussian fit')
        plt.show()
        
        plt.scatter(x_positions, max_fluxes)
        plt.ylabel('amplitude of gaussian fit')
        plt.xlabel('x position')
        plt.title('binned pseudo-spectrum')
        plt.show()
        
            
    else:
        print('No trace_method specified that works (possibly none at all)')
        print('trace_method:', trace_method)
    
    #seeing_band = np.sum(np.copy(target_band[:,seeing_range[0]:seeing_range[1]]),axis=1)
    #seeing_popt, seeing_pcov = fit_gaussian_curve(y_pos, seeing_band, seeing_p0, trace_band_width, plot_all=True, bounds = see_fit_bounds)
    #seeing_sigma = seeing_popt[2]
    #print seeing_popt
    #print "Seeing sigma: ", seeing_popt[2]
    polynomial_fit= np.polyfit(x_positions,y_positions,poly_degree)
    print(polynomial_fit)
    print(polynomial_fit.shape)
    #20190605 added this offsetting part
    polynomial_fit[-1]=polynomial_fit[-1]+trace_offset
    plotting_x_coords= np.indices(img_data.shape)[1,1]
    poly_curve_y = np.polyval(polynomial_fit, plotting_x_coords)
    def bkg_trace(x_positions, sign='minus'):
        #return  np.int_(poly_curve_y[x_positions]+bkg_shift)
        if sign=='minus':
            return  spt.discrete_int(poly_curve_y[x_positions])-bkg_shift
        elif sign=='plus':
            return spt.discrete_int(poly_curve_y[x_positions])+bkg_shift
    std_dev = np.std(np.polyval(polynomial_fit, x_positions)-y_positions)
    #plt.imshow(np.log10(img_data),cmap = 'hot', interpolation = 'none')
    plt.imshow(np.sqrt(img_data),cmap = 'hot', interpolation = 'none')
    #plt.plot(plotting_x_coords, poly_curve_y, color = 'blue', label  = 'polynomial fit')
    plt.plot(x_positions,y_positions, color = 'black', label = 'max values', linestyle = 'none', marker = '*')
    #plt.plot(plotting_x_coords, poly_curve_y+core_sides, color = 'blue', linestyle = '--')
    #plt.plot(plotting_x_coords, poly_curve_y-core_sides, color = 'blue', linestyle= '--')
    plt.plot(plotting_x_coords, spt.discrete_int(poly_curve_y), color = 'blue', label  ='discrete polynomial fit')
    plt.plot(plotting_x_coords, spt.discrete_int(poly_curve_y)+core_sides, color = 'blue', linestyle = '--')
    plt.plot(plotting_x_coords, spt.discrete_int(poly_curve_y)-core_sides, color = 'blue', linestyle= '--')
    #plt.plot(x_positions, np.int_(poly_curve_y+bkg_shift), color = 'cyan', label = 'background')
    #plt.plot(x_positions, np.int_(poly_curve_y+bkg_shift-core_sides), color = 'cyan', linestyle= '--')
    #plt.plot(x_positions, np.int_(poly_curve_y+bkg_shift+core_sides), color = 'cyan', linestyle = '--')
    plt.plot(plotting_x_coords, bkg_trace(plotting_x_coords), color = 'cyan', label = 'background')
    plt.plot(plotting_x_coords,bkg_trace(plotting_x_coords)-bkg_core_sides, color = 'cyan', linestyle= '--')
    plt.plot(plotting_x_coords, bkg_trace(plotting_x_coords)+bkg_core_sides, color = 'cyan', linestyle = '--')
    plt.plot(plotting_x_coords, bkg_trace(plotting_x_coords, sign='plus'), color = 'cyan', label = 'background')
    plt.plot(plotting_x_coords,bkg_trace(plotting_x_coords, sign='plus')-bkg_core_sides, color = 'cyan', linestyle= '--')
    plt.plot(plotting_x_coords, bkg_trace(plotting_x_coords, sign='plus')+bkg_core_sides, color = 'cyan', linestyle = '--')
    
    #print('\n\n\n=======\n')
    #print(filename)
    #print('bkg_core_sides right after image plot', bkg_core_sides)
    #print('\n=========\n\n\n')
    plt.legend()
    plt.show()
    
    plt.plot(plotting_x_coords, poly_curve_y, color = 'blue', label  = 'polynomial fit')
    plt.plot(x_positions,y_positions, color = 'black', label = 'max values', linestyle = 'none', marker = '*')
    plt.plot(plotting_x_coords, poly_curve_y+core_sides, color = 'blue', linestyle = '--')
    plt.plot(plotting_x_coords, poly_curve_y-core_sides, color = 'blue', linestyle= '--')
    plt.plot(plotting_x_coords, spt.discrete_int(poly_curve_y), color = 'green', label  ='discrete polynomial fit')
    plt.plot(plotting_x_coords, spt.discrete_int(poly_curve_y)+core_sides, color = 'green', linestyle = '--')
    plt.plot(plotting_x_coords, spt.discrete_int(poly_curve_y)-core_sides, color = 'green', linestyle= '--')
    #plt.plot(x_positions, np.int_(poly_curve_y+bkg_shift), color = 'cyan', label = 'background')
    #plt.plot(x_positions, np.int_(poly_curve_y+bkg_shift-core_sides), color = 'cyan', linestyle= '--')
    #plt.plot(x_positions, np.int_(poly_curve_y+bkg_shift+core_sides), color = 'cyan', linestyle = '--')

    plt.legend()
    plt.ylabel('y pixel')
    plt.xlabel('x pixel')
    plt.show()
    #plt.imshow(target_band[:,seeing_range[0]:seeing_range[1]], cmap='hot')
    #plt.show()
    #seeing_band = np.sum(np.copy(target_band[:,seeing_range[0]:seeing_range[1]]),axis=1)
    #seeing_p0[1]=np.argmax(seeing_band)
    #seeing_popt, seeing_pcov = fit_gaussian_curve(y_pos, seeing_band, seeing_p0, trace_band_width, plot_all=True, bounds = see_fit_bounds)
    #seeing_sigma = seeing_popt[2]
    #print seeing_popt
    #print "Seeing sigma: ", seeing_popt[2]
    target_light= np.array([])
    bkg_light= np.array([])
    lamp_light= np.array([])
    print(target_light.shape)
    for x_pos in plotting_x_coords:
        xsum= np.sum(target_med[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides+1),x_pos])
        target_light= np.append(target_light,[xsum])
        #bkg_sum= np.sum(target_med[np.int_(poly_curve_y[x_pos]+bkg_shift-core_sides):np.int_(poly_curve_y[x_pos]+bkg_shift+core_sides+1),x_pos])
        bkg_sum= np.sum(target_med[bkg_trace(x_pos)-core_sides:bkg_trace(x_pos)+core_sides+1, x_pos])
        #if bkg_method=='avg':
        #elif bkg_method=='poly':
            #pass
        #else:
            #print('No valid bkg_method selected')
            #print('bkg_method:', bkg_method)
            #pass
        bkg_light= np.append(bkg_light,[bkg_sum])
        #lamp_sum= np.sum(lamp_im[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides+1),x_pos])
        lamp_sum= np.average(lamp_im[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides+1),x_pos])
        lamp_light= np.append(lamp_light,[lamp_sum])
    if do_wavelengths:
        plt.plot(plotting_x_coords,target_light,'-')
        plt.xlabel('x (pixel)')
        plt.ylabel('Counts')
        plt.title('Target Spectrum')
        plt.show()
    else:
        pass
    target_light= target_light-bkg_light
    if do_wavelengths:
        global need_offset
        if need_offset:
            for x_spot in line_x_checks:
                plt.axvline( x= x_spot, color = 'r',linestyle='--')
            #for x_spot in np.array(WaveList_Fe_930_12_24[0])/2.:
                #plt.axvline( x= x_spot, color = 'r')
            plt.plot(plotting_x_coords,lamp_light,'-')
            plt.xlabel('x (pixel)')
            plt.ylabel('Counts')
            plt.title('Lamp Spectrum (record corresponding dotted line and emission pixels)')
            #plt.yscale('log')
            plt.show()
            #offset = 0
            try:
                dotted_pixel = float(raw_input("dotted line pixel>>>"))
                emission_pixel= float(raw_input("emission line pixel>>>"))
            except NameError as error:
                print("NameError:", error)
                print('Presumably in Python 3 instead of 2 because "raw_input" is not recognized')
                dotted_pixel = float(input("dotted line pixel>>>"))
                emission_pixel= float(input("emission line pixel>>>"))
            #dotted_pixel=0
            #emission_pixel=0
            offset = emission_pixel-dotted_pixel
            global store_offset
            store_offset= offset
            #global need_offset
            need_offset=False
            #print "skipping offsetting. Change lines 261 - 264 if you want otherwise."
        else:
            print("Using stored offset because we've already done it")
            offset=store_offset
    
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
        for x_spot in line_x_checks2:
            plt.axvline( x= x_spot, color = 'k',linestyle = '--')
        plt.plot(plotting_x_coords, lamp_light, label = 'lamp data', color = 'blue')

        for lamp_line_guess,lamp_line_wave in zip( line_x_checks2,lamp_lines):
            try:
                lamp_params, lamp_cov = fit_gaussian_curve(plotting_x_coords, lamp_light, [lamp_p0[0], lamp_line_guess, lamp_p0[2], lamp_p0[3]], line_search_width, bounds= lamp_bounds)
                #if ((np.abs(lamp_params[0]) > 1.) and (np.abs(lamp_params[2])< 20) and (lamp_params[0] > 0) and (np.abs(lamp_line_guess-lamp_params[1]) < line_search_width) and  (np.abs(lamp_params[2])> 1)):
                if ((np.abs(lamp_params[0]) > 1.)  and (lamp_params[0] > 0) and (np.abs(lamp_line_guess-lamp_params[1]) < line_search_width) and  (np.abs(lamp_params[2])> 1e-3)):
                    print('Going one layer deeper on the fit')
                    lamp_params, lamp_cov = fit_gaussian_curve(plotting_x_coords, lamp_light, lamp_params, line_search_width, bounds= lamp_bounds)
                    if ((np.abs(lamp_params[0]) > 1.)  and (lamp_params[0] > 0) and (np.abs(lamp_line_guess-lamp_params[1]) < line_search_width) and  (np.abs(lamp_params[2])> 1e-3)):
                        print('survived the second cut')
                        peaks_found.append(lamp_params[1])
                        wave_peaks_found.append(lamp_line_wave)
                        #plt.plot(plotting_x_coords, lamp_light, label = 'lamp data', color = 'blue')
                        plt.plot(plotting_x_coords, gaussian_curve(plotting_x_coords, lamp_params[0], lamp_params[1], lamp_params[2], lamp_params[3]), color = 'r', label = 'Gaussian Fit')
                        plt.title("guess: " + str(lamp_line_guess) + ' fit:' + str(lamp_params[1]))
                        #for x_spot in line_x_checks2:
                            #plt.axvline( x= x_spot, color = 'k',linestyle = '--')
                        plt.axvline(x = lamp_line_guess, color = 'r', linestyle= '--')
                    else:
                        print('did not survive second cut')
                        
                        #commented 2021-06-11 and replaced with the second level of hell above.
                    #peaks_found.append(lamp_params[1])
                    #wave_peaks_found.append(lamp_line_wave)
                    ##plt.plot(plotting_x_coords, lamp_light, label = 'lamp data', color = 'blue')
                    #plt.plot(plotting_x_coords, gaussian_curve(plotting_x_coords, lamp_params[0], lamp_params[1], lamp_params[2], lamp_params[3]), color = 'r', label = 'Gaussian Fit')
                    #plt.title("guess: " + str(lamp_line_guess) + ' fit:' + str(lamp_params[1]))
                    ##for x_spot in line_x_checks2:
                        ##plt.axvline( x= x_spot, color = 'k',linestyle = '--')
                    #plt.axvline(x = lamp_line_guess, color = 'r', linestyle= '--')
                    
                    
                    
                    #plt.xlabel('Pixel')
                    #plt.ylabel('Counts')
                    #plt.legend()
                    #plt.show()
                else:
                    print( "Gaussian too flat, flipped, or narrow (or not within the actual fitting region...):", lamp_params)
                    #plt.plot(plotting_x_coords, lamp_light, label = 'lamp data', color = 'blue')
                    #plt.plot(plotting_x_coords, gaussian_curve(plotting_x_coords, lamp_params[0], lamp_params[1], lamp_params[2], lamp_params[3]), color = 'r', label = 'Gaussian Fit')
                    #plt.title("guess: " + str(lamp_line_guess) + ' fit:' + str(lamp_params[1]))
                    #for x_spot in line_x_checks2:
                        #plt.axvline( x= x_spot, color = 'k',linestyle = '--')
                    #plt.axvline(x = lamp_line_guess, color = 'r', linestyle= '--')
                    #plt.xlabel('Pixel')
                    #plt.ylabel('Counts')
                    #plt.legend()
                    #plt.show()
            except RuntimeError as error:
                print(error)
        def get_core_name(filename=filename):
            #leading_characters=5 #number of characters before the actual "name" of the object, this will be the 0123_ portion of the filename. I'm going to split out the stuff on either side of the decimals
            interior_name=filename.split('.')[1]
            #core_name=interior_name[leading_characters:] #nevermind on removing the numerical characters because if I break up a given target into more than one run it will only retain the extraction from a single run. So while it will be messier to have files output in such a way that it appears only one specific frame (e.g. 0123_objectname_setup) was measured for residuals etc, it will allow retention of more information if you know how to look.
            core_name=interior_name
            return core_name
        core_name=get_core_name()
        
        plt.xlabel('Pixel')
        plt.ylabel('Counts')
        plt.title('all fitting results')
        plt.show()
        peaks_found = np.array(peaks_found)
        wave_peaks_found = np.array(wave_peaks_found)
        if not os.path.exists(cp.wave_sol_dir):
            os.makedirs(cp.wave_sol_dir)
        else:
            pass
        np.savetxt(cp.wave_sol_dir+core_name+'_measured_pixel_coords.txt', np.append([peaks_found],[wave_peaks_found],axis=0).T, delimiter='\t')
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
        poly_coeffs_lamp =np.polyfit(peaks_found, wave_peaks_found, lamp_poly_degree)
        def x_to_wavelength(plotting_x_coords):
            poly_curve_wavelength=np.polyval(poly_coeffs_lamp, plotting_x_coords)
            return poly_curve_wavelength
        poly_curve_wavelength= x_to_wavelength(plotting_x_coords)
        plt.plot(poly_curve_wavelength,target_light,'-')
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.ylabel('Counts')
        plt.title('Target Spectrum')
        #plt.yscale('log')
        #plt.ylim(10,200)
        plt.show()

        plt.plot(plotting_x_coords, poly_curve_wavelength,  label = 'wavelength solution', color ='blue')
        plt.plot(peaks_found, wave_peaks_found, marker= '*', linestyle = 'none', label = 'fitted values', color = 'red' )
        plt.plot(line_x_checks2, lamp_lines, label = 'input points', color = 'green', marker = '*', linestyle = 'none')
        plt.title("wavelength to pixel position")
        plt.legend()
        plt.show()

        plt.axhline(y=0 ,  label = 'wavelength solution', color ='blue')
        plt.plot(peaks_found, wave_peaks_found-x_to_wavelength(peaks_found), marker= '*', linestyle = 'none', label = 'fitted values', color = 'red' )
        plt.plot(line_x_checks2, lamp_lines-x_to_wavelength(line_x_checks2), label = 'input points', color = 'green', marker = '*', linestyle = 'none')
        plt.title("wavelength to pixel position Residuals")
        plt.xlabel('Pixel')
        plt.ylabel(r'Wavelength Residual $\AA$')
        plt.legend(loc= 'best')
        plt.show()
        
        plt.axhline(y=0 ,  label = 'wavelength solution', color ='blue')
        plt.plot(wave_peaks_found, wave_peaks_found-x_to_wavelength(peaks_found), marker= '*', linestyle = 'none', label = 'fitted values', color = 'red' )
        #plt.plot(line_x_checks2, lamp_lines-x_to_wavelength(line_x_checks2), label = 'input points', color = 'green', marker = '*', linestyle = 'none')
        plt.title("wavelength position Residuals")
        plt.xlabel('Wavelength')
        plt.ylabel(r'Wavelength Residual $\AA$')
        plt.legend(loc= 'best')
        plt.show()
        
        #So upon reflection, I think it makes more sense to record the residuals of wavelengths in the solution
        residuals_array=np.append([wave_peaks_found],[wave_peaks_found-x_to_wavelength(peaks_found)], axis=0).T
        np.savetxt(cp.wave_sol_dir+core_name+'_wavelength_residuals.txt', residuals_array, delimiter='\t')
    else:
        print("skipping wavelength fitting\npresumably because you already did it, so there's no need to do it again.")
        pass
    return [polynomial_fit, poly_coeffs_lamp], seeing_sigma
#######3


######3 Flat handling

print ("skip_flat=", skip_flat)
if skip_flat:
    example_hdu= fits.open(speclist[1])
    example_im= example_hdu[0].data
    normed_flat= np.ones(example_im.shape)
    print ("skipping flat-fielding.\nYes, I know it just went to all the trouble of calculating the flat stuff.")
else:
    #normed_flat= normalize_flat(plot_all = True)
    i= fits.open(masterflatfile)
    header = fits.getheader(masterflatfile)
    master_flat= i[0].data
    master_flat_err= i[1].data #normalized sigma values to the original counts values. 
    readnoise = header['RDNOISE']
    normed_flat=master_flat
print ("\n==================\n")
print ("skip_flat=", skip_flat)
print ("\n\n")
print ("\n\n==================\n")


#####

last_file_lamp = False
target_stack = []
association_index = -1
polynomial_list = [] #should eventually be [[trace_polynomial, wavelength_fit_polynomial],[ trace...,wavelength...]]
sigma_list= [] #to be appended to the headers for 
seeing_list = []
#need to determine if the given image is a lamp or a target spectrum
do_wavelengths=True
poly_coeffs_lamp=[0]
for counter, img in enumerate(speclist):
    filename= glob(img)[0]
    #if '_fe' in filename.lower():
    if '_fe.' in filename.lower():
        print ('Lamp file detected: ', filename)
        print ('Updating lamp reference image.')
        lamp_i = fits.open(filename)
        lamp_header = fits.getheader(filename)
        lamp_im= lamp_i[0].data
        if last_file_lamp:
            #since the previous file was a lamp, that would make this a new run, so we'd want to use this lamp file, right?
            print ("Double lamp detected, so it must be a new run.")
            
            
        else:
            lamp_i = fits.open(filename)
            lamp_header = fits.getheader(filename)
            lamp_im= lamp_i[0].data
            #last_file_lamp = True
        last_file_lamp = True #since the image has to be a lamp
        association_index+=1
    else:
        #the filename doesn't contain a lamp indicator, so it must be a target spectrum
        print( "Target file detected: ", filename)
        i= fits.open(filename)
        header = fits.getheader(filename)
        img_data= np.copy(i[0].data)
        if not skip_flat:
            img_data= img_data/normed_flat #division by the flat. The noise from the flat is not accounted for currently
            print('doing a division by a flat, but it should be 1s if skip_flat=True.')
        else:
            print('actually skipping the flat and not dividing by anything.')
        target_stack.append(img_data)
        last_file_lamp = False
        #if '_fe' in speclist[counter+1].lower():
        if '_fe.' in speclist[counter+1].lower():
            print("Next file is a lamp, so we're going to do the trace and wavelength calibration.")
            print("Using last lamp file as calibration lamp")
            prev_filename=filename
            filename= speclist[counter+1]
            lamp_i = fits.open(filename)
            lamp_header = fits.getheader(filename)
            lamp_im= lamp_i[0].data
            target_med = np.nanmedian(target_stack, axis = 0)
            #new_coeffs, seeing_sig= get_trace_waves(target_med, lamp_im)
            if expedited_wavecals:
                new_coeffs, seeing_sig= get_trace_waves(target_med, lamp_im, do_wavelengths=do_wavelengths, poly_coeffs_lamp=poly_coeffs_lamp,filename=prev_filename)
                poly_coeffs_lamp=new_coeffs[1]
                do_wavelengths=False
            else:
                new_coeffs, seeing_sig= get_trace_waves(target_med, lamp_im,filename=prev_filename)
            sigma_list.append(seeing_sig)
            seeing_list.append(2*np.sqrt(2*np.log(2))*seeing_sig) #assuming normal distribution for that
            polynomial_list.append(new_coeffs)
            print("Resetting the target_stack.")
            target_stack = [] #

        else:
            print("Next file is not a lamp.")
        
print(polynomial_list)
        #last_file_lamp= False #since this image isn't a lamp


last_file_lamp = False
target_stack = []
association_index = -1
new_filelist =[]
#need to determine if the given image is a lamp or a target spectrum
for counter, img in enumerate(speclist):
    filename= glob(img)[0]
    #if '_fe' in filename.lower():
    if '_fe.' in filename.lower():
        print('Lamp file detected: ', filename)
        
        if last_file_lamp:
            #since the previous file was a lamp, that would make this a new run, so we'd want to use this lamp file, right?
            print("Double lamp detected, so it must be a new run.")
            
            
        else:
            association_index+=1
            print("association_index increased: ", association_index)
            #last_file_lamp = True
            
        last_file_lamp = True #since the image has to be a lamp
        
    else:
        #the filename doesn't contain a lamp indicator, so it must be a target spectrum
        print("Target file detected: ", filename)
        
        i= fits.open(filename)
        header = fits.getheader(filename)
        img_data= np.copy(i[0].data)
        #As of 2024-03-05, the below line was running in the code.
        #img_data= img_data/normed_flat #also have to divide it here... since the other place was for seeing
        #2024-03-05 I'm changing it to be inside an if statement about skip_flat like it should have been the whole time
        if skip_flat:
            pass #don't do anything if you're skipping the flat
        else:
            img_data= img_data/normed_flat #also have to divide it here... since the other place was for seeing
        filename = 'w' + filename
        new_filelist.append(filename)
        polynomials = polynomial_list[association_index]
        seeing_sig = sigma_list[association_index]
        
        #20190624 changed
        
        core_sides, bkg_core_sides= seeing_window(seeing_sig)
        
        #end of 20190624 changed stuff
        
        #outsourced spectral extaction function
        #ctarget_light, cbkg_light, cnoise_spec= spext.extract_spectrum(img_data, header, polynomials=polynomials, core_sides= core_sides, bkg_core_sides=bkg_core_sides, bkg_shift=bkg_shift, bkg_method=bkg_method, bkg_poly_deg= bkg_poly, n_biases= n_biases)
        #print('\n\n\n=======\n')
        #print(filename)
        #print('bkg_core_sides that are input into the actual extraction', bkg_core_sides)
        #print('\n=========\n\n\n')
        target_light, bkg_light, noise_spectrum= spext.extract_spectrum(img_data, header, polynomials=polynomials, core_sides= core_sides, bkg_core_sides=bkg_core_sides, bkg_shift=bkg_shift, bkg_method=bkg_method, bkg_poly_deg= bkg_poly, n_biases= n_biases)
        
        ##2024-10-07 bug tracking
        #plt.plot(target_light,label="target_light")
        #plt.plot(bkg_light,label="bkg_light")
        #plt.title("in for loop back in wave_cal.py")
        #plt.legend()
        #plt.show()
        
        seeing_FWHM = seeing_list[association_index]
        band_inds= np.indices(img_data.shape)
        x_positions= band_inds[1,1]
        #target_light= np.array([])
        #target_noise2_list= np.array([])
        #bkg_light= np.array([])
        #bkg_noise2_list= np.array([])
        #bkg_up_comb=np.array([])
        #bkg_down_comb= np.array([])
        poly_curve_y = np.polyval(polynomials[0], x_positions)
        poly_curve_wavelength= np.polyval(polynomials[1], x_positions)
        upper_edges= np.polyval(polynomials[1], x_positions+0.5) #upper wavelength edges
        lower_edges= np.polyval(polynomials[1], x_positions-0.5) #lower_wavelength edges
        dlambda_vals= upper_edges-lower_edges
        
        #def bkg_trace(x_positions, sign='minus'):
            ##return  np.int_(poly_curve_y[x_positions]+bkg_shift)
            #if sign=='minus':
                #return  spt.discrete_int(poly_curve_y[x_positions])-bkg_shift
            #elif sign=='plus':
                #return spt.discrete_int(poly_curve_y[x_positions])+bkg_shift
        #for x_pos in x_positions:
            ##trace_vals=img_data[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides+1),x_pos]
            #trace_vals=img_data[spt.discrete_int(poly_curve_y[x_pos])-core_sides:spt.discrete_int(poly_curve_y[x_pos])+core_sides+1,x_pos]
            #xsum= np.sum(trace_vals)
            ##xsum= np.sum(img_data[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides+1),x_pos]) #old way 2018-10-31
            #target_light= np.append(target_light,[xsum])
            #target_noise2 = np.copy(xsum+trace_vals.shape[0]*header['RDNOISE']**2+trace_vals.shape[0]*header['RDNOISE']**2/n_biases)
            #target_noise2_list= np.append(target_noise2_list, [target_noise2])
            ##up_bkg=img_data[np.int_(poly_curve_y[x_pos]+bkg_shift-bkg_core_sides):np.int_(poly_curve_y[x_pos]+bkg_shift+bkg_core_sides+1),x_pos]
            ##bkg_up_coords= 
            #up_bkg=img_data[bkg_trace(x_pos, sign='plus')-bkg_core_sides:bkg_trace(x_pos, sign='plus')+bkg_core_sides+1,x_pos]
            ##down_bkg= img_data[np.int_(poly_curve_y[x_pos]-bkg_shift-bkg_core_sides):np.int_(poly_curve_y[x_pos]-bkg_shift+bkg_core_sides+1),x_pos]
            #down_bkg= img_data[bkg_trace(x_pos, sign='minus')-bkg_core_sides:bkg_trace(x_pos, sign='minus') +bkg_core_sides+1,x_pos]
            
            #plt.legend()
            #bkg_comb= np.append(up_bkg, down_bkg)
            ##print "trace_vals.shape", trace_vals.shape
            ##print "bkg_comb.shape", bkg_comb.shape
            #bkg_noise2= trace_vals.shape[0]*np.copy(np.mean(bkg_comb)/bkg_comb.shape[0]+header['RDNOISE']**2/bkg_comb.shape[0]+header['RDNOISE']**2/(bkg_comb.shape[0]*n_biases))
            ##bkg_sum= np.sum(img_data[np.int_(poly_curve_y[x_pos]+bkg_shift-core_sides):np.int_(poly_curve_y[x_pos]+bkg_shift+core_sides+1),x_pos])
            #bkg_sum= trace_vals.shape[0]*np.copy(np.mean(bkg_comb)) #take the mean of the bkg portion of the sky
            #bkg_light= np.append(bkg_light,[bkg_sum])
            #bkg_noise2_list= np.append(bkg_noise2_list, [bkg_noise2]) #list of noise values for a single pixel (resulting from the mean of the sky) for a given column
            #bkg_up_comb= np.append(bkg_up_comb, [np.sum(up_bkg)])
            #bkg_down_comb= np.append(bkg_down_comb, [np.sum(down_bkg)])
        ##plt.plot(bkg_up_comb, label='bkg_up_comb')
        ##plt.plot(bkg_down_comb, label='bkg_down_comb')
        ##plt.plot(target_light, label='target_light')
        ##plt.title('Before sky subtraction')
        ##plt.legend()
        ##plt.xlabel('pixel')
        ##plt.ylabel('counts')
        ##plt.show()
        ##plt.plot(x_positions,target_light,'-')
        ##plt.xlabel('x (pixel)')
        ##plt.ylabel('Counts')
        ##plt.title('Target Spectrum')
        ##plt.show()
        ##noise_spectrum = np.copy(np.sqrt(target_light + bkg_light + y_trace_width*header['RDNOISE'])) #old way 2018-10-31
        #noise_spectrum= np.copy(np.sqrt(target_noise2_list+bkg_noise2_list)) #combination of noises of the background pixels and the target trace pixels.
        #print "noise_spectrum.shape", noise_spectrum.shape
        #target_light= target_light-bkg_light
        #noise_spectrum = noise_spectrum/target_light #normalized noise values by the target spectrum, so now unitless.
        #target_light= target_light/header['EXPTIME'] #converting to counts/s
        #bkg_light= bkg_light/header['EXPTIME'] #converting to counts/s
        
        #plt.plot(target_light, label='target avg bkg')
        #plt.plot(ctarget_light, label='spext bkg '+bkg_method, linestyle='--')
        #plt.plot(bkg_light, label='avg bkg')
        #plt.plot(cbkg_light, label='spext bkg '+bkg_method, linestyle='--')
        #plt.legend()
        #plt.title('comparison of bkg methods')
        #plt.show()
        
        
        header.append(card=('bkgtype', bkg_method, 'method of bkg determination'))
        if bkg_method=='poly':
            header.append(card=('bkgdeg', bkg_poly, 'bkg polynomial degree'))
        else:
            pass
        
        
        if (do_airglow_corr and setup_dict['air_corr']) :
            #airline_array= np.genfromtxt(cp.line_list_dir+ cp.airline_name, names = True, delimiter='\t')
            airline_array= Table.read(cp.line_list_dir+cp.airline_name, format='ascii.tab')
            #print(airline_array)
            use_array=np.int_(airline_array['use'])
            good_airlines= np.copy(airline_array[np.where(use_array==1)])
            air_waves = np.float_(good_airlines['User'])
            #air_names= good_airlines['Name']+good_airlines['Name2']
            for air_wave, name, name2 in zip(air_waves, good_airlines['Name'], good_airlines['Name2']):
                #print(name+name2, type(name))
                air_name=name+name2
                plt.axvline(x=air_wave, linestyle='--', color=cp.airline_color[air_name[:2]])
                plt.text(air_wave, np.nanmax(bkg_light), air_name, color=cp.airline_color[air_name[:2]], rotation=90)
            print("need to put dlambda into this part again since we're about to move around wavelengths in the future.")
            plt.plot(poly_curve_wavelength, bkg_light)
            plt.xlabel(r'Wavelength ($\AA$)')
            plt.ylabel('Counts/s')
            plt.title('Sky with predicted line wavelengths marked')
            plt.show()
            coll_air_offsets= []
            coll_air_lam_offsets=[]
            for air_wave, name, name2 in zip(air_waves, good_airlines['Name'], good_airlines['Name2']):
                #print(name+name2, type(name))
                air_name=name+name2
                offset, offset_lambda= find_skyline_offset(x_positions, bkg_light, air_wave, polynomials[1], header, plot_all=False)
                coll_air_offsets.append(offset)
                coll_air_lam_offsets.append(offset_lambda)
            print('coll_air_offsets:', coll_air_offsets)
            print('coll_air_lam_offsets:', coll_air_lam_offsets)
            #avg_air_offset=np.mean(coll_air_offsets)
            if (air_off_type== 'lambda'):
                print('airglow offsetting will be done with wavelength zeropoint')
                avg_air_offset= np.median(coll_air_lam_offsets)
                poly_curve_wavelength= np.polyval(polynomials[1], x_positions)+avg_air_offset
                upper_edges= np.polyval(polynomials[1], x_positions+0.5)+avg_air_offset #upper wavelength edges
                lower_edges= np.polyval(polynomials[1], x_positions-0.5) +avg_air_offset#lower_wavelength edges
                dlambda_vals= upper_edges-lower_edges
                header.append(card=('airofstd', np.std(coll_air_offsets), 'std dev of sky emission angstrom offsets'))
                header.append(card=('airofavg', avg_air_offset, 'median of sky emission angstrom offsets, this was applied'))
            elif(air_off_type=='pixel'):
                avg_air_offset=np.median(coll_air_offsets)
                print('airglow offsetting will be done with pixel zeropoint')
                poly_curve_wavelength= np.polyval(polynomials[1], x_positions+avg_air_offset)
                upper_edges= np.polyval(polynomials[1], x_positions+avg_air_offset+0.5) #upper wavelength edges
                lower_edges= np.polyval(polynomials[1], x_positions+avg_air_offset-0.5) #lower_wavelength edges
                dlambda_vals= upper_edges-lower_edges
                header.append(card=('airofstd', np.std(coll_air_offsets), 'std dev of sky emission pixel offsets'))
                header.append(card=('airofavg', avg_air_offset, 'median of sky emission pixel offsets, this was applied'))
            else:
                print('\n\nYou seem to have not chosen a valid "air_off_type"; you used:', air_off_type)
                print('So no airglow correction is actually going to happen\n\n')
            header.append(card=('AIROFTYP', air_off_type, 'method of zeropoint correction for skylines'))
            for air_wave, name, name2 in zip(air_waves, good_airlines['Name'], good_airlines['Name2']):
                #print(name+name2, type(name))
                air_name=name+name2
                plt.axvline(x=air_wave, linestyle='--', color=cp.airline_color[air_name[:2]])
                plt.text(air_wave, np.nanmax(bkg_light), air_name, color=cp.airline_color[air_name[:2]], rotation=90)
            plt.plot(poly_curve_wavelength, bkg_light)
            plt.xlabel(r'Wavelength ($\AA$)')
            plt.ylabel('Counts/s')
            plt.title('Sky with Correction of'+ str(avg_air_offset)+air_off_type+'  applied')
            plt.show()
            
        else:
            pass
        
        
        plt.plot(x_positions,target_light,'-')
        plt.xlabel('x (pixel)')
        #plt.ylabel('Counts')
        plt.ylabel('Counts/s')
        plt.title('Target Spectrum after background subtraction')
        plt.show()
        header= to_barycenter(header) #append the BMJD_TDB value
        header.append(card= ("pix_scal", 0.3, ' "/pixel'))
        header.append(card = ('see_sig', seeing_sig, 'Sigma of Gauss seeing fit (pixels)'))
        header.append(card = ('see_FWHM', seeing_FWHM, 'Seeing (pixels)'))
        header.append(card = ('skipflat', skip_flat, 'flatfielding skipped or not'))
        header.append(card=('trc_mthd',trace_method, 'trace_method Method of identifying trace'))
        header.append(card=('trc_off', trace_offset, 'offset of ext trace from fit'))
        header.append(card=('width', core_sides*2+1, 'width of extracted region for trace'))
        header.append(card=('bkgwidth', bkg_core_sides*2+1, 'width of bkg regions'))
        header.append(card=('bkgshift', bkg_shift, 'shift of bkg region from center of trace'))
        header.append(card=('waveoff',store_offset,'pix offset input for lamp file'))
        
        if do_save_wavesoln:
            wave_soln_name=save_wavesoln(filename, polynomials[1])
            header.append(card=('wavesoln', wave_soln_name, 'filename of wavelength solution poly coeffs'))
        else:
            pass
        if do_save_tracetrack:
            tracetrack_name=save_tracetrack(filename, polynomials[0])
            header.append(card=('trctrack', tracetrack_name, 'filename of trace tracking poly coeffs'))
        else:
            pass
        #poly_curve_wavelength= barycentric_vel_corr(header, poly_curve_wavelength) #correction of Earth's orbital motion
        header['units']= 'Counts/s'
        
        ###2024-10-07 bug tracking
        ##plt.plot(target_light,label="target_light")
        ##plt.plot(bkg_light,label="bkg_light")
        ##plt.title("right before writing to fits file")
        ##plt.legend()
        ##plt.show()
        
        hdu = fits.PrimaryHDU(poly_curve_wavelength, header = header)
        hdu1= fits.ImageHDU(target_light)
        hdu2= fits.ImageHDU(bkg_light)
        hdu3 = fits.ImageHDU(noise_spectrum)
        hdu4=fits.ImageHDU(dlambda_vals)
        hdulist= fits.HDUList([hdu, hdu1, hdu2, hdu3, hdu4])
        hdulist.writeto(filename, overwrite= True)
        #target_stack.append(img_data)
        last_file_lamp = False
        if '_fe' in speclist[counter+1].lower():
            print("Next file is a lamp")
            pass
            #print "Next file is a lamp, so we're going to do the trace and wavelength calibration."
            

        else:
            print("Next file is not a lamp.")
            
            
            
        
