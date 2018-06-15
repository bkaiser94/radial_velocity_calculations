"""
This should read in the raw image files for the spectra and the iron lamp that corresponds.

listComb should only include the spectra for one standard star at a time.
listFe is the reference lamp spectrum you want to use. IT HAS TO HAVE THE SAME BINNING AS THE
STANDARD SPECTRUM!

This is intended to produce the 1-d extracted spectra that are also median-combined  for a given standard star to be used in flux calibration later.

Outputs are used by flux_calibration.py to produce sensitivity curves

Outputs wcmtb.*       named files.

STEP 4 of Reduction

PART OF STEP-BY-STEP REDUCTION
"""

#I need something to indicate that this is supposed to be a standard star or the actual target in order to fix the 


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


c= 2.998E5  #km/s
zerolistname= 'listZero'
combinelistname = 'listComb'
#flatlistname = 'listFlat'
#speclistname= 'listSpec'
lamplistname = 'listFe'
#output_filename = 'radial_velocities.txt'


linefilename = 'JJ_FeAr_lines.txt'
#linefilename = 'FeAr_3650to5250_lines_GOODMAN.txt'
#flatlist= np.genfromtxt(flatlistname,dtype = 'str' )
zerolist = np.genfromtxt(zerolistname, dtype ='str')
combinelist = np.genfromtxt(combinelistname, dtype = 'str')
lamplist = np.genfromtxt(lamplistname, dtype = 'str')
#speclist= np.genfromtxt(speclistname, dtype = 'str')

#########3
slit_ystart = 1   #The beginning of the image that has light from outside
slit_yend= 199     #The end of the image with same
trace_xstart = 9
trace_xend = 2055
#bkg_width= 10   #How many pixel rows should be sampled on each edge of the slit
trace_band_mid= 100   #y-pixel that's about the center of the bulge of the galaxy
trace_band_width = 20 #pixel width to determine the centroid of the galaxy
poly_degree= 3
core_sides=  5
bkg_width= core_sides
bkg_shift= 25


lamp_sigma_guess = 2
line_search_width= 3
balmer_sigma_guess= 14
lamp_p0 = [100, 500,  lamp_sigma_guess, 0]
#balmer_p0= [-1, 500, balmer_sigma_guess,balmer_line_sides[0],0]
balmer_p0= [-100, 500, balmer_sigma_guess,0]


fear_array= np.genfromtxt(linefilename, names = True)
line_x_checks = np.copy(fear_array['Pixel']) +90
lamp_lines = np.copy(fear_array['User'])
line_sides = np.ones(line_x_checks.shape[0])*line_search_width

def make_image_stack(imagelist, times= True):
    """
    
    """
    images = []
    timestamps = []
    exptimes = []
    expstarts = []
    filename1 = glob(imagelist[0])
    image=fits.open(filename1[0])
    header1 = fits.getheader(filename1[0])
    copy_header = fits.Header.copy(header1)
    for img in imagelist:
        filename = glob(img)[0]
        i= fits.open(img)
        header = fits.getheader(img)
        img_data= i[0].data
        images.append(img_data)
        gain =header['GAIN']
        readnoise = header['RDNOISE']
        if times:
            #starttime = header['OPENTIME']
            #startdate = header['OPENDATE']
            #expstart = Time(str(startdate)+'T'+str(starttime), format = 'isot', scale = 'utc')
            #expstart = str(startdate)+'T' +str(starttime)
            expstart = header['DATE-OBS']
            expstarts.append(expstart)
            exptime= header['EXPTIME']
            exptimes.append(exptime)
            #timestamp = expstart+exptime/2
            #print timestamp
            #timestamps.append(timestamp)
    if times:
        expstarts = Time(expstarts, format = 'isot', scale = 'utc')
        timestamps = expstarts + np.array(exptimes)*u.s/2.
        
        
    return np.array(images),gain, readnoise, timestamps,copy_header
#################

bias_stack,gain, readnoise, bias_times, bias_header = make_image_stack(zerolist, times= False)
bias_med = np.nanmedian(bias_stack, axis=0)



lampfile = glob(np.array([lamplist])[0])[0]
lamp_hdu = fits.open(lampfile)
lamp_im = np.copy(lamp_hdu[0].data)

target_stack, gain, readnoise, target_times, target_header = make_image_stack(combinelist)

### bias subtraction #####
target_stack = target_stack - bias_med
lamp_im = lamp_im - bias_med

#### Trimming images ####
target_stack = target_stack[:, :, trace_xstart:trace_xend]
lamp_im = lamp_im[:, trace_xstart:trace_xend]


#### cosmic ray subtraction #### 
for target_frame in target_stack:
    target_cosmic= cosmics.cosmicsimage(target_frame, gain=gain, readnoise=readnoise, sigclip = 5.0, sigfrac = 0.3, objlim = 5.0)
    target_cosmic.run(maxiter= 4)
    target_frame= np.copy(target_cosmic.cleanarray)
    
### convert to electron counts
target_stack = target_stack * gain
    

##### Median combining target frames####

print target_stack.shape
target_med = np.median(target_stack, axis =0)
print target_med.shape

##### polynomial fitting #### 

target_band= target_med[trace_band_mid-trace_band_width/2:trace_band_mid+trace_band_width/2,:]
print 'target_band.shape',target_band.shape
band_inds= np.indices(target_band.shape)
x_positions= band_inds[1,1]
print 'xpositionsshape', x_positions.shape
y_positions= np.argmax(target_band,axis=0)+(trace_band_mid-trace_band_width/2)
print 'yshape', y_positions.shape
polynomial_fit= np.polyfit(x_positions,y_positions,poly_degree)
print polynomial_fit
print polynomial_fit.shape
poly_curve_y= polynomial_fit[3]+polynomial_fit[2]*x_positions + polynomial_fit[1]*(x_positions**2)+ polynomial_fit[0]*(x_positions**3)

plt.imshow(np.log(target_med),cmap = 'hot', interpolation = 'none')
plt.plot(x_positions, poly_curve_y, color = 'blue', label  = 'polynomial fit')
plt.plot(x_positions,y_positions, color = 'black', label = 'max values', linestyle = 'none', marker = '*')
plt.legend()
plt.show()


target_light= np.array([])
bkg_light= np.array([])
lamp_light= np.array([])
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

####################
### background subtraction ####

target_light= target_light-bkg_light


##### Wavelength calibration ####


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


def gaussian_curve(x, a, x0, sigma,b):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+b

def fit_gaussian_curve(x_pixels, light_values, p0_list, search_width):
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


peaks_found=[]
wave_peaks_found = []
for lamp_line_guess,lamp_line_wave in zip( line_x_checks,lamp_lines):
    try:
        lamp_params, lamp_cov = fit_gaussian_curve(x_positions, lamp_light, [lamp_p0[0], lamp_line_guess, lamp_p0[2], lamp_p0[3]], line_search_width)
        if ((np.abs(lamp_params[0]) > 1.) and (np.abs(lamp_params[2])< 20) and (lamp_params[0] > 0) and (np.abs(lamp_line_guess-lamp_params[1]) < line_search_width) ):
            peaks_found.append(lamp_params[1])
            wave_peaks_found.append(lamp_line_wave)
            plt.plot(x_positions, lamp_light, label = 'lamp data', color = 'blue')
            plt.plot(x_positions, gaussian_curve(x_positions, lamp_params[0], lamp_params[1], lamp_params[2], lamp_params[3]), color = 'r', label = 'Gaussian Fit')
            plt.title("guess: " + str(lamp_line_guess) + ' fit:' + str(lamp_params[1]))
            for x_spot in line_x_checks:
                plt.axvline( x= x_spot, color = 'k',linestyle = '--')
            plt.axvline(x = lamp_line_guess, color = 'r', linestyle= '--')
            plt.legend()
            plt.show()
        else:
            print "Gaussian too flat or flipped (or not within the actual fitting region...):", lamp_params
    except RuntimeError as error:
        print error

######## Line locations
peaks_found = np.array(peaks_found)
wave_peaks_found = np.array(wave_peaks_found)
print "line_x_checks:"
print line_x_checks
print "peaks found"
print peaks_found
print "wave_peaks_found"
print wave_peaks_found
for line,peak,wave in zip(line_x_checks, peaks_found, wave_peaks_found):
    print line, peak, wave
#polynomial fitting
#poly_coeffs_lamp= np.polyfit(centroids,lamp_lines,2)
poly_coeffs_lamp =np.polyfit(peaks_found, wave_peaks_found, 5)


def x_to_wavelength(x_positions):
	#poly_curve_wavelength= poly_coeffs_lamp[2]+poly_coeffs_lamp[1]*x_positions + poly_coeffs_lamp[0]*(x_positions**2)
	poly_curve_wavelength= poly_coeffs_lamp[-1]+poly_coeffs_lamp[-2]*x_positions + poly_coeffs_lamp[-3]*(x_positions**2)+poly_coeffs_lamp[-4]*(x_positions**3)+poly_coeffs_lamp[-5]*(x_positions**4)+poly_coeffs_lamp[-6]*(x_positions**5)
	return poly_curve_wavelength
poly_curve_wavelength= x_to_wavelength(x_positions)


print poly_curve_wavelength.shape
print poly_coeffs_lamp
print "smallest wavelength on plot: ",poly_curve_wavelength[0]
print "largest wavelength on plot: ", poly_curve_wavelength[-1]

#### diagnostic plots #### 
plt.plot(poly_curve_wavelength,target_light,'-')
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Counts')
plt.title('Target Spectrum')
#plt.yscale('log')
#plt.ylim(10,200)
plt.show()

plt.plot(x_positions, poly_curve_wavelength,  label = 'wavelength solution', color ='blue')
plt.plot(peaks_found, wave_peaks_found, marker= '*', linestyle = 'none', label = 'fitted values', color = 'red' )
plt.plot(line_x_checks, lamp_lines, label = 'input points', color = 'green', marker = '*', linestyle = 'none')
plt.title("wavelength to pixel position")
plt.legend()
plt.show()

plt.axhline(y=0 ,  label = 'wavelength solution', color ='blue')
plt.plot(peaks_found, wave_peaks_found-x_to_wavelength(peaks_found), marker= '*', linestyle = 'none', label = 'fitted values', color = 'red' )
plt.plot(line_x_checks, lamp_lines-x_to_wavelength(line_x_checks), label = 'input points', color = 'green', marker = '*', linestyle = 'none')
plt.title("wavelength to pixel position Residuals")
plt.legend()
plt.show()


######## outputting the reduced 1-d spectra
#0 index should be the wavelength solution for the file
# 1 index should be the reduced spectrum
#2 index should be the sky spectrum


output_filename = ''.join(combinelist[0].split('_')[1:])
output_filename = "wcmtb." + output_filename #wavelength-calibrated, cosmic-subtracted, median-combined, trimmed, bias-subtracted
#output_array = np.vstack([poly_curve_wavelength, target_light, bkg_light])
#print output_array.shape
hdu = fits.PrimaryHDU(poly_curve_wavelength, header = target_header)
hdu1 = fits.ImageHDU(target_light)
hdu2= fits.ImageHDU(bkg_light)
hdulist = fits.HDUList([hdu, hdu1, hdu2])
hdulist.writeto(output_filename, overwrite = True)
#fits.writeto(output_filename, output_array, header= target_header, overwrite = True)

