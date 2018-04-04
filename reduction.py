"""
Written by Ben Kaiser (UNC - Chapel Hill), with some additions of code from the red_cam_pipeline,
which contains code from the ZZCeti_pipeline from Josh Fuchs github , which was written by J Meza, updated by J Fuchs, and further updated by Ben Kaiser
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


c= 2.998E5  #km/s
zerolistname= 'listZero'
flatlistname = 'listFlat'
speclistname= 'listSpec'

output_filename = 'radial_velocities.txt'


#linefilename = 'FeAr_3650to5250_lines.txt'
linefilename = 'JJ_FeAr_lines.txt'
#linefilename = 'FeAr_3650to5250_lines_GOODMAN.txt'
flatlist= np.genfromtxt(flatlistname,dtype = 'str' )
zerolist = np.genfromtxt(zerolistname, dtype ='str')
speclist= np.genfromtxt(speclistname, dtype = 'str')

#########3
slit_ystart = 1   #The beginning of the image that has light from outside
slit_yend= 199     #The end of the image with same
trace_xstart = 9
trace_xend = 2055
#bkg_width= 10   #How many pixel rows should be sampled on each edge of the slit
trace_band_mid= 105   #y-pixel that's about the center of the bulge of the galaxy
trace_band_width = 20 #pixel width to determine the centroid of the galaxy
poly_degree= 3
core_sides=  5
bkg_width= core_sides
bkg_shift= 50

#lamp_lines= np.array([4358.3, 5460.7, 5769., 5790.6, 6965.4]) #HgAr wavelengths in angstroms
#line_sides= np.array([4,4,4,4,4])
#line_x_checks= np.array([175, 728, 883, 893, 1490]) #x-values that are approximately in the middle of the lamp lines so we know where to look for centroiding
balmer_rest_waves= np.array([4101.734, 4340.472, 4861.35]) #Only 3 Balmer Lines seem to be in this image.
balmer_x_checks= np.array([610 ,890, 1510])   #wavelength values that should be checked for the Balmer emission lines in our spectrum (1360 line should be 30 pixels on either side)
#balmer_line_sides= np.array([30,30, 31]) #Number of pixels on either side of the guessed peak that should be included in the centroiding effort for the balmer lines (those are based on attempted lines).
balmer_line_sides = 30.
###
lamp_sigma_guess = 3
line_search_width= 5
balmer_sigma_guess= 14
lamp_p0 = [100, 500,  lamp_sigma_guess, 0]
#balmer_p0= [-1, 500, balmer_sigma_guess,balmer_line_sides[0],0]
balmer_p0= [-100, 500, balmer_sigma_guess,0]

###
n_trials= 1000
########
#iron lines from red_cam_pipeline

WaveList_Fe_930_12_24= np.array([ [431.795, 1057.76, 1194.97, 1315.35, 
                                   1381.1, 1444.58,  1630.61, 
                                   1682.99, 1726.03, 1779.61, 
                                   1893.19, 2132.53, 2210.85, 
                                   2279.64, 2361.34, 2443.06, 2468.22, 
                                   2515.14, 2630.53, 2795.45, 2886.45, 
                                   2985.15, 3085.52, 3162.56,  
                                   3367.86, 3795.76, 3845.65, 
                                   3907.57], 
                                   
                                   [3729.3087, 3994.7918, 4052.9208, 4103.9121, 
                                    4131.7235, 4158.5905,  4237.2198, 
                                    4259.3619, 4277.5282, 4300.1008, 
                                    4348.064, 4448.8792, 4481.8107, 
                                    4510.7332, 4545.0519, 4579.3495, 4589.8978, 
                                    4609.5673, 4657.9012, 4726.8683, 4764.8646, 
                                    4806.0205, 4847.8095, 4879.8635, 
                                    4965.0795, 5141.7827, 5162.2846, 
                                    5187.7462] ]) 
                                   
#WaveList_Fe_930_12_24= np.array([ [431.795, 1057.76,  
                                   #1444.58,  1630.61, 
                                   #1726.03, 
                                   #2132.53,  
                                   #2279.64, 2361.34,  
                                   #2886.45, 
                                   #2985.15, 3085.52,   
                                   #3367.86, 3795.76, 3845.65, 
                                   #3907.57], 
                                   
                                   #[3729.3087, 3994.7918,   
                                    #4158.5905,  4237.2198, 
                                    #4277.5282, 
                                    #4448.8792,  
                                    #4510.7332, 4545.0519, 
                                    #4764.8646, 
                                    #4806.0205, 4847.8095,  
                                    #4965.0795, 5141.7827, 5162.2846, 
                                    #5187.7462] ]) 
#line_x_checks = WaveList_Fe_930_12_24[0]/2.
#lamp_lines = WaveList_Fe_930_12_24[1]
#line_sides = np.ones(line_x_checks.shape[0])*4.

#fear_array = np.genfromtxt(linefilename, names= True)
#line_x_checks = np.int_(fear_array['Pixel']/8000.*(trace_xend-trace_xstart)) #rescaled to match the dimensions of GOODMAN
#lamp_lines = fear_array['Fit']
#good_fits = np.where(np.abs(fear_array['Residual'] < 0.01))
#line_x_checks = line_x_checks[good_fits]
#lamp_lines= lamp_lines[good_fits]
#line_sides = np.ones(line_x_checks.shape[0])*4.

#fear_array = np.genfromtxt(linefilename, names=True)
#all_line_x_checks = np.float_(fear_array['GOODMAN_Pixel'])
#good_fits = np.where(fear_array['Use']==True)
#print fear_array['Use']
#line_x_checks = np.copy(all_line_x_checks[good_fits])
#lamp_lines= np.copy(fear_array['User'][good_fits])
#line_sides= np.ones(line_x_checks.shape[0])*line_search_width

#For JJ's line list
fear_array= np.genfromtxt(linefilename, names = True)
line_x_checks = np.copy(fear_array['Pixel']) +91
lamp_lines = np.copy(fear_array['User'])
line_sides = np.ones(line_x_checks.shape[0])*4.


def make_image_stack(imagelist, times= True):
    """
    
    """
    images = []
    timestamps = []
    exptimes = []
    expstarts = []
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
        timestamps = expstarts + np.array(exptimes)*u.s
        
        
    return np.array(images),gain, readnoise, timestamps



###############




bias_stack,gain, readnoise, bias_times = make_image_stack(zerolist, times= False)
bias_med = np.nanmedian(bias_stack, axis=0)

#flat_stack, gain, readnoise, flat_times = make_image_stack(flatlist)
#flat_stack = flat_stack - bias_med #bias subtraction
#flat_med = np.nanmedian(flat_stack, axis=0)
#flat_norm = flat_med/np.nanmedian(flat_med)


target_stack, gain, readnoise, target_times = make_image_stack(speclist)
print target_stack.shape
target_stack = target_stack - bias_med #bias subtraction
#target_stack = target_stack/flat_norm

target_stack = target_stack[:, :, trace_xstart:trace_xend]

first_lamp =np.copy(target_stack[0])
first_lamp_time = target_times[0]
end_lamp = np.copy(target_stack[-1])
end_lamp_time = target_times[-1]
target_stack= target_stack[1:-1] #removing the lamps
target_times= target_times[1:-1]
print target_stack.shape

#lamp_im = np.nanmean([first_lamp, end_lamp], axis=0)
lamp_time = first_lamp_time
lamp_im = first_lamp
#lamp_im= end_lamp
#lamp_time = end_lamp_time
target_med = np.nanmedian(target_stack, axis =0)



##########################Polynomial Curve Fitting
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

#plt.plot(x_positions,target_light,'-')
#plt.xlabel('x (pixel)')
#plt.ylabel('Counts')
#plt.title('Target Spectrum with Sky Background Subtraction')
#plt.show()

plt.plot(x_positions,target_light,'-')
plt.xlabel('x (pixel)')
plt.ylabel('Counts')
plt.title('DO NOT EXIT THIS PLOT UNTIL YOU HAVE INPUT ALL 3 OF THE RAW INPUT VALUES FOR THE BALMER LINES')
plt.show()

Balmer2= float(raw_input("Rightmost Balmer line center pixel>>>"))
Balmer1= float(raw_input("2nd Rightmost Balmer line center pixel>>>"))
Balmer0= float(raw_input("3rd Rightmost Balmer line center pixel>>>"))

balmer_x_checks=np.array([Balmer0,Balmer1,Balmer2])


#plt.plot(x_positions,bkg_light,'-')
#plt.xlabel('x (pixel)')
#plt.ylabel('Counts')
#plt.title('Sky Background Spectrum')
#plt.show()
####

def impose_floor(input_light, floor_cutoff):
    """
    basically shifts the whole spectrum down in brightness by a constant value and then sets all negative values to zero
    """
    output_light = np.copy(input_light)-floor_cutoff
    ignore = np.where(output_light < 0.)
    output_light[ignore]= 0.
    return output_light
######

#lamp_light= impose_floor(lamp_light, 0.125*np.max(lamp_light))


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

low_values = np.where(np.abs(x_positions-500)<300)
plt.plot(x_positions[low_values],lamp_light[low_values],linestyle= '-')
#plt.ylim((0,1000))
plt.xlabel('x (pixel)')
plt.ylabel('Counts')
plt.title('Lamp Spectrum (record corresponding dotted line and emission pixels)')
plt.yscale('log')
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
####
####

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



def line_centroiding(lamp_light,line_x_checks,line_sides):
	centroids= np.array([])
	for i in range(0,line_x_checks.shape[0]):
		xmin= line_x_checks[i]-line_sides[i]
		xmax= line_x_checks[i]+line_sides[i]+1
		line_centroid= (np.sum(x_positions[xmin:xmax]*lamp_light[xmin:xmax])/np.float_(np.sum(lamp_light[xmin:xmax])))
		#Running the centroiding a second time to make the fit slightly better.
		xmin= np.int_(line_centroid-line_sides[i])
		xmax= np.int_(line_centroid+line_sides[i]+1)
		line_centroid= (np.sum(x_positions[xmin:xmax]*lamp_light[xmin:xmax])/np.float_(np.sum(lamp_light[xmin:xmax])))
		centroids= np.append(centroids,[line_centroid])
	return centroids

#centroids = line_centroiding(lamp_light, line_x_checks,line_sides)
peaks_found=[]
wave_peaks_found = []
for lamp_line_guess,lamp_line_wave in zip( line_x_checks,lamp_lines):
    try:
        lamp_params, lamp_cov = fit_gaussian_curve(x_positions, lamp_light, [lamp_p0[0], lamp_line_guess, lamp_p0[2], lamp_p0[3]], line_search_width)
        if ((np.abs(lamp_params[0]) > 1.) and (np.abs(lamp_params[2])< 20) and (lamp_params[0] > 0) ):
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
            print "Gaussian too flat or flipped:", lamp_params[0], lamp_params[2]
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


#print "test x-positions: ", line_x_checks
#print "centroids (x-positions): ", centroids
#print "Wavelengths (anstroms): ", lamp_lines

print poly_curve_wavelength.shape
print poly_coeffs_lamp
print "smallest wavelength on plot: ",poly_curve_wavelength[0]
print "largest wavelength on plot: ", poly_curve_wavelength[-1]

print "Guessed Line Wavelengths: ", x_to_wavelength(balmer_x_checks)
print "Balmer Rest wavelengths (anstroms): ", balmer_rest_waves


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

def get_redshift(rest_lam, obs_lam):
    return (obs_lam-rest_lam)/rest_lam

def get_radial_velocity(redshift):
    return redshift*c
#balmer_centroids= line_centroiding(target_light, balmer_x_checks, balmer_line_sides)

rv_list= []
rv_low=[]
rv_high = []
for target_frame in target_stack:
    target_cosmic= cosmics.cosmicsimage(target_frame, gain=gain, readnoise=readnoise, sigclip = 5.0, sigfrac = 0.3, objlim = 5.0)
    target_cosmic.run(maxiter= 4)
    target_frame= target_cosmic.cleanarray
    target_light= np.array([])
    bkg_light= np.array([])
    print target_light.shape
    for x_pos in x_positions:
        xsum= np.sum(target_frame[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides),x_pos])
        target_light= np.append(target_light,[xsum])
        bkg_sum= np.sum(target_frame[np.int_(poly_curve_y[x_pos]+bkg_shift-core_sides):np.int_(poly_curve_y[x_pos]+bkg_shift+core_sides),x_pos])
        bkg_light= np.append(bkg_light,[bkg_sum])
    #plt.plot(x_positions,target_light,'-')
    #plt.xlabel('x (pixel)')
    #plt.ylabel('Counts')
    #plt.title('Target Spectrum')
    #plt.show()
    balmer_centers = []
    balmer_sigmas = []
    for balmer_line_x, balmer_wave in zip(balmer_x_checks, balmer_rest_waves):
        try:
            balmer_params, balmer_cov = fit_gaussian_curve(x_positions, target_light, [balmer_p0[0], balmer_line_x, balmer_p0[2], balmer_p0[3]], balmer_line_sides)
            balmer_centers.append(balmer_params[1])
            balmer_sigmas.append(balmer_params[2])
            plt.plot(x_positions, target_light, label = 'observed', color = 'blue')
            plt.plot(x_positions, gaussian_curve(x_positions,balmer_params[0], balmer_params[1], balmer_params[2], balmer_params[3]), color = 'r', label = 'Gaussian Fit')
            plt.title("guess: " + str(balmer_line_x) + ' fit:' + str(balmer_params[1])+' restwave:' + str(balmer_wave))
            plt.legend()
            plt.show()
        except RuntimeError as error:
            print error
    balmer_centers = np.array(balmer_centers)
    balmer_sigmas= np.array(balmer_sigmas)
    #balmer_centroids_waves= x_to_wavelength(balmer_centroids)
    balmer_sigma_up = x_to_wavelength(np.copy(balmer_centers+balmer_sigmas))
    balmer_sigma_down = x_to_wavelength(np.copy(balmer_centers-balmer_sigmas))
    balmer_centers = x_to_wavelength(balmer_centers)
    print "Balmer_lines: ", balmer_rest_waves
    print "Target_Balmer_lines", balmer_centers
    #print "Target_Balmer_lines ", balmer_centroids_waves
    #redshifts = get_redshift(balmer_rest_waves, balmer_centroids_waves)
    redshifts= get_redshift(balmer_rest_waves, balmer_centers)
    print "redshifts: ", redshifts
    radial_velocities = get_radial_velocity(redshifts)
    balmer_sigma_down= get_radial_velocity(get_redshift(balmer_rest_waves, balmer_sigma_down))
    balmer_sigma_up= get_radial_velocity(get_redshift(balmer_rest_waves, balmer_sigma_up))
    print "radial_velocities:", radial_velocities
    print np.mean(radial_velocities)
    rv_list.append(radial_velocities)
    rv_low.append(balmer_sigma_down)
    rv_high.append(balmer_sigma_up)
    print '-------------------'

for rv1, rv_val, rv2 in zip(rv_low, rv_list, rv_high):
    print rv1, "|", rv_val, "|", rv2

rv_list = np.array(rv_list).T
target_times= target_times.mjd
target_times = np.array([target_times])
print target_times.shape
print rv_list.shape
print target_times.dtype
print rv_list.dtype
comb_array = np.vstack([target_times, rv_list])
print comb_array.shape
response = raw_input("Should the data be output?>>> ")
response_num = int(raw_input("What number should be tacked onto the end of the filename?>>>"))
if response.startswith('y'):
    #try:
        #previous_array = np.genfromtxt(output_filename)
        #print "====="
        #print previous_array.shape
        #print comb_array.T.shape
        #new_comb_array = np.copy(np.vstack([previous_array,comb_array.T]))
        #print new_comb_array.shape
        #print new_comb_array
        #np.savetxt(output_filename, new_comb_array)
    #except IOError as error:
        #print error
        #print comb_array.shape
    split_out = output_filename.split('.')
    output_filename= split_out[0]+str(response_num)+'.'+split_out[1]
    np.savetxt(output_filename, comb_array.T)
else:
    print "As you wish. Data will not be saved...exiting"

