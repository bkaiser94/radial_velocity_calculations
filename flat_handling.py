"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-11-04

This script should take all of the flats and generate an overall flat file that should be used to correct the individual exposures for low level variations in the pixels.

straight_reduction.py should already have done the bias-subtraction, cosmic ray removal, and trimming of the overscan regions.

Need to fit a polynomial to the underlying spectral features of the lamp to remove them and leave behind only the low-level variations of the instrument itself.

THIS IS ACTUALLY STEP 2-B FOR THE REDUCTION PROCESS. THIS NEEDS TO BE EXECUTED AFTER straight_reduction.py and before wave_cal.py
"""

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from astropy.io import fits
import scipy.interpolate as scinterp


#flatlistname='listFlat'
#zerolistname= 'listZero'

#flatlist = np.genfromtxt(flatlistname, dtype ='str')

plot_all= True
master_flat_name= 'mctb.master_flat.fits'

poly_degree= 50

def make_image_stack(imagelist):
    """
    
    """
    images = []
    print imagelist
    for img in imagelist:
        filename = glob(img)[0]
        i= fits.open(img)
        header = fits.getheader(img)
        img_data= i[0].data
        images.append(img_data)
        readnoise = header['RDNOISE']
        
        
    return np.array(images), readnoise, header



def normalize_flat(input_flat):
    single_projection= np.nanmean(input_flat, axis=0)
    print('single_projection.shape', single_projection.shape)
    x_positions=np.arange(0, single_projection.shape[0])
    flat_poly_coeffs= np.polyfit( x_positions, single_projection, poly_degree)
    flat_poly_vals= np.polyval(flat_poly_coeffs,x_positions)
    
    cs = scinterp.CubicSpline(x_positions, single_projection)
    
    
    plt.plot(single_projection, label='single_projection')
    plt.plot(flat_poly_vals, label='flat_poly_vals')
    plt.plot(x_positions, cs(x_positions), label='CubicSpline')
    plt.title('Flat projection and polynomial of degree ' + str(poly_degree))
    plt.legend()
    plt.show()
    
    plt.plot(single_projection/flat_poly_vals, label='divide by polynomial')
    plt.plot(single_projection/cs(x_positions), label='divide by CubicSpline')
    plt.title('Divided residuals of single_projection')
    plt.legend()
    plt.show()
    return



def make_masterflat(flatlistname):
    """
    read in all of the flat images from a list. These should already be bias-subtracted and cosmic-ray removed.
    
    produce an averaged version of the flat that isn't normalized
    
    produce a corresponding frame of the uncertainties for each pixel
    
    """
    flatlist= np.genfromtxt(flatlistname, dtype= 'str')
    
    flatstack, readnoise, header = make_image_stack(flatlist)
    #flatstack= flatstack*gain #straight_reduction.py already outputs images in counts, so this is unnecessary as a step
    flatsum= np.sum(flatstack, axis=0) #summed the flats across all frames
    flatstack_err2= flatstack+ readnoise**2 #creating the variance version of the flatstack
    flatsum_err= np.copy(np.sqrt(np.sum(flatstack_err2,axis=0))) #summed errors
    frame_sums = np.sum(np.sum(flatstack, axis=-1),axis=-1)
    if plot_all:
        projection = flatsum[100]
        plt.plot(projection)
        plt.show()
        projection_err= flatsum_err[100]
        plt.plot(projection_err)
        plt.show()
        plt.plot(frame_sums, marker = 'o', linestyle = 'none')
        plt.show()
    else:
        pass
    flatsum_err= flatsum_err/flatsum #normalized noise values.
    header.append(card = ('sum_ext', 0, 'extension of summed flat'))
    header.append(card= ('err_ext', 1, 'extension of scaled sigma values'))
    hdu= fits.PrimaryHDU(flatsum, header= header)
    hdu1=fits.ImageHDU(flatsum_err)
    hdulist= fits.HDUList([hdu,hdu1])
    hdulist.writeto(master_flat_name, overwrite= True)
    
    normalize_flat(flatsum)
    return flatsum



flatsum= make_masterflat('listCTBflat')
