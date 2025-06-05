"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-06-03

This should take a 1x2 binned image and rebin it to be a 2x2 so it can be rammed through other code.


"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
#import cosmics
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
import scipy.interpolate as scinterp
import os

imlistname = 'listRebin'

imlist= np.genfromtxt(imlistname, dtype = 'str')


def rebin_image(input_file):
    header=fits.getheader(input_file)
    hdu=fits.open(input_file)
    image=hdu[0].data
    print('shape',image.shape)
    rowlist=[]
    for index in np.arange(0,image.shape[1],2):
        added=image[:,index]+image[:,index+1]
        rowlist.append(added)
    
    out_image=np.array(rowlist).T
    #plt.title('Original Image')
    #plt.imshow(np.log10(image))
    #plt.show()
    #plt.title('binned image')
    #plt.imshow(np.log10(out_image))
    #plt.show()
    
    new_name=input_file.split('.')[0]+'_rebin.fits'
    new_hdu=fits.PrimaryHDU(out_image,header=header)
    new_hdu.writeto(new_name)
    return


for filename in imlist:
    filename=glob(filename)[0]
    rebin_image(filename)
    
