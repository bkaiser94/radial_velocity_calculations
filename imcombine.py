"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-09-19


Should be a relatively simple script to combine 2-d images together and output the result as a FITS file.
Files should be read from some sort of input filelist which I think I'll make be named listIM. But that does have the
drawback that i's don't distinguish from L's when capital in many fonts, including the one in the editor I'm using 
right now...


"""
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
import cosmics
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const
from astropy.modeling import models as asmodels
from astropy.modeling import fitting as asfitting
from astropy.table import Table, Column


imlist_name= 'listIM'


imlist= np.genfromtxt(imlist_name, dtype='str')

im_stack=[]

for im_name in imlist:
    i= fits.open(im_name)
    header = fits.getheader(im_name)
    img_data= np.copy(i[0].data)
    im_stack.append(img_data)
    
    
med_im= np.nanmedian(im_stack, axis=0)
print('median image generated')

og_name= imlist[0]

dec_parts = og_name.split('.')
mid_sub_parts= dec_parts[1].split('_')

mid_new_name= '_'.join(mid_sub_parts[1:])
tail_end= ''.join(dec_parts[2:])
full_new_name= '.'.join(['med_'+dec_parts[0], mid_new_name, tail_end])

print('full_new_name:', full_new_name)

plt.imshow(med_im)
plt.show()

print('saving image')
hdu= fits.PrimaryHDU(med_im, header=header)
hdu.writeto(full_new_name, overwrite=True)
print('image saved.')

