"""
This script should create a master Zero frame, which is used for bias subtraction of all science frames, and then it should also trim all science frames and remove cosmic rays from all of them. Outputs a 2-d spectrum that is in electron counts.

STEP 2 of Reduction
This is the first script in the real reduction of the step-by-step process.

Requires 'listZero' file  of all the desired Zero frames and 'listSpec' file of all the science frames to be used, meaning on-target or the lamp frames that bracket.

Creates 'listCTB' at the end for use in the next step of the process.
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

#import cal_params as cp
import get_cal_params as gcp

zerolistname= 'listZero'
speclistname = 'listSpec'
flatlistname= 'listFlat'


zerolist = np.genfromtxt(zerolistname, dtype ='str')
speclist= np.genfromtxt(speclistname, dtype = 'str')
flatlist= np.genfromtxt(flatlistname, dtype='str')

####
slit_ystart = 1   #The beginning of the image that has light from outside
slit_yend= 199     #The end of the image with same
#im_xstart = 9
#im_xend = 2055

#red cam

#im_params= {'Blue':{'im_xstart': 9, 'im_xend':  2055},'Red':{'im_xstart':26,'im_xend':2071}}

#technically there doesn't appear to be any end trimming neaded with red...




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


bias_stack,gain, readnoise, bias_times = make_image_stack(zerolist, times= False)
bias_med = np.nanmedian(bias_stack, axis=0)

filename= glob(speclist[1])[0]
print filename
#i = fits.open(speclist[1])
header= fits.getheader(speclist[1])
#im_xstart=im_params[header['INSTCONF']]['im_xstart']
#im_xend= im_params[header['INSTCONF']]['im_xend']



new_file_list= []
for img in speclist:
    filename= glob(img)[0]
    print filename
    i = fits.open(img)
    header= fits.getheader(img)
    setup_dict=gcp.get_cal_params(header)
    trim_regions=setup_dict['trimregions']
    gain =header['GAIN']
    readnoise = header['RDNOISE']
    img_data= i[0].data
    img_data= img_data-bias_med #bias subtraction
    img_data = img_data[trim_regions['y'][0]:trim_regions['y'][1],trim_regions['x'][0]:trim_regions['x'][1]] #trimming the edges
    target_cosmic= cosmics.cosmicsimage(img_data, gain=gain, readnoise=readnoise, sigclip = 5.0, sigfrac = 0.3, objlim = 5.0)
    target_cosmic.run(maxiter= 4)
    img_data= target_cosmic.cleanarray
    new_filename= 'ctb.' + filename
    img_data = img_data * gain
    header.append(card =( 'Counts', 'True', 'if spectrum in counts'))
    new_file_list.append(new_filename)
    new_hdu = fits.PrimaryHDU(img_data, header = header)
    #new_hdu.verify('fix')
    #new_hdu.writeto(new_filename, overwrite = True)
    new_hdu.writeto(new_filename, output_verify= 'fix', overwrite = True)
    
new_file_array = np.array(new_file_list)
print new_file_array

np.savetxt('listCTB', new_file_array, fmt = '%s')

########## Now doing the flats

new_file_list= []
for img in flatlist:
    filename= glob(img)[0]
    print filename
    i = fits.open(img)
    header= fits.getheader(img)
    setup_dict=gcp.get_cal_params(header)
    trim_regions=setup_dict['trimregions']
    gain =header['GAIN']
    readnoise = header['RDNOISE']
    img_data= i[0].data
    img_data= img_data-bias_med #bias subtraction
    #img_data = img_data[slit_ystart:slit_yend,im_xstart:im_xend] #trimming the edges
    img_data = img_data[trim_regions['y'][0]:trim_regions['y'][1],trim_regions['x'][0]:trim_regions['x'][1]] #trimming the edges
    target_cosmic= cosmics.cosmicsimage(img_data, gain=gain, readnoise=readnoise, sigclip = 5.0, sigfrac = 0.3, objlim = 5.0)
    target_cosmic.run(maxiter= 4)
    img_data= target_cosmic.cleanarray
    new_filename= 'ctb.' + filename
    img_data = img_data * gain
    header.append(card =( 'Counts', 'True', 'if spectrum in counts'))
    new_file_list.append(new_filename)
    new_hdu = fits.PrimaryHDU(img_data, header = header)
    #new_hdu.verify('fix')
    #new_hdu.writeto(new_filename, overwrite = True)
    new_hdu.writeto(new_filename, output_verify= 'fix', overwrite = True)
    
new_file_array = np.array(new_file_list)
print new_file_array

np.savetxt('listCTBflat', new_file_array, fmt = '%s')
    
    
    
    
