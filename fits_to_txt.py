"""
Created by Ben Kaiser 2019-05-14 (UNC-Chapel Hill)

Read in spectra in the FITS format typical of my own 'radial_velocity_calcuations/' spectra and then output
them as .txt files that can be used by PyHammer...hopefully.

The output filenames should be identical in everyway except that they'll be .txt files instead
ACtually they're going to be CSV files I think. I should probably make it capable of outputting txt 
files too, but the priority will be CSV files


"""
from __future__ import print_function
import numpy as np
#import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
#from astropy.time import Time
#from astropy import coordinates as coords
#from astropy import units as u
#from astropy import constants as const
#from astropy import convolution as conv
#import scipy.interpolate as scinterp
#import time
#start = time.time()



#print start
import spec_plot_tools as spt


#input_string= 'avg_fwctb*fits'
#input_string=''
input_filename=sys.argv[1]
print('input_filename:',input_filename)

hdu=fits.open(input_filename)
wavelengths=hdu[0].data
flux=hdu[1].data
dlambda=hdu[4].data

output_array=np.vstack([wavelengths,flux,dlambda]).T

output_name_base=input_filename[:-4]
print('output_name_base',output_name_base)
output_filename=output_name_base+'csv'
print("saving", output_filename)
np.savetxt(output_filename,output_array, delimiter=',')
print("saved", output_filename)

#input_filenames= glob(input_string)


#new_array=np.genfromtxt(output_filename,delimiter=',')
#print(output_array)
#print(new_array)
#print(output_array-new_array)

#print(new_array[0])
