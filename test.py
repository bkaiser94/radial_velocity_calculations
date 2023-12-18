"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-08-05


Random script for testing random things without having to get a whole thing working in the actual reduction 
process.



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


#input_file= 'ctb.0312_EG274_400m2.fits'

#i= fits.open(input_file)
#header = fits.getheader(input_file)
#img_data= np.copy(i[0].data)

#spt.rebin_image(img_data, rebin_axis=1, rebin_num=10, plot_all=True)

glob_string='*.fits'

filenames=glob(glob_string)
filenames=sorted(filenames)

filenum_list=[]
header_present_list=[]
exp_time_list=[]
marker_list=[]
for i,filename in enumerate(filenames):
    print(filename)
    header =fits.getheader(filename)
    filenum_list.append(i)
    try:
        exp_time_list.append(header['exptime'])
        timeval=header['closetim']
        header_present_list.append(1)
        marker_list.append('g')
    except KeyError as error:
        print('KeyError:',error)
        header_present_list.append(0)
        marker_list.append('r')
        
        
plt.plot(filenum_list,header_present_list,marker='o')
plt.title('whether or not there is an opentime header')
plt.show()



#plt.plot(filenum_list, exp_time_list, marker='o')
plt.scatter(filenum_list, exp_time_list, c=marker_list)
plt.title("Exp times")
plt.yscale('log')
plt.show()
