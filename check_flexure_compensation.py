"""
Created by Ben Kaiser (UNC-Chapel Hill) on 2022-06-01

This script is intended to check if flexure compensation on the Goodman Spectrograph is 
working, or to at least make plots that let a user check if flexure compensation was working on 
the spectrograph.


"""

from __future__ import print_function

import numpy as np
import sys
import glob
#from glob import glob
import matplotlib.pyplot as plt
from astropy.io import fits
#from astropy import units as u
#from astropy import constants as const
from astropy.time import Time
from astropy.table import Table
#import scipy.interpolate as scinterp
from astropy import coordinates as coords

import cal_params as cp

import plot_spec as ps

#print(glob.__version__)

filename_base="*.fits"
#listfile='listTest'
#filelist=np.genfromtxt(listfile,dtype='str')



#plot_title="Sheila's Group observing 2022-05-25"
plot_title='Ben observing 2022-05-21'

filelist=glob.glob(filename_base)
filelist=sorted(filelist)
#print(filelist)
def plot_delta_cam(file_list, x_val='default'):
    delta_cam_list=[]
    rot_angle_list=[]
    counter=0
    for filename in file_list:
        header=fits.getheader(filename)
        actual=header['cam_ang']
        target=header['cam_targ']
        delta_cam=actual-target
        delta_cam_list.append(delta_cam)
        rot_angle=header['rotator']
        #if (rot_angle>180.):
            #rot_angle=rot_angle-360.
        #else:
            #pass
        if x_val=='default':
            rot_angle_list.append(rot_angle)
        else:
            rot_angle_list.append(counter)
        #plt.text(rot_angle,delta_cam, filename)
        counter+=1
    
    plt.scatter(rot_angle_list,delta_cam_list,color='b')        
    plt.ylim(-0.04,0.04)
    plt.grid()
    if x_val=='default':
        #plt.xlim(-180,180)
        plt.xlim(0,360)
        plt.xlabel('rotator')
        xvals=np.linspace(0,360,10000)
        yvals=-0.033*np.sin(xvals/180*np.pi)
        
        plt.plot(xvals,yvals)
    else:
        plt.xlabel('file order')
    plt.ylabel(r'$\Delta$Camera Angle (cam_ang-cam_targ)')
    
    plt.show()
    return


def plot_head_v_order(file_list,header_name):
    counter=0
    header_list=[]
    counter_list=[]
    for filename in file_list:
        header=fits.getheader(filename)
        header_list.append(header[header_name])
        counter_list.append(counter)
        counter+=1
        
        
    plt.scatter(counter_list, header_list)
    plt.ylabel(header_name)
    plt.xlabel('file order')
    plt.show()
    
    return



header_list=[
    ['rotator','cam_ang']
    
    
    
    
    ]


#plt.title('Ben observing 2022-05-21')
print_rot_vals=np.arange(0,375,15.)
print_ang_vals=-0.033*np.sin(print_rot_vals/180.*np.pi)
print('Rotator', 'Actual-Target')
for rot_val, ang_val in zip(print_rot_vals,print_ang_vals):
    print(rot_val, np.round(ang_val,2))

plt.title(plot_title)
plot_delta_cam(filelist)

plt.title(plot_title)
plot_delta_cam(filelist,x_val='something')

plt.xlim(-180,180)
#plt.title('Ben observing 2022-05-21')
plt.title(plot_title)
ps.plot_head_2_head(filelist,'rotator','cam_targ')

plt.title(plot_title)
plot_head_v_order(filelist, 'cam_targ')

plt.title(plot_title)

plot_head_v_order(filelist,'exptime')









