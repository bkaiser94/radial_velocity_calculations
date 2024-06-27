"""
Created by Ben Kaiser (UNC-Chapel Hill) 2024-04-28


This script should basically replicate the specfoc function in IRAF from IC7. I realize that PyRAF can almost certainly do that, but given that it requires some miniconda environment, I suspect that it will require reinstallation of conda, which screws up the real python path, so I'm going to try making a simple code to fit FWHM values to a CCD image.

By the way, I'm coding this to work with a slit image taken with bessel-R for the time-being, (so only a single peak to fit FWHM to) which is why I'm bothtering to try this manually.



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


input_files=sorted(glob('*besselR*focus*.fits'))
#input_files=sorted(glob('*besselr*focus*.fits')) #need this one for 2024-04-07 only


glob_search_width=30
#center_search_val=2067 #2024-04-17
#center_search_val=2135. #2024-02-07
#center_search_val=2130. #2024-02-10
#center_search_val=2134. #2024-03-04. 
#center_search_val=2130. #2024-03-12. 
#center_search_val=2130. #2024-04-07. 
#center_search_val=2126. #2024-05-07
#center_search_val=2111. #2024-05-08
#center_search_val=2116. #2024-05-14
center_search_val=2115. #2024-06-04




p0_slit_list=[1500,center_search_val,3,500]

refined_focus_width=900.
refined_focus_add=300.




################################

def gaussian_curve(x, a, x0, sigma,b):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+b

def fit_parabola(focus_list, fwhm_list, plot_all=False):
    focus_list=np.array(focus_list)
    fwhm_list=np.array(fwhm_list)
    def parabola(x, a, b, c):
        return a*x**2+b*x+c
    
    popt, pcov=sciop.curve_fit(parabola, focus_list,fwhm_list)
    print('popt', popt)
    min_focus=-1.*popt[1]/(2.*popt[0])
    print('min focus:', min_focus)
    print('fwhm for min_focus: ', parabola(min_focus,popt[0],popt[1],popt[2]))
    if plot_all:
        fine_xvals=np.linspace(np.min(focus_list),np.max(focus_list),100)
        fine_yvals=parabola(fine_xvals, popt[0],popt[1],popt[2])
        plt.plot(fine_xvals,fine_yvals,label='first parabolic fit')
        plt.axvline(x=min_focus, color='k', linestyle='--', label='first min focus: '+str(min_focus))
    else:
        pass
    
    #second_focus_inds=np.where(np.abs(focus_list-min_focus)<refined_focus_width)
    #second_focus_list=focus_list[second_focus_inds]
    #second_fwhm_list=fwhm_list[second_focus_inds]
    #popt2, pcov2=sciop.curve_fit(parabola,second_focus_list,second_fwhm_list)
    #second_min_focus=-1*popt2[1]/(2.*popt2[0])
    #second_fwhm=parabola(second_min_focus, popt2[0],popt2[1],popt2[2])
    def improve_solution(focus_list, min_focus):
        second_focus_inds=np.where(np.abs(focus_list-min_focus)<refined_focus_width)
        second_focus_list=focus_list[second_focus_inds]
        second_fwhm_list=fwhm_list[second_focus_inds]
        print('second_fwhm_list', second_fwhm_list)
        print('second_fwhm_list', second_fwhm_list)
        if len(second_focus_list < 3):
            print('too few indices included in search')
            print('tuple version', second_focus_inds)
            second_focus_inds=np.where(np.abs(focus_list-min_focus)<refined_focus_width+refined_focus_add)
            print('expanded inds', second_focus_inds)
            second_focus_list=focus_list[second_focus_inds]
            second_fwhm_list=fwhm_list[second_focus_inds]
            #second_focus_inds=list(second_focus_inds[0],second_focus_inds[1])
            #print('second',second_focus_inds)
            #second_focus_inds.append(2)
            #new_focus_inds=second_focus_inds
            #print('new',new_focus_inds)
            #second_focus_list=focus_list[new_focus_inds]
            #second_fwhm_list=focus_list[new_focus_inds]
        else:
            pass
        print('second_focus_list',second_focus_list)
        print('second_fwhm_list', second_fwhm_list)
        popt2, pcov2=sciop.curve_fit(parabola,second_focus_list,second_fwhm_list)
        second_min_focus=-1*popt2[1]/(2.*popt2[0])
        second_fwhm=parabola(second_min_focus, popt2[0],popt2[1],popt2[2])
        
        return second_min_focus, second_fwhm, popt2, second_focus_list
    second_min_focus, second_fwhm, popt2, second_focus_list=improve_solution(focus_list, min_focus)
    print('second min focus', second_min_focus)
    print('fwhm for second_min_focus',second_fwhm )
    print('fwhm in arcseconds', second_fwhm*0.15)
    if plot_all:
        fine_xvals=np.linspace(np.min(second_focus_list),np.max(second_focus_list),100)
        fine_yvals=parabola(fine_xvals, popt2[0],popt2[1],popt2[2])
        plt.plot(fine_xvals,fine_yvals,label='second parabolic fit')
        plt.axvline(x=second_min_focus, color='r', linestyle=':', label='second min focus: '+str(second_min_focus))
    else:
        pass
    
    third_min_focus, third_fwhm, popt3, third_focus_list=improve_solution(focus_list, second_min_focus)
    print('third min focus:', third_min_focus)
    print('fwhm third_min_focus', third_fwhm)
    print('fwhm in arcseconds', third_fwhm*0.15)
    if plot_all:
        fine_xvals=np.linspace(np.min(third_focus_list),np.max(third_focus_list),100)
        fine_yvals=parabola(fine_xvals, popt3[0],popt3[1],popt3[2])
        plt.plot(fine_xvals,fine_yvals,label='third parabolic fit')
        plt.axvline(x=third_min_focus, color='b', linestyle=':', label='third min focus: '+str(third_min_focus))
    else:
        pass
    return


def fit_gaussian_curve(x_pixels, light_values, p0_list, search_width=glob_search_width, plot_all = False, bounds = (-np.inf, np.inf), fixed_width=True):
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

        #print('\n\n\n=======\n')
        #print(filename)
        #print('bkg_core_sides right before cross section plot', bkg_core_sides)
        #print('\n=========\n\n\n')


        plt.legend()
        plt.show()
    else:
        pass
    return popt, pcov


#######################################


fwhm_list=[]
focus_list=[]





for input_file in input_files:
    hdu= fits.open(input_file)
    header=fits.getheader(input_file)
    focus_list.append(header['cam_foc'])
    image=hdu[0].data
    #collapsed_spec=np.mean(image,axis=0)
    collapsed_spec=np.mean(image[90:110],axis=0)
    print('collapsed_spec.shape',collapsed_spec.shape)
    x_pixels=np.indices(collapsed_spec.shape)[0]
    print(x_pixels[0])
    popt, pcov=fit_gaussian_curve(x_pixels, collapsed_spec, p0_slit_list, plot_all=True)
    sigma=popt[2]
    fwhm=2*np.sqrt(2*np.log(2))*sigma
    fwhm_list.append(fwhm)
    
fit_parabola(focus_list, fwhm_list,plot_all=True)
plt.plot(focus_list,fwhm_list, marker='o')
plt.xlabel('Focus')
plt.ylabel('FWHM unbinned pixels')
plt.legend()
plt.show()

plt.plot(focus_list,np.array(fwhm_list)*0.15)
plt.xlabel('Focus')
plt.ylabel('FWHM arcseconds')
plt.show()




dates=[
    '2024-02-07',
    '2024-02-10',
    '2024-03-04',
    '2024-03-12',
    '2024-04-07',
    '2024-04-17',
    '2024-05-07',
    '2024-05-08',
    '2024-05-14',
    '2024-06-04'
    ]


focii=[
    0.458,
    0.44,
    0.424,
    0.401,
    0.5091,
    0.633,
    0.583,
    0.627,
    0.398,
    0.401
    ]


col_swaps=[
    '2024-03-06',
    '2024-05-10'
    ]

#old method using astropy dates/time 
#dates=Time(dates)
#plt.plot(dates.mjd, focii,marker='o')

#switched to numpy datetime per matplotlib documentation
#dates=np.datetime64(np.array(dates))
dates=np.array(dates, dtype='datetime64')
plt.plot(dates, focii,marker='o')
#plt.axvline(x=Time('2024-03-06').mjd, color='k', linestyle='--', label='Collimator swap (2024-03-06)')
for coldate in col_swaps:
    #plt.axvline(x=Time(coldate).mjd, color='k', linestyle='--', label='Collimator swap '+coldate)
    plt.axvline(x=np.datetime64(coldate), color='k', linestyle='--', label='Collimator swap '+coldate)

plt.axhline(y=0.45, label='0.45" slit width', color='r', linestyle=':')
plt.legend()
plt.ylabel('FWHM in arcseconds of best focus from focus tests')
#plt.xlabel('Date (MJD)')
plt.xlabel('Date')
plt.xticks(rotation=70)
plt.grid()
plt.title('bessel-R image of 0.45" slit FWHM (in arcseconds) implied by focus minimization')
plt.show()






