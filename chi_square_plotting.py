"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-07-16

It should read in the output from the chi-squared measurements for the model atmosphere fitting and then make 
whatever plots I want without having to rerun the entire fitting code everytime. This should have been done a 
long time ago.

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
from astropy import convolution as conv
import scipy.interpolate as scinterp



import wdatmos
import spec_plot_tools as spt
import kernel_builder


#input_filename= 'chi_square_values.csv'
#input_filename= 'chi_square_values_noca.csv'
#input_filename= 'chi_square_values_quad.csv'
#input_filename= 'chi_square_values_balm.csv'
#input_filename= 'chi_square_values_noerr.csv'
#input_filename= 'chi_square_values_norm_.csv'
input_filename= 'chi_square_values_norm_metalmask.csv'



#clevels= np.arange(0,1.5, 0.01)
clevels= np.arange(0,1.5, 0.05)
clevels= np.arange(0,100, 1)

dclevels= np.array([1,2.3, 3.53,4.72,5.89,7.04])

input_array = np.genfromtxt(input_filename, names= True, delimiter= ',')
teff_array = input_array['teff']
logg_array = input_array['logg']
rv_array = input_array['rv']
dist_array = input_array['chi_square']

def rescale_chi_square(chi_square_vals):
    chi_square_min = np.min(chi_square_vals)
    chi_square_vals= chi_square_vals/chi_square_min


def chi_square_dots(teff_array, logg_array, dist_array):
    min_index = np.argmin(dist_array)
    print "Teff and logg min chi-squared values: ", teff_array[min_index],logg_array[min_index], "|chi-sq:", dist_array[min_index]
    contour_array = np.vstack([teff_array,logg_array, dist_array])
    #plt.imshow(contour_array, aspect= 100)
    #plt.contour(teff_array, logg_array, dist_array)
    marker_scale = 1/dist_array* dist_array.min() *40.
    #plt.scatter(teff_array, logg_array, s= 1./dist_array*30, c = 1./dist_array*20)
    plt.scatter(teff_array, logg_array, s=marker_scale, c = marker_scale)
    plt.plot(teff_array[min_index],logg_array[min_index], marker = '*', markersize = 14)
    plt.xlabel('T_eff')
    plt.ylabel('logg')
    #plt.show()
    
def chi_square_contours(teff_array, logg_array, dist_array, clevels= clevels, get_bounds= False):
    min_index = np.argmin(dist_array)
    print "Teff and logg min chi-squared values: ", teff_array[min_index],logg_array[min_index], "|chi-sq:", dist_array[min_index]
    #contour_array = np.vstack([teff_array,logg_array, dist_array])
    #plt.imshow(contour_array, aspect= 100)
    #plt.contour(teff_array, logg_array, dist_array)
    #marker_scale = 1/dist_array* dist_array.min() *40.
    #plt.scatter(teff_array, logg_array, s= 1./dist_array*30, c = 1./dist_array*20)
    #plt.scatter(teff_array, logg_array, s=marker_scale, c = marker_scale)
    logg_array_d2 = np.copy(logg_array.reshape(11,36))
    teff_array_d2= np.copy(teff_array.reshape(11,36))
    dist_array_d2= np.copy(dist_array.reshape(11,36))
    #plt.contour(teff_array, logg_array, dist_array)
    contours= plt.contour(teff_array_d2, logg_array_d2, dist_array_d2, levels = clevels, colors= 'black')
    
        #print thing.get_segments()
       
    plt.clabel(contours, inline=True)
    plt.plot(teff_array[min_index],logg_array[min_index], marker = '*', markersize = 14)
    plt.xlabel('T_eff')
    plt.ylabel('logg')
    plt.show()
    if get_bounds:
        print "\n************\nOnly works for a single level input\n*************\n"
        for thing in contours.collections:
            segments = thing.get_segments() #tuple of the line segments that comprise each layer
            try:
                total_array=np.vstack([segments[0],segments[1],segments[2]])
            except IndexError:
                try:
                    total_array=np.vstack([segments[0],segments[1]])
                except IndexError:
                    total_array= segments
            uppers=np.max(total_array, axis=0)
            lowers= np.min(total_array, axis=0)
            return uppers, lowers
            #for seg in segments:
                #print seg.shape
                #print seg
                #uppers= np.max(seg, axis=0)
                #lowers = np.min(seg, axis=0)
                #return uppers, lowers

print logg_array.shape[0] / 11.
#logg_array_d2 = np.copy(logg_array.reshape(11,36))
#teff_array_d2= np.copy(teff_array.reshape(11,36))
#dist_array_d2= np.copy(dist_array.reshape(11,36))
##plt.contour(teff_array, logg_array, dist_array)
#contours= plt.contour(teff_array_d2, logg_array_d2, dist_array_d2, levels = clevels, colors= 'black')
#plt.clabel(contours, inline=True)
#plt.xlabel('T_eff')
#plt.ylabel('logg')
#plt.show()
#chi_square_contours(teff_array, logg_array, dist_array, get_bounds= True)
chi_square_dots(teff_array,logg_array, dist_array)
#chi_square_contours(teff_array, logg_array, dist_array-np.min(dist_array), clevels= dclevels)
upper_bound, lower_bound= chi_square_contours(teff_array, logg_array, dist_array-np.min(dist_array), clevels= [1], get_bounds= True)
print "Upper bounds: ", upper_bound
print "Lower bounds: ", lower_bound


marker_scale = (logg_array-logg_array.min())*10
plt.scatter(teff_array+(np.random.rand(teff_array.shape[0])*50), dist_array, s=marker_scale)
plt.xlabel("Teff (K)")
plt.ylabel(r"red $\chi^2$")
plt.show()

teff_scale = -1*(teff_array-np.min(teff_array))
plt.scatter(logg_array+(np.random.rand(logg_array.shape[0])*0.1), dist_array, c = teff_scale)
#plt.scatter(logg_array+np.linspace(0,0.1, 20), dist_array, c = teff_scale)
plt.xlabel("log(g)")
plt.ylabel(r"red $\chi^2$")
plt.show()
