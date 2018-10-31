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
#input_filename= 'chi_square_values_norm_metalmask.csv'
#input_filename=  'chi_square_values_norm_wholespec.csv'
#input_filename= 'chi_square_values_davchange.csv'
#input_filename='chi_square_values_20180910.csv'
input_filename= 'chi_square_values_20180911small.csv'
#input_filename= 'chi_square_values_20180911.csv'
#input_filename='20181008_run1.csv'
#input_filename= '20181010_run1.csv'
input_filename= '20181021B_31spec.csv'

#clevels= np.arange(0,1.5, 0.01)
clevels= np.arange(0,1.5, 0.05)
clevels= np.arange(0,100, 1)
delta_chi2_bound = 2.3
parabola_fit_width = 3

dclevels= np.array([1,2.3, 3.53,4.72,5.89,7.04])

input_array = np.genfromtxt(input_filename, names= True, delimiter= ',')
teff_array = input_array['teff']
logg_array = input_array['logg']
rv_array = input_array['rv']
dist_array = input_array['chi_square']
rescale_dist= input_array['revised_chi_square']

def rescale_chi_square(chi_square_vals):
    chi_square_min = np.min(chi_square_vals)
    chi_square_vals= chi_square_vals/chi_square_min

def extract_teffs_loggs(teff_array, logg_array, input_array, axis=0):
    """
    NOT DONE
    
    Clunky implementation of a version of np.argmin that can act along a single axis of an array
    """
    logg_array_d2 = np.copy(logg_array.reshape(11,36))
    teff_array_d2= np.copy(teff_array.reshape(11,36))
    dist_array_d2= np.copy(dist_array.reshape(11,36))
    min_chisq= np.nanmin(dist_array_d2)
    min_loc= np.where(dist_array_d2= min_chisq)
    return


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
    
def chi_square_parabolas(teff_array, logg_array, dist_array):
    """
    dist_array should be a delta chi-square array.
    """
    def parabola_func(x,a,b):
        return a*x**2+b
    min_index= np.argmin(dist_array)
    min_teff= teff_array[min_index]
    min_logg= logg_array[min_index]
    min_chisq= dist_array[min_index]
    print "Teff and logg min chi-squared values: ", teff_array[min_index],logg_array[min_index], "|chi-sq:", dist_array[min_index]
    #logg_array_d2 = np.copy(logg_array.reshape(11,36))
    #teff_array_d2= np.copy(teff_array.reshape(11,36))
    #dist_array_d2= np.copy(dist_array.reshape(11,36))
    #min_teff_dists= np.nanmin(dist_array_d2, axis=0)
    teff_loggs= np.copy(logg_array[np.where(teff_array== teff_array[min_index])])
    teff_chisq= np.copy(dist_array[np.where(teff_array== teff_array[min_index])])
    logg_teffs= np.copy(teff_array[np.where(logg_array== logg_array[min_index])])
    logg_chisq= np.copy(dist_array[np.where(logg_array== logg_array[min_index])])
    min_teff_loggs_i = np.argmin(teff_chisq)
    min_logg_teffs_i = np.argmin(logg_chisq)
    teff_loggs= teff_loggs[min_teff_loggs_i-parabola_fit_width: min_teff_loggs_i+parabola_fit_width+1]
    logg_teffs= logg_teffs[min_logg_teffs_i-parabola_fit_width: min_logg_teffs_i+parabola_fit_width+1]
    teff_chisq= teff_chisq[min_teff_loggs_i-parabola_fit_width: min_teff_loggs_i+parabola_fit_width+1]
    logg_chisq= logg_chisq[min_logg_teffs_i-parabola_fit_width: min_logg_teffs_i+parabola_fit_width+1]
    print min_index
    print min_teff_loggs_i
    print teff_loggs.shape
    print logg_teffs.shape
    print teff_chisq.shape
    print logg_chisq.shape
    popt_teff, cov= sciop.curve_fit(parabola_func, teff_loggs-min_logg, teff_chisq, bounds= ([0., -np.inf],[np.inf, delta_chi2_bound]))
    popt_logg, cov= sciop.curve_fit(parabola_func, logg_teffs-min_teff, logg_chisq, bounds= ([0., -np.inf],[np.inf, delta_chi2_bound]))
    delta_teff = np.sqrt((delta_chi2_bound-popt_logg[1])/popt_logg[0])
    delta_logg= np.sqrt((delta_chi2_bound-popt_teff[1])/popt_teff[0])
    print "popt_teff:", popt_teff
    print "popt_logg:", popt_logg
    print "Best fit: Teff=", min_teff, "+/-", delta_teff, "log(g)=", min_logg, "+/-", delta_logg, "| rescale_chi_square=", dist_array[min_index]
    plt.scatter(teff_loggs, teff_chisq)
    plt.plot(teff_loggs, parabola_func(teff_loggs-min_logg, popt_teff[0], popt_teff[1]))
    plt.xlabel('Log(g)')
    plt.ylabel('Chi-square')
    plt.show()
    plt.scatter(logg_teffs, logg_chisq)
    plt.plot(logg_teffs, parabola_func(logg_teffs-min_teff, popt_logg[0], popt_logg[1]))
    plt.xlabel('Teff')
    plt.ylabel('Chi-square')
    plt.show()
    return
    
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
    #print "dist_array_d2:",dist_array_d2
    #print "logg_array_d2:", logg_array_d2
    #print "teff_array_d2:", teff_array_d2
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

#print logg_array.shape[0] / 11.
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
chi_square_dots(teff_array,logg_array, rescale_dist)
#chi_square_contours(teff_array, logg_array, dist_array-np.min(dist_array), clevels= dclevels)
#upper_bound, lower_bound= chi_square_contours(teff_array, logg_array, dist_array-np.min(dist_array), clevels= [1], get_bounds= True)
#upper_bound, lower_bound= chi_square_contours(teff_array, logg_array, dist_array-np.min(dist_array), clevels= [2.3], get_bounds= True)
upper_bound, lower_bound= chi_square_contours(teff_array, logg_array, rescale_dist-np.min(rescale_dist), clevels= [2.3], get_bounds= True)
print "Upper bounds: ", upper_bound
print "Lower bounds: ", lower_bound

chi_square_parabolas(teff_array, logg_array,  rescale_dist-np.nanmin(rescale_dist))
marker_scale = (logg_array-logg_array.min())*10
#plt.scatter(teff_array+(np.random.rand(teff_array.shape[0])*50), dist_array, s=marker_scale)
plt.scatter(teff_array+(np.random.rand(teff_array.shape[0])*50), rescale_dist, s=marker_scale)
plt.xlabel("Teff (K)")
plt.ylabel(r"red $\chi^2$")
plt.show()

teff_scale = -1*(teff_array-np.min(teff_array))
plt.scatter(logg_array+(np.random.rand(logg_array.shape[0])*0.1), dist_array, c = teff_scale)
#plt.scatter(logg_array+np.linspace(0,0.1, 20), dist_array, c = teff_scale)
plt.xlabel("log(g)")
plt.ylabel(r"red $\chi^2$")
plt.show()
