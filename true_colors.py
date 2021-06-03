"""
Created by Ben Kaiser (UNC-Chapel Hill) 2021-06-03


Take an optical spectrum and generate an RGB circle (or file) that shows the "true" color of the
spectrum. This could probably be generalized to spectra that do not cover the full visible range 
and the code infers the slope at those other parts, but initially it's just going to to work for my
Goodman data of J1644, and it will probably actually use the stitched spectrum that I used to 
make the figure 1 of the Science paper.

The spectrum will be rescaled to the CIE 1931 XYZ system
(https://en.wikipedia.org/wiki/CIE_1931_color_space ;wikipedia link because it's plain-ish language, but I didn't use the math here), using simple
analytic color-matching functions of equation 2 from Wyman et al. 2013
(http://jcgt.org/published/0002/02/01/).

I'm then using skimage.color.xyz2rgb() to convert these XYZ colors to the sRGB color space, 
which I'm pretty sure is the endzone for the color transformation.

Then I just have to make a circle with that color. I'm 80% sure this last step will inexplicably be 
the most difficult for me based on my extreme difficulties with simple tasks *shrug*.


"""



from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
import skimage.color


import spec_plot_tools as spt
import cal_params as cp


def x1931(wave):
    """
    from equation 2 of Wyman et al. 2013
    
    wave should be in nanometers, I'm 99% sure; I didn't see Wyman et al. explicitly state to do 
    so, but they discuss wavelengths in nm throughout the paper, so it's not a crazy leap.
    
    """
    
    return 1.065*np.exp(-0.5*((wave-595.8)/33.33)**2)+0.366 * np.exp(-0.5 * ((wave-446.8)/19.44)**2)

def y1931(wave):
    """
    from equation 2 of Wyman et al. 2013
    
    wave should be in nanometers, I'm 99% sure; I didn't see Wyman et al. explicitly state to do 
    so, but they discuss wavelengths in nm throughout the paper, so it's not a crazy leap.
    
    """

    return  1.014 * np.exp(-0.5*((np.log(wave)-np.log(556.3))/0.075)**2)




def z1931(wave):
    """
    from equation 2 of Wyman et al. 2013
    
    wave should be in nanometers, I'm 99% sure; I didn't see Wyman et al. explicitly state to do 
    so, but they discuss wavelengths in nm throughout the paper, so it's not a crazy leap.
    
    """
    
    return 1.839 * np.exp(-0.5*((np.log(wave)-np.log(449.8))/0.051)**2)






