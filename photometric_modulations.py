"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-08-14

Implementation of the ellipsoidal modulation equation from Morris and Naftilan 1993... or at least 
an attempt to do so.

Ideally, this will expand to include the other types of photometric modulation in the future.



"""
import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
from astropy.units import cds
import scipy.stats as scistats
import scipy.optimize as sciop



## Actually all of the phi values that are input here need to be converted to be in radians instead of 
## fractions of phase since I'm giving them to python functions
def constant_terms(u1, tau1, R1, A, q, i):
    first_half= 1+(15+u1)*(1+tau1)*(R1/A)**3*(2+5*q)*(2-3*np.sin(i)**2)/(60*(3-u1))
    second_half= 9*(1-u1)*(3+tau1)*(R1/A)**5*q*(8-40*np.sin(i)**2 + 35*np.sin(i)**4)/(256*(3-u1))
    return first_half + second_half

def first_cos(u1, tau1, R1, A, q, i):
    first_half= 15*u1*(2+tau1)*(R1/A)**4*q*(4*np.sin(i)-5*np.sin(i)**3)/(32*(3-u1))
    return first_half

def second_cos(u1, tau1, R1, A, q, i):
    first_half= -3*(15+u1)*(1+tau1)*(R1/A)**3*q*np.sin(i)**2/(20*(3-u1))
    second_half= -15*(1-u1)*(3+tau1)*(R1/A)**5 * q* (6*np.sin(i)**2 -7 *np.sin(i)**4)/(64*(3-u1))
    return first_half+second_half

def third_cos(u1, tau1, R1, A, q, i):
    first_half= -25*u1*(2+tau1)*(R1/A)**4*q*np.sin(i)**3/(32*(3-u1))
    return first_half

def fourth_cos(u1, tau1, R1, A, q, i):
    first_half= 105*(1-u1)*(3+tau1)*(R1/A)**5*q*np.sin(i)**4/(256*(3-u1))
    return first_half


def ellipsoidal_brightness(L0, phi, u1, tau1, R1, A, q, i):
    """
    return the brightness of the star at a given phase. The phase is in terms of phi=0 when the
    distorted star is farthest from the observer, which means it is 0.5 offset from the phi I use in the
    other parts of this repository because that's how astronomy is.
    """
    return L0*(constant_terms(u1, tau1, R1, A, q, i)+first_cos(u1, tau1, R1, A, q, i)*np.cos(phi)+second_cos(u1, tau1, R1, A, q, i)*np.cos(2*phi)+ third_cos(u1, tau1, R1, A, q, i)*np.cos(3*phi)+fourth_cos(u1, tau1, R1, A, q, i)*np.cos(4*phi))


phi_range = np.linspace(0,2*np.pi,1000)

brightnesses= ellipsoidal_brightness(1,phi_range, 0.1, 0.1, 0.1, 100, 0.12, np.pi/3)

plt.plot(phi_range, brightnesses)
plt.show()
               
               
               
               
            
