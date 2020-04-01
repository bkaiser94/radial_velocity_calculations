"""
Created by Ben Kaiser (UNC-Chapel Hill) 2020-01-06

This is essentially a Python implementation of some of the routines/equations presented in Hurley, Pols, and Tout 2000 with an emphasis (and probably entirety of implementation) focused on being used for white dwarf progenitor MS lifetimes.


"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy import constants as const
from astropy.time import Time
from astropy.table import Table
import scipy.interpolate as scinterp



import cal_params as cp
from a_coeffs import a_coeffs

#nothing

default_z=0.0001 #allegedly thick disk value
#default_z=0.02 #approximately solar


def make_match(value, array):
    try:
        output_val=np.ones(array.shape)*value
    except AttributeError:
        output_val=value
    return output_val

def zeta(z):
    return np.log10(z/0.02)

def a(z,index):
    """
    the a_n equation from the beginning of Appendix A of Hurley et al. 2000
    
    I just left out the mu term because none of these variables have a mu value, and I didn't bother putting it in the
    a_coeffs list either as a result.
    """
    zeta_val=zeta(z)
    a_index=a_coeffs[index]
    a_out=a_index[0]+a_index[1]*zeta_val+a_index[2]*zeta_val**2+a_index[3]*zeta_val**3
    #print('a_'+str(index),a_out)
    return a_out

def get_t_bgb(z,mass):
    #print('a_1',a(z,1))
    #print('a(z,2)*mass**4',a(z,2)*mass**4)
    #print('a(z,3)*mass**5.5 ',a(z,3)*mass**5.5 )
    #print('mass**7',mass**7)
    #print('(a(z,1)+a(z,2)*mass**4+a(z,3)*mass**5.5 + mass**7)',(a(z,1)+a(z,2)*mass**4+a(z,3)*mass**5.5 + mass**7))
    #print('a(z,4)*mass**2',a(z,4)*mass**2)
    #print('a(z,5)*mass**7',a(z,5)*mass**7)
    #print('(a(z,4)*mass**2 + a(z,5)*mass**7)',(a(z,4)*mass**2 + a(z,5)*mass**7))
    return(a(z,1)+a(z,2)*mass**4+a(z,3)*mass**5.5 + mass**7)/(a(z,4)*mass**2 + a(z,5)*mass**7)

def get_x(z):
    first_val=make_match(0.95,z)
    last_val=make_match(0.99, z)
    return np.max([first_val,
                   np.min([
                       0.95-0.03*(zeta(z)+0.30103),
                       last_val
                       ],axis=0)],axis=0)

def get_mu(z,mass):
    #print(np.max([a(z,6)/(mass**a(z,7)),a(z,8)+(a(z,9)/mass**(a(z,10)))],axis=0))
    first_val=make_match(0.5, mass)
    return np.max([
        first_val,
        1.0-0.01 * np.max([
            a(z,6)/(mass**a(z,7)),
            a(z,8)+(a(z,9)/(mass**a(z,10)))
            ],axis=0)],axis=0)


def get_t_ms(mass, z=default_z):
    """
    Default is what I'm pretty sure is solar metallicity
    """
    t_bgb=get_t_bgb(z,mass)
    return np.max([
        get_mu(z,mass)*t_bgb,
        get_x(z)*t_bgb
        ],axis=0)*1e-3



def get_B(mass):
    """
    un-numbered equation below equation(38)
    """
    return np.max([3e4, 500+1.75e4 * mass**0.6], axis=0)


def get_D(mass, zeta, mass_HeF=2.5):
    """
    
    un-numbered equation above equation (39) 
    
    I'm pretty sure I need to add in a linear interpolation element to cover the gap (if there is one) between mass_HeF and 2.5
    """
    D_lo=5.37+0.135*zeta
    D_hi=np.max([-1.0,0.975*D_lo-0.18*mass, 0.5*D_lo-0.06*mass],axis=0)
    lo_inds=np.where(mass <= mass_HeF)
    hi_inds= np.where(mass >= 2.5)
    D_hi[lo_inds]=D_lo
    return 10**D_hi

def get_core_luminosity_rel(mass, mass_core, mass_HeF=2.5):
    B=get_B(mass)
    p_array=np.ones(mass.shape)
    q_array=np.ones(mass.shape)
    low_inds= np.where(mass <= mass_HeF)
    high_inds= np.where(mass > mass_HeF)
    p_array[low_inds]=6
    p_array[high_inds]=5
    q_array[low_inds]=3
    q_array[high_inds]=2
    
    
    return


def get_t_hems(mass):
    """
    Eq. (79) from Hurley et al. 2000
    """
    return (0.4129 + 18.81 * mass**4+1.853*mass**6)/mass**6.5

def get_t_he():
    
    
    return
#test_mass=1.
#print('zeta', zeta(default_z))
#print('t_bgb', get_t_bgb(default_z, test_mass))
#print('mu', get_mu(default_z, test_mass))
#print('x', get_x(default_z))
#print('t_ms',get_t_ms(test_mass))


#mass_vals=np.linspace(0.1,13.,100)
#log_mass_vals=np.linspace(-0.4, 1.9, 100)
#mass_vals= 10.**log_mass_vals

#t_ms_vals= get_t_ms(mass_vals)
#t_bgb_vals=get_t_bgb(default_z, mass_vals)


#plt.plot(np.log10(mass_vals), np.log10(t_ms_vals))
##plt.plot(np.log10(mass_vals), np.log10(t_bgb_vals))
#plt.plot(np.log10(mass_vals), np.log10(get_t_ms(mass_vals,z=0.0001)))
#plt.xlim(-0.4, 1.9)
#plt.show()


