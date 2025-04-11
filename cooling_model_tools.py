"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-04-10


Want this to basically create the cooling model tables and be imported by other scripts to 
be used. This should also contain the cooling model functions probably. It definitely should 
contain the functions to get back to mass from Teff and logg as that is the primary reason 
I'm making this file.



"""


from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from astropy.table import Table, Column
import scipy.interpolate as scinterp

import time


import spec_plot_tools as spt
import cal_params as cp

default_atm_type='thinH'
default_atm_type='thickH'

default_cooling_modeler='bedard2020'

default_interp_kind='quintic'

default_n=1e6

def get_cooling_model_table(cooling_modeler=default_cooling_modeler, atm_type=default_atm_type):
    if cooling_modeler== 'bedard2020':
        if atm_type=='thinH':
            cooling_model_file='bedard2020_seq_thinH.csv'
            interp_kind='cubic'
            print("\n\n*****\n\nbedard2020 with thinH atmosphere selected.\nSetting interp_kind='cubic'\nBecause I know it works correctly-ish.")
        elif atm_type=='thickH':
            cooling_model_file='bedard2020_seq_thickH.csv'
            interp_kind=default_interp_kind
    elif cooling_modeler=='fontaine2001':
        if atm_type=='thinH':
            cooling_model_file='COModel_ThinH.csv'
            print("\n\n*****\n\nfontaine2001 with thinH atmosphere selected.\nSetting default_interp_kind='cubic'\nBecause I know it works correctly-ish.")
            interp_kind='cubic'

        elif atm_type=='thickH':
            cooling_model_file='COModel_ThickH.csv'
            print("\n\n fontaine2001 with thickH atmosphere doesn't actually work correctly for some reason with any interpolation, so... we're gonna stop here. Pick a different setting. We also have the bedard2020 cooling_modeler option, which I have hopefully gotten running by the time you see this message...\n\n")
            sys.exit()
    else:
        print('\n\n**********\n\nInvalid cooling_modeler and atm_type selected:')
        print('cooling_modeler:',cooling_modeler,'atm_type:',atm_type,'\n********\n\n')

    cooling_model_file=cp.ref_dir+'WD_cooling_models/'+cooling_model_file
    cooling_table= Table.read(cooling_model_file)

    
    
    return cooling_table, interp_kind

def loggteff_to_m(teff, logg, cooling_table):
    loggteff_to_m_interp = scinterp.SmoothBivariateSpline(cooling_table['Teff'], cooling_table['logg'], cooling_table['Mass'])
    #print('doing loggteff_to_m')
    #print(loggteff_to_m_interp(teff,logg),loggteff_to_m_interp(teff,logg)[0])
    return loggteff_to_m_interp(teff,logg)[0]

def operate_on_dist(dist1, dist2, cooling_table,function):
    """
    dist1 is the distribution of the first input to 'function'
    
    dist2 is the distribution of the second input to 'function'
    
    function is the most likely scinterp.interp2d() function output or object or whatever that you've defined that you'd 
    like to feed distributions through and be able to track which inputs yielded them.
    
    For whatever reason, scinterp.interp2d() -generated objects seem to completely randomize the indices of the 
    outputs relative to the inputs. I have no idea why. It also wants to take the 2 input (10,) shape arrays and make 
    a (10,10) output, whose indices have no correlation to the indices of either input (10,) array. I know, that 
    seems ridiculous, but it's what I've been experiencing, so I just made a damn for-loop.
    
    """
    output_dist= []
    for el1, el2 in zip(dist1, dist2):
        out_el= function(el1, el2, cooling_table)
        output_dist.append(out_el)
    output_dist=np.array(output_dist).T[0]
    return output_dist


if __name__=='Main':
    print('running script directly')
    cooling_table, interp_kind=get_cooling_model_table()
    #mass_dist=operate_on_dist(np.random.normal(loc=5000., scale=250., size=int(default_n),np.random.normal(loc=8.0., scale=0.1., size=int(default_n)))


