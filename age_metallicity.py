"""
Created by Ben Kaiser (UNC-Chapel Hill) 2020-01-09


Take [Fe/H] and output [Ca/Fe] (and any other Z/Fe desired ratio) to be used for a given conversion. Hopefully in 
the future this will contain an actual age-metallicity relation additionally so it can all be self-contained,
but alas I dream...


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


GALAH_thick_mfile='Lin2020_GALAH_thickdisc_m_clean.txt'

GALAH_thin_mfile='Lin2020_GALAH_thindisc_m_clean.txt'

#GALAH_thick_mfile=cp.abundance_dir+GALAH_thick_mfile
#GALAH_thin_mfile=cp.abundance_dir+GALAH_thin_mfile




print('\n\n****************\nTHIS IS WRONG!!!!')
print('The Lin et al. 2020 functions from GALAH are supposed to be [X/Fe] vs. "Age" (Gyr), but I made this script plot versus [Fe/H], which is not what the coefficients were for.\n***********************\n\n')



def get_b_file(m_file):
    return m_file.replace('m','b')

def metal_function(FeH, element,m_table, b_table):
    bound_lists= m_table.colnames[1:]
    new_list=[]
    element_m = m_table.loc[element]
    element_b= b_table.loc[element]
    print(element_m)
    print(element_b)
    is_array=False
    try:
        output_val=np.ones(FeH.shape)
        is_array=True
    except AttributeError:
        pass
    for thing in bound_lists:
        thing=thing.replace('[', '')
        thing=thing.replace(']','')
        thing1, thing2= thing.split(',')
        thing1=float(thing1)
        thing2=float(thing2)
        #new_list.append(list(thing))
        new_list.append([thing1,thing2])
    for number, bound in enumerate(new_list):
        m_val= element_m[bound_lists[number]]
        b_val=element_b[bound_lists[number]]
        if is_array:
            inplay=np.where((FeH > bound[0])& (FeH <= bound[1]))
            output_val[inplay]=FeH[inplay]*m_val+b_val
        else:
            if ((FeH > bound[0])& (FeH <= bound[1])):
                output_val=FeH*m_val+b_val
            else:
                pass
        
    return output_val



def get_metal_vals(FeH,m_file, element='Ca'):
    b_file=get_b_file(m_file)
    m_table=Table.read(cp.abundance_dir+m_file, format='ascii')
    b_table=Table.read(cp.abundance_dir+b_file, format='ascii')
    m_table.add_index(m_table.colnames[0])#making Element a searchable index in the table.
    b_table.add_index(b_table.colnames[0]) #making Element a searchable index in the table.
    metal_ratios= metal_function(FeH, element, m_table, b_table)
    return metal_ratios



FeH_vals= np.linspace(-0.49, 0.5, 100)
thick_CaFe= get_metal_vals(FeH_vals, GALAH_thick_mfile, element='Ca')
thin_CaFe = get_metal_vals(FeH_vals, GALAH_thin_mfile, element='Ca')



plt.plot(FeH_vals, thin_CaFe, label='Thin Disk')
plt.plot(FeH_vals, thick_CaFe, label='Thick Disk')
plt.xlabel('[Fe/H]')
plt.ylabel('[Ca/Fe]')
plt.legend()
plt.show()



print('\n\n****************\nTHIS IS WRONG!!!!')
print('The Lin et al. 2020 functions from GALAH are supposed to be [X/Fe] vs. "Age" (Gyr), but I made this script plot versus [Fe/H], which is not what the coefficients were for.\n***********************\n\n')



