"""
Created by Ben Kaiser (UNC-Chapel Hill) 2022-06-08

Import nuclear reaction cross-sections (empirical excitation functions) from Read & Viola 1984 
Table VIII and make them available to use in computing relative abundances of nuclear reaction 
products in a Python-friendly format. Probably astropy tables.




"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
import time

import cal_params as cp





def make_cross_table(projectile='p', target='C12'):
    filename=projectile+'_'+target+'.csv'
    filename=cp.nuke_cross_dir+filename
    full_table=Table.read(filename)
    
    return full_table


projectile_list=['p','He4']
target_list=['C12', 'N14', 'O16']
product_list=['A6', 'A7','A9','A10','A11']

for projectile in projectile_list:
    for target in target_list:
        cross_table=make_cross_table(projectile=projectile, target=target)
        for product in product_list:
            #plt.plot(cross_table['E_MeV'], cross_table[product],label=r''+projectile+'+'+target+r'$\rightarrow$'+product)
            if projectile=='p':
                plt.plot(cross_table['E_MeV'], cross_table[product],label=product)
                plt.xlabel('E (MeV)')
            else:
                plt.plot(cross_table['E_MeV']/4., cross_table[product],label=product)
                plt.xlabel('E/nucleon (MeV)')
        plt.title(projectile+'+'+target)
        plt.ylabel('Cross Section (mb)')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()
        plt.show()
