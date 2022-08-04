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



projectile_list=['p','He4']
target_list=['C12', 'N14', 'O16']
product_list=['A6', 'A7','A9','A10','A11']


def make_cross_table(projectile='p', target='C12'):
    filename=projectile+'_'+target+'.csv'
    filename=cp.nuke_cross_dir+filename
    full_table=Table.read(filename)
    
    return full_table

p_C12_table=make_cross_table(projectile='p', target='C12')
p_N14_table=make_cross_table(projectile='p', target='N14')
p_O16_table=make_cross_table(projectile='p', target='O16')
He4_C12_table=make_cross_table(projectile='He4', target='C12')
He4_N14_table=make_cross_table(projectile='He4', target='N14')
He4_O16_table=make_cross_table(projectile='He4', target='O16')

def select_cross_table(projectile='p', target='C12'):
    if projectile=='p':
        if target=='C12':
            cross_table=p_C12_table
        elif target=='N14':
            cross_table=p_N14_table
        elif target=='O16':
            cross_table=p_O16_table
        else:
            print('target not recognized', target)
    elif projectile=='He4':
        if target=='C12':
            cross_table=He4_C12_table
        elif target=='N14':
            cross_table=He4_N14_table
        elif target=='O16':
            cross_table=He4_O16_table
        else:
            print('target not recognized', target)
    else:
        print('projectile not recognized', projectile)
    
    return cross_table


def get_cross_section(energy, projectile='p', target='C12', product='A6'):
    """
    energy: float value or array in MeV (but without units attached)
    
    """
    print('Getting cross section of '+projectile+' + ' + target+ ' -> ' + product)
    cross_table=select_cross_table(projectile=projectile, target=target)
    cross_section=np.interp(energy, cross_table['E_MeV'], cross_table[product])
    
    return cross_section

def get_rel_cross_section(energy, projectile='p', target='C12', product1='A6', product2="A7"):
    """
    energy: float value or array in MeV (but without units attached)
    
    returns cross section for product 1 over product 2.
    
    
    """
    
    cross_table=select_cross_table(projectile=projectile, target=target)
    cross_section1=np.interp(energy, cross_table['E_MeV'], cross_table[product1])
    cross_section2=np.interp(energy, cross_table['E_MeV'], cross_table[product2])
    return cross_section1/cross_section2


if __name__ =='__main__':
    energy_range=np.linspace(0,1e4,int(1e5))
    for projectile in projectile_list:
        for target in target_list:
            A7_over_A9=get_rel_cross_section(energy_range, projectile=projectile, target=target, product1='A7', product2='A9')
            plt.plot(energy_range, A7_over_A9, label=projectile+'+'+target)
    plt.xlabel('Energy (MeV)')
    plt.legend()
    plt.ylabel(r'$\sigma_{Li-7}/\sigma_{Be-9}$')
    plt.title('Relative Cross sections of Li-7 and Be-9 (assuming all 7-mass isotopes are Li-7 and all 9-mass isotopes are Be-9')
    plt.axhline(y=1, linestyle='--',color='k')
    plt.axvline(x=1,linestyle='--',color='k')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

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
