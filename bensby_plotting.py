"""
Created by Ben Kaiser (UNC-Chapel Hill) 2020-03-02

Import the pre-joined Bensby et al. 2014 and Bensby and Lind 2018 table of abundances and then generate arrays to be plotted.

Probably also will just have a function to plot abundances that are requested. I don't know.

"""
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
import periodictable as pt
from astropy import units as u
from astropy import constants as const

import cal_params as cp


bensby_file= 'Bensby_2014_2018_abunds.csv'

bensby_file=cp.abundance_dir+bensby_file
lodders_file=cp.abundance_dir+'Lodders2009_solarsystem_abundances.csv'

bensby_table=Table.read(bensby_file)
lodders_table=Table.read(lodders_file)

lodders_table.add_index('element')

def get_rel_abund(el1, el2):
    """
    el1 and el2 as strings that need to match the element capitalization as well.
    get a new relative abundance from the [el/Fe] abundances that are contained in the table by default.
    
    This cannot be done with lithium. Nor can it be done with Iron. Actually I suppose it could be done with Iron, but the answer will be 0 for every single star since they'll be [Fe/H] minus [Fe/H].
    """
    
    el1el2= bensby_table[el1+'/Fe']-bensby_table[el2+'/Fe']
    el1el2_err= np.sqrt(bensby_table['e_'+el1+'/Fe']**2+bensby_table['e_'+el2+'/Fe']**2)
    return el1el2, el1el2_err

def get_ages():
    low_error=bensby_table['Age']- bensby_table['b_Age']
    hi_error=bensby_table['B_Age']-bensby_table['Age']
    return bensby_table['Age'], [low_error,hi_error]



def plot_el1el2_FeH(el1,el2,error_bars=True):
    el1el2, el1el2_err= get_rel_abund(el1,el2)
    if error_bars:
        plt.errorbar(bensby_table['Fe/H'], el1el2, xerr=bensby_table['e_Fe/H'], yerr=el1el2_err, linestyle='None', capsize=0, marker='o')
    else:
        plt.plot(bensby_table['Fe/H'], el1el2,  linestyle='None',  marker='o')
    plt.ylabel('['+el1+'/'+el2+']')
    plt.xlabel('[Fe/H]')
    return
    
    
def plot_el1el2_age(el1, el2, error_bars=True, mask_err_free=True):
    age, age_err=get_ages()
    el1el2, el1el2_err= get_rel_abund(el1,el2)
    if mask_err_free:
        #el_mask=np.where(~np.isnan(el1el2_err))
        el_mask=el1el2_err.mask
        print('el_mask',el_mask)
        print(type(el1el2_err))
        age.mask=el_mask
        age_err[0].mask=el_mask
        age_err[1].mask=el_mask
        #el1el2_err=el1el2_err
        #age_err[0]=age_err[0][el_mask]
        #age_err[1]=age_err[1][el_mask]
        #el1el2_err=el1el2_err[el_mask]
    else:
        pass
    if error_bars:
        plt.errorbar(age, el1el2, xerr=age_err, yerr=el1el2_err, linestyle='None', capsize=0, marker='o')
    else:
        plt.plot(age, el1el2,  linestyle='None',  marker='o')
    plt.ylabel('['+el1+'/'+el2+']')
    plt.xlabel('Age (Gyr)')
    return

def get_lica():
    lica=bensby_table['ALi']-bensby_table['Ca/Fe']-bensby_table['Fe/H']-lodders_table.loc['Li']['A_el']
    return lica

def plot_lica_FeH_pop():
    lica=get_lica()
    thick_disk_stars=np.where((bensby_table['td/d']>1) & (bensby_table['td/h']> 1))

    plt.plot(bensby_table['Fe/H'][thick_disk_stars], lica[thick_disk_stars], label='Bensby Thick Disk', linestyle='None', marker='o')
    
    thick_disk_stars=np.where(bensby_table['td/d']<1)

    plt.plot(bensby_table['Fe/H'][thick_disk_stars], lica[thick_disk_stars], label='Bensby Thin Disk', linestyle='None', marker='o')
    
    thick_disk_stars=np.where(bensby_table['td/h']<1)

    plt.plot(bensby_table['Fe/H'][thick_disk_stars], lica[thick_disk_stars], label='Bensby Halo', linestyle='None', marker='o' )
    plt.xlabel('[Fe/H]')
    plt.ylabel('[Li/Ca]')

    return


def plot_lica_age_pop():
    lica=get_lica()
    thick_disk_stars=np.where((bensby_table['td/d']>1) & (bensby_table['td/h']> 1))

    plt.plot(bensby_table['Age'][thick_disk_stars], lica[thick_disk_stars], label='Bensby Thick Disk', linestyle='None', marker='o')
    
    thick_disk_stars=np.where(bensby_table['td/d']<1)

    plt.plot(bensby_table['Age'][thick_disk_stars], lica[thick_disk_stars], label='Bensby Thin Disk', linestyle='None', marker='o')
    
    thick_disk_stars=np.where(bensby_table['td/h']<1)

    plt.plot(bensby_table['Age'][thick_disk_stars], lica[thick_disk_stars], label='Bensby Halo', linestyle='None', marker='o' )
    plt.xlabel('Age (Gyr)')
    plt.ylabel('[Li/Ca]')

    return

if __name__ == '__main__':
    #bensby_table.pprint()

    naca, naca_err= get_rel_abund('Na', 'Ca')
    lica=get_lica()

    #plt.errorbar(bensby_table['Fe/H'], naca, xerr=bensby_table['e_Fe/H'], yerr=naca_err, linestyle='None', capsize=0, marker='o')
    #plt.show()

    #plot_el1el2_FeH('Na','Ca', error_bars=True)
    #plt.show()

    #plot_el1el2_FeH('Na','Ca', error_bars=False)
    #plt.show()

    age, age_err=get_ages()

    #plt.errorbar(age, naca, xerr=age_err, yerr=naca_err, linestyle='None', capsize=0, marker='o')
    #plt.show()

    #plot_el1el2_age('Na', 'Ca',error_bars=False)
    #plot_el1el2_age('Na', 'Ca', mask_err_free=False)
    #plt.show()

    #plt.scatter(bensby_table['td/d'], bensby_table['td/h'])
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.show()

    #plt.scatter(bensby_table['td/d'], naca)
    #plt.xlabel('P(thick disk)/P(thin disk)')
    #plt.ylabel('[Na/Ca]')
    #plt.xscale('log')
    #plt.show()


    #plt.scatter(bensby_table['td/d'], bensby_table['Ca/Fe'])
    #plt.xlabel('P(thick disk)/P(thin disk)')
    #plt.ylabel('[Ca/Fe]')
    #plt.xscale('log')
    #plt.show()

    thick_disk_stars=np.where((bensby_table['td/d']>1) & (bensby_table['td/h']> 1))



    plt.plot(bensby_table['Fe/H'][thick_disk_stars], bensby_table['Ca/Fe'][thick_disk_stars], label='Thick Disk', linestyle='None', marker='o')
    #plt.xlabel('[Fe/H]')
    #plt.ylabel('[Ca/Fe]')
    #plt.show()

    thick_disk_stars=np.where(bensby_table['td/d']<1)



    plt.plot(bensby_table['Fe/H'][thick_disk_stars], bensby_table['Ca/Fe'][thick_disk_stars], label='Thin Disk', linestyle='None', marker='o')
    #plt.xlabel('[Fe/H]')
    #plt.ylabel('[Ca/Fe]')
    #plt.show()


    thick_disk_stars=np.where(bensby_table['td/h']<1)



    plt.plot(bensby_table['Fe/H'][thick_disk_stars], bensby_table['Ca/Fe'][thick_disk_stars], label='Halo', linestyle='None', marker='o' )
    plt.xlabel('[Fe/H]')
    plt.ylabel('[Ca/Fe]')
    #plt.yscale('log')
    #plt.xscale('log')
    plt.legend()
    plt.show()


    thick_disk_stars=np.where((bensby_table['td/d']>1) & (bensby_table['td/h']> 1))



    plt.plot(bensby_table['Fe/H'][thick_disk_stars], lica[thick_disk_stars], label='Thick Disk', linestyle='None', marker='o')
    #plt.xlabel('[Fe/H]')
    #plt.ylabel('[Ca/Fe]')
    #plt.show()

    thick_disk_stars=np.where(bensby_table['td/d']<1)



    plt.plot(bensby_table['Fe/H'][thick_disk_stars], lica[thick_disk_stars], label='Thin Disk', linestyle='None', marker='o')
    #plt.xlabel('[Fe/H]')
    #plt.ylabel('[Ca/Fe]')
    #plt.show()


    thick_disk_stars=np.where(bensby_table['td/h']<1)



    plt.plot(bensby_table['Fe/H'][thick_disk_stars], lica[thick_disk_stars], label='Halo', linestyle='None', marker='o' )
    plt.xlabel('[Fe/H]')
    plt.ylabel('[Li/Ca]')
    #plt.yscale('log')
    #plt.xscale('log')
    plt.legend()
    plt.show()

    plot_lica_FeH_pop()
    plt.legend()
    plt.show()

    plot_lica_age_pop()
    plt.legend()
    plt.show()
