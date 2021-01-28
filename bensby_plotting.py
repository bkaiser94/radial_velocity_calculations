"""
Created by Ben Kaiser (UNC-Chapel Hill) 2020-03-02

Import the pre-joined Bensby et al. 2014 and Bensby and Lind 2018 table of abundances and then generate arrays to be plotted.

Probably also will just have a function to plot abundances that are requested. I don't know.

"""
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.table import join as ATjoin
from astropy.table import hstack as AThstack
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



#Bensby et al. 2014 probability cuts for populations. From their paper.

thin_disk_bound_tdd= 0.5 #below this value is likely thin disk
thick_disk_bound_tdd=2.0 #potential thick disk stars
halo_bound_tdh= 1. #only cut here is TD/H < 1 for likely halo. That's it.

inbetween_range_tdd= [thin_disk_bound_tdd, thick_disk_bound_tdd]

def get_rel_abund(el1, el2):
    """
    el1 and el2 as strings that need to match the element capitalization as well.
    get a new relative abundance from the [el/Fe] abundances that are contained in the table by default.
    
    This cannot be done with lithium. Nor can it be done with Iron. Actually I suppose it could be done with Iron, but the answer will be 0 for every single star since they'll be [Fe/H] minus [Fe/H].
    """
    try:
        el1el2= bensby_table[el1+'/Fe']-bensby_table[el2+'/Fe']
        el1el2_err= np.sqrt(bensby_table['e_'+el1+'/Fe']**2+bensby_table['e_'+el2+'/Fe']**2)
    except KeyError as error:
        print('\n\n****************\nKeyError:',error,'\n****************\n\n')
        el1el2= bensby_table[el1+'/Fe']
        el1el2_err= bensby_table['e_'+el1+'/Fe']
    return el1el2, el1el2_err

def get_ages():
    low_error=bensby_table['Age']- bensby_table['b_Age']
    hi_error=bensby_table['B_Age']-bensby_table['Age']
    return bensby_table['Age'], [low_error,hi_error]

def get_ALi():
    low_error=bensby_table['ALi']- bensby_table['b_ALi']
    hi_error=bensby_table['B_ALi']-bensby_table['ALi']
    print(hi_error)
    print(type(hi_error))
    #return bensby_table['ALi'], np.array([low_error,hi_error])
    #return bensby_table['ALi'],[low_error.data,hi_error.data]
    return bensby_table['ALi'],low_error, hi_error





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

def get_lica(with_errors=False, sol_norm=True):
    if sol_norm:
        lica=bensby_table['ALi']-bensby_table['Ca/Fe']-bensby_table['Fe/H']-lodders_table.loc['Li']['A_el']
        if with_errors:
            pass
            ALi, ALi_lo_error, ALi_hi_error= get_ALi()
            def combine_errors(this_error):
                return np.sqrt(this_error**2+ bensby_table['e_Ca/Fe']**2+bensby_table['e_Fe/H']**2+lodders_table.loc['Li']['A_el_err']**2)
            lica_lo_error=combine_errors(ALi_lo_error)
            lica_hi_error=combine_errors(ALi_hi_error)
            lica=bensby_table['ALi']-bensby_table['Ca/Fe']-bensby_table['Fe/H']-lodders_table.loc['Li']['A_el']
            return lica, [lica_lo_error, lica_hi_error]
        else:
            return lica
    else:
        lica=bensby_table['ALi']-bensby_table['Ca/Fe']-bensby_table['Fe/H']-lodders_table.loc['Ca']['A_el']
        if with_errors:
            pass
            ALi, ALi_lo_error, ALi_hi_error= get_ALi()
            def combine_errors(this_error):
                return np.sqrt(this_error**2+ bensby_table['e_Ca/Fe']**2+bensby_table['e_Fe/H']**2+lodders_table.loc['Ca']['A_el_err']**2)
            lica_lo_error=combine_errors(ALi_lo_error)
            lica_hi_error=combine_errors(ALi_hi_error)
            lica=bensby_table['ALi']-bensby_table['Ca/Fe']-bensby_table['Fe/H']-lodders_table.loc['Ca']['A_el']
            return lica, [lica_lo_error, lica_hi_error]
        else:
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

def plot_FeH_age():
    star_age, star_age_error= get_ages()
    pop_id_array=np.int_(np.zeros(bensby_table['td/d'].shape) )#array to be comprised of numbers that represent the population to which each star should belong; going to be added to as we go here.
    
    thick_disk_stars=np.where(bensby_table['td/d']<=thin_disk_bound_tdd)
    FeH=bensby_table['Fe/H']
    FeH_error=bensby_table['e_Fe/H']
    pop_id_array[thick_disk_stars]=0
    #plt.plot(bensby_table['Age'][thick_disk_stars], lica[thick_disk_stars], label='Thin Disk', linestyle='None', marker=marker, color=colors[0])

    plt.errorbar(star_age[thick_disk_stars], FeH[thick_disk_stars],xerr=FeH_error[thick_disk_stars],yerr=star_age_error[thick_disk_stars], label='Thin Disk', linestyle='None', marker=marker, color=colors[0], )
    
    thick_disk_stars=np.where((bensby_table['td/d']>=thick_disk_bound_tdd) & (bensby_table['td/h']>= halo_bound_tdh))
    
    pop_id_array[thick_disk_stars]=1

    plt.errorbar(star_age[thick_disk_stars], FeH[thick_disk_stars],xerr=FeH_error[thick_disk_stars],yerr=star_age_error[thick_disk_stars], label='Thick Disk', linestyle='None', marker=marker, color=colors[1])
    
    
    thick_disk_stars=np.where(bensby_table['td/h']<halo_bound_tdh)
    
    pop_id_array[thick_disk_stars]=2


    plt.errorbar(star_age[thick_disk_stars], FeH[thick_disk_stars],xerr=FeH_error[thick_disk_stars],yerr=star_age_error[thick_disk_stars], label='Halo', linestyle='None', marker=marker, color=colors[2])
    
    thick_disk_stars=np.where((bensby_table['td/d']<thick_disk_bound_tdd)&(bensby_table['td/d']> thin_disk_bound_tdd))
    
    pop_id_array[thick_disk_stars]=3

    plt.errorbar(star_age[thick_disk_stars], FeH[thick_disk_stars],xerr=FeH_error[thick_disk_stars],yerr=star_age_error[thick_disk_stars], label='In Between', linestyle='None', marker=marker, color=colors[3])
    return


def plot_lica_age_pop(colors=['b','b','b','b'],marker='o', rep_errors=False, sol_norm=True,markersize=8):
    if rep_errors:
        lica, lica_error= get_lica(with_errors=True, sol_norm=sol_norm)
    else:
        lica=get_lica(sol_norm=sol_norm)
    star_age, star_age_error= get_ages()
    pop_id_array=np.int_(np.zeros(bensby_table['td/d'].shape) )#array to be comprised of numbers that represent the population to which each star should belong; going to be added to as we go here.
    
    thick_disk_stars=np.where(bensby_table['td/d']<=thin_disk_bound_tdd)
    
    pop_id_array[thick_disk_stars]=0
    #plt.plot(bensby_table['Age'][thick_disk_stars], lica[thick_disk_stars], label='Thin Disk', linestyle='None', marker=marker, color=colors[0])

    plt.plot(bensby_table['Age'][thick_disk_stars], lica[thick_disk_stars], label='Thin Disk', linestyle='None', marker=marker, color=colors[0],markersize=markersize) #with labels
    
    #plt.plot(bensby_table['Age'][thick_disk_stars], lica[thick_disk_stars], label='Neighborhood Stars', linestyle='None', marker=marker, color=colors[0],markersize=markersize) #single label
    
    thick_disk_stars=np.where((bensby_table['td/d']>=thick_disk_bound_tdd) & (bensby_table['td/h']>= halo_bound_tdh))
    
    pop_id_array[thick_disk_stars]=1

    plt.plot(bensby_table['Age'][thick_disk_stars], lica[thick_disk_stars], label='Thick Disk', linestyle='None', marker=marker, color=colors[1],markersize=markersize) #with lables
    
    #plt.plot(bensby_table['Age'][thick_disk_stars], lica[thick_disk_stars], linestyle='None', label='', marker=marker, color=colors[1],markersize=markersize) #single label
    
    
    thick_disk_stars=np.where(bensby_table['td/h']<halo_bound_tdh)
    
    pop_id_array[thick_disk_stars]=2


    plt.plot(bensby_table['Age'][thick_disk_stars], lica[thick_disk_stars], label='Halo', linestyle='None', marker=marker, color=colors[2],markersize=markersize) #with labels
    
    #plt.plot(bensby_table['Age'][thick_disk_stars], lica[thick_disk_stars],  label='', linestyle='None', marker=marker, color=colors[2],markersize=markersize) #single label
    
    thick_disk_stars=np.where((bensby_table['td/d']<thick_disk_bound_tdd)&(bensby_table['td/d']> thin_disk_bound_tdd))
    
    pop_id_array[thick_disk_stars]=3

    plt.plot(bensby_table['Age'][thick_disk_stars], lica[thick_disk_stars], label='In Between', linestyle='None', marker=marker, color=colors[3],markersize=markersize) #with labels
    
    #plt.plot(bensby_table['Age'][thick_disk_stars], lica[thick_disk_stars], label='', linestyle='None', marker=marker, color=colors[3],markersize=markersize) #single label
    
    #I'm going to hold off on making new versions of the plotting itself until I've determined it looks different conclusively
    #rep_range=[8,9.5] #Age range that the representative point will be pulled from; it will be the max [Li/Ca] from that group.
    rep_range=[3.5,4.5]
    
    if rep_errors:
        med_age_error=np.nanmedian(star_age_error, axis=1)
        med_lica_error=np.nanmedian(lica_error, axis=1)
        sorted_order= np.argsort(lica)
        indices= np.indices(lica.shape)[0] #need an array of indices because I'm going to have to cut it down and need to retain the original indices.
        print("star_age.shape", star_age.shape)
        in_range=np.where((star_age > rep_range[0]) & (star_age < rep_range[1]))
        print('in_range', in_range)
        print(indices.shape)
        
        #max_arg= np.argmax(lica[in_range])
        max_arg= np.argmin(lica[in_range]) #I realize this is actually the minimum now, but it's easier to not rename all of the following variables
        
        max_index= indices[in_range][max_arg]
        print('max_index', max_index)
        print('max_lica', lica[max_index])
        print('max age', star_age[max_index])
        print('med_age_error', med_age_error)
        print('med_lica_error', med_lica_error)
        print(med_age_error.shape)
        print(np.array(star_age[max_index]).shape)
        print(colors[pop_id_array[max_index]])
        #plt.errorbar(np.array([star_age[max_index]]), np.array([lica[max_index]]), xerr=np.array([med_age_error]), marker=marker, linestyle='None', color=colors[pop_id_array[max_index]])
        #plt.errorbar(star_age[max_index], lica[max_index], xerr=np.array([med_age_error]).T, marker=marker, linestyle='None', color=colors[pop_id_array[max_index]])
        plt.errorbar(star_age[max_index], lica[max_index], xerr=np.array([med_age_error]).T, yerr=np.array([med_lica_error]).T , marker=marker, linestyle='None', color=colors[pop_id_array[max_index]],markersize=markersize)
    else:
        pass
    
    
    plt.xlabel('Age (Gyr)')
    plt.ylabel('[Li/Ca]')

    return



def plot_ALi_FeH(colors=['b','b','b','b'],marker='o',markersize=8):
    #ALi, ALi_err= get_ALi()
    ALi, ALi_lo, ALi_hi= get_ALi()
    #ALi_err=ALi_err.T
    print(type(ALi))
    #print(type(ALi_err))
    #print(ALi)
    #print(ALi_err)
    print('ALi.shape', ALi.shape)
    #print('ALi_err.shape', ALi_err.shape)
    #thick_disk_stars=np.where(bensby_table['td/d']<=thin_disk_bound_tdd)

    #plt.errorbar(bensby_table['Fe/H'][thick_disk_stars], ALi[thick_disk_stars], xerr=bensby_table['e_Fe/H'][thick_disk_stars],yerr=ALi_err[thick_disk_stars],label='Thin Disk', linestyle='None', marker=marker, color=colors[0])
    
    #thick_disk_stars=np.where((bensby_table['td/d']>=thick_disk_bound_tdd) & (bensby_table['td/h']>= halo_bound_tdh))

    #plt.errorbar(bensby_table['Fe/H'][thick_disk_stars], ALi[thick_disk_stars], xerr=bensby_table['e_Fe/H'][thick_disk_stars],yerr=ALi_err[thick_disk_stars],label='Thick Disk', linestyle='None', marker=marker, color=colors[1])
    
    
    #thick_disk_stars=np.where(bensby_table['td/h']<halo_bound_tdh)

    #plt.errorbar(bensby_table['Fe/H'][thick_disk_stars], ALi[thick_disk_stars], xerr=bensby_table['e_Fe/H'][thick_disk_stars],yerr=ALi_err[thick_disk_stars],label='Halo', linestyle='None', marker=marker, color=colors[2])
    
    #thick_disk_stars=np.where((bensby_table['td/d']<thick_disk_bound_tdd)&(bensby_table['td/d']> thin_disk_bound_tdd))

    #plt.errorbar(bensby_table['Fe/H'][thick_disk_stars], ALi[thick_disk_stars], xerr=bensby_table['e_Fe/H'][thick_disk_stars],yerr=ALi_err[thick_disk_stars],label='In Between', linestyle='None', marker=marker, color=colors[3])
    
    #######
    
    
    thick_disk_stars=np.where((bensby_table['td/d']<thick_disk_bound_tdd)&(bensby_table['td/d']> thin_disk_bound_tdd))

    plt.errorbar(bensby_table['Fe/H'][thick_disk_stars], ALi[thick_disk_stars], xerr=bensby_table['e_Fe/H'][thick_disk_stars],yerr=[ALi_lo[thick_disk_stars], ALi_hi[thick_disk_stars]],label='In Between', linestyle='None', marker=marker, color=colors[3], markersize=markersize)
    thick_disk_stars=np.where(bensby_table['td/d']<=thin_disk_bound_tdd)

    #print(ALi_err[thick_disk_stars])


    plt.errorbar(bensby_table['Fe/H'][thick_disk_stars], ALi[thick_disk_stars], xerr=bensby_table['e_Fe/H'][thick_disk_stars],yerr=[ALi_lo[thick_disk_stars], ALi_hi[thick_disk_stars]],label='Thin Disk', linestyle='None', marker=marker, color=colors[0],markersize=markersize)
    
    thick_disk_stars=np.where((bensby_table['td/d']>=thick_disk_bound_tdd) & (bensby_table['td/h']>= halo_bound_tdh))

    plt.errorbar(bensby_table['Fe/H'][thick_disk_stars], ALi[thick_disk_stars], xerr=bensby_table['e_Fe/H'][thick_disk_stars],yerr=[ALi_lo[thick_disk_stars], ALi_hi[thick_disk_stars]],label='Thick Disk', linestyle='None', marker=marker, color=colors[1],markersize=markersize)
    
    
    thick_disk_stars=np.where(bensby_table['td/h']<halo_bound_tdh)

    plt.errorbar(bensby_table['Fe/H'][thick_disk_stars], ALi[thick_disk_stars], xerr=bensby_table['e_Fe/H'][thick_disk_stars],yerr=[ALi_lo[thick_disk_stars], ALi_hi[thick_disk_stars]],label='Halo', linestyle='None', marker=marker, color=colors[2],markersize=markersize)
    
    
    ###########
    #plt.errorbar(bensby_table['Fe/H'], ALi, xerr=bensby_table['e_Fe/H'], yerr=ALi_err, linestyle='None', marker='o', label='BL2018')
    
    plt.xlabel('[Fe/H]')
    plt.ylabel('A(Li)')
    return

if __name__ == '__main__':
    #bensby_table.pprint()
    
    
    plot_el1el2_FeH('Ca','Fe', error_bars=True)
    plt.show()
    
    lica, lica_error= get_lica(with_errors=True)
    plt.hist(lica_error[0], bins=20)
    plt.xlabel('[Li/Ca] error lower bound')
    plt.show()
    
    plt.hist(lica_error[1], bins=20)
    plt.xlabel('[Li/Ca] error upper bound')
    plt.show()
    
    plot_lica_age_pop(rep_errors=True)
    plt.show()
    
    
    naca, naca_err= get_rel_abund('Na', 'Ca')
    lica=get_lica()
    
    plot_ALi_FeH()
    plt.show()

    plt.errorbar(bensby_table['Fe/H'], naca, xerr=bensby_table['e_Fe/H'], yerr=naca_err, linestyle='None', capsize=0, marker='o')
    plt.show()

    plot_el1el2_FeH('Na','Ca', error_bars=True)
    plt.show()

    plot_el1el2_FeH('Na','Ca', error_bars=False)
    plt.show()

    age, age_err=get_ages()
    age_err_array=np.array(age_err)
    print('age_err_array.shape', age_err_array.shape)
    med_age_errs= np.nanmedian(age_err_array, axis=1)
    print('med_age_errs', med_age_errs)

    plt.errorbar(age, naca, xerr=age_err, yerr=naca_err, linestyle='None', capsize=0, marker='o')
    plt.xlabel('Age')
    plt.ylabel('[Na/Ca]')
    plt.show()
    
    plt.hist(age_err_array[0])
    plt.xlabel('lower errors on age')
    plt.show()
    
    plt.hist(age_err_array[1])
    plt.xlabel('upper errors on age')
    plt.show()

    plot_el1el2_age('Na', 'Ca',error_bars=False)
    plot_el1el2_age('Na', 'Ca', mask_err_free=False)
    plt.show()

    plt.scatter(bensby_table['td/d'], bensby_table['td/h'])
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

    plt.scatter(bensby_table['td/d'], naca)
    plt.xlabel('P(thick disk)/P(thin disk)')
    plt.ylabel('[Na/Ca]')
    plt.xscale('log')
    plt.show()


    plt.scatter(bensby_table['td/d'], bensby_table['Ca/Fe'])
    plt.xlabel('P(thick disk)/P(thin disk)')
    plt.ylabel('[Ca/Fe]')
    plt.xscale('log')
    plt.show()

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
