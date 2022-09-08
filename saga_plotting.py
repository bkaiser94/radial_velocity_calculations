"""
Created by Ben Kaiser (UNC-Chapel Hill) 2022-09-08

Read in the "recommended" Milky Way dataset from the SAGA database of extremely metal poor stars and 
produce plots of various abundances using subsets of the sample.

This script will mostly be called by other scripts most likely in order to produce the backgrounds of plots.

Database URL: http://sagadatabase.jp

I'm going to write this to work in Python3 as its main goal. If it happens to be compatible with Python 2.7 
that's a happy accident.


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


RGB_max_teff=6000. #max temperature for RGB when paired with the max logg (to be RGB both criteria must be satisfied
RGB_max_logg=3.5

saga_file='SAGA_MetalPoor_txt_table_recommended.tsv'

saga_file=cp.abundance_dir+saga_file

converters={'*/H*': np.float64, 'd(*':np.float64}
saga_table=Table.read(saga_file, delimiter='\t',format='csv',converters=converters)
#saga_table=saga_table[:50]
#saga_table.pprint()

#saga_first_inds=np.where((saga_table['[Na/H]']-saga_table['[Fe/H]'])>0)
#saga_table=saga_table[saga_first_inds]


def plot_el_abunds(x,y,markersize=4, alpha=1, color='k', errorbars=True, MS_only=True,require_uncertainties=True):
    """
    Should be solar-normalized abundances that you want plotted. I'm going to default to using 
    the abundances that are not for specific ionization states. Not sure if this is the right call or 
    not though
    
    Inputs____
    x: abundance for the x-axis. Does not have to be an existing abundance in the data table, but 
    must have its components in the data table
    
    y:abundance for the y-axis. Does not have to be an existing abundance in the data table, but 
    must have its components in the data table
    
    MS_only: boolean of whether to implement the cuts specified by the SAGA database on what 
    counts as RGB vs. Main sequence.
    
    
    
    """
    print("MS_only:", MS_only)
    
    if MS_only:
        saga_subinds=np.where((saga_table['Teff']>RGB_max_teff) | (saga_table["log g"] > RGB_max_logg))
        saga_subtable=saga_table[saga_subinds]
    else:
        saga_subtable=saga_table
    
    def get_el_vals(el_ratio):
        #need to handle Li though too...
        try:
            print(saga_subtable[el_ratio][0])
            el_vals=saga_subtable[el_ratio]
            err_el_string=el_ratio.split('/')[0].replace('[','')
            el_err=saga_subtable['d('+err_el_string+')']
        except KeyError as error:
            print("KeyError:",error)
            el1, el2=el_ratio.split('/')
            el2=el2.replace(']','')
            el2='['+el2
            try:
                el1_H=saga_subtable[el1+'/H]']
            except KeyError as error:
                print('el1',el1,'hopefully is Li')
                el1_H=saga_subtable['A('+el1.replace('[','')+')']-3.26#the solar A(Li)=3.26 from Asplund et al. (2009) cited by SAGA as its reference value to get [Li/H]
            try:
                el2_H=saga_subtable[el2+'/H]']
            except KeyError as error:
                print('el2',el2,'hopefully is Li')
                el2_H=saga_subtable['A('+el2.replace('[','')+')']-3.26#subtract the solar A(Li)=3.26 from Asplund et al. (2009) cited by SAGA as its reference value to get [Li/H]
            #nans_el1=np.sum(np.isnan(el1))
            #nans_el2=np.sum(np.isnan(el2))
            el_vals=el1_H-el2_H
            #nans_difference=np.sum(np.isnan(el_vals))
            #print(nans_el1,nans_el2,nans_difference)
            el1_err=saga_subtable['d('+el1.replace('[','')+')']
            el2_err=saga_subtable['d('+el2.replace('[','')+')']
            #print('el1_err', el1_err)
            #print('el2_err',el2_err)
            
            #plt.hist(el1_err, alpha=0.2,label=el1,bins=np.arange(0,3,0.01))
            #plt.hist(el2_err,alpha=0.2,label=el2,bins=np.arange(0,3,0.01))
            #plt.legend()
            #plt.show()
            el_err=np.sqrt(el1_err**2.+el2_err**2.)
            #print('el_err', el_err)
        return el_vals,el_err
    #try:
        #print(saga_subtable[x][0])
        #x_vals=saga_subtable[x]
    #except KeyError as error:
        #print("KeyError:",error)
        #x_el1, x_el2=x.split('/')
        #x_el1_H=saga_subtable[x_el1+'/H]']
        #x_el2_H=saga_subtable[x_el2+'/H]']
        #x_vals=x_el1_H-x_el2_H
        
    x_vals,x_errs=get_el_vals(x)
    y_vals,y_errs=get_el_vals(y)
    
   
    if (require_uncertainties or errorbars):
        #require only points that have uncertainties for all elements can be included, written to impose this condition if you're plotting errorbars so no points without errorbars on both axes can be plotted.
        x_vals.mask=x_errs.mask
        y_vals.mask=y_errs.mask
    else:
        pass
    
    #print('max',x,np.nanmax(x_vals))
    #print(saga_subtable[np.argmax(x_vals)])
    #print('max',y,np.nanmax(y_vals))
    #print(saga_subtable[np.argmax(y_vals)])
    if errorbars:
        
        #plot errorbars now
        plt.errorbar(x_vals,y_vals,xerr=x_errs, yerr=y_errs, linestyle='none',marker='o',markersize=markersize,color=color,alpha=alpha,markeredgewidth=0)
    else:
        plt.errorbar(x_vals,y_vals,linestyle='none',marker='o',markersize=markersize,color=color,alpha=alpha,markeredgewidth=0)
    plt.xlabel(x)
    plt.ylabel(y)
    #plt.show()
    
    return



if __name__=='__main__':
    plot_el_abunds('[Fe/H]', '[Na/Ca]',MS_only=False,color='r')
    plot_el_abunds('[Fe/H]', '[Na/Ca]',MS_only=True)
    #plt.show()
    
    
    plot_el_abunds('[Fe/H]', '[Na/Ca]',MS_only=False,color='b',errorbars=False,alpha=0.4)
    plot_el_abunds('[Fe/H]', '[Na/Ca]',MS_only=True,errorbars=False,alpha=0.4,color='g')
    plt.show()
    
    #plot_el_abunds('[Fe/H]', '[Na/Ca]',MS_only=True,errorbars=False)
    #plot_el_abunds('[Fe/H]', '[Li/Fe]',MS_only=False)

    #plot_el_abunds('[Li/Ca]', '[Na/Ca]',MS_only=False)
    #plot_el_abunds('[Na/Ca]','Li/Ca]', MS_only=False)
    
    
    plot_el_abunds('[Fe/H]', '[Na/Fe]',MS_only=False,color='r')
    plot_el_abunds('[Fe/H]', '[Na/Fe]',MS_only=True)
    plt.show()
    
    plot_el_abunds('[Na/H]', '[Na/Fe]',MS_only=True)
    plt.show()
    
    plot_el_abunds('[Fe/H]', '[Ca/Fe]',MS_only=False,color='r')
    plot_el_abunds('[Fe/H]', '[Ca/Fe]',MS_only=True)
    plt.show()
    
else:
    pass














