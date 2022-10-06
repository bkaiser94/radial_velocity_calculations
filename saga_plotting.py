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
Crich_cut=0.7
EMP_cut=-2.5
CEMPs_cut=0.5





classifier_dict={
    'MP':['teal','P'],
    'EMP': ['gold','*'],
    'CEMP':['purple','o'],
    'C-rich':['b','s'],
    'CEMP-s':['salmon','^'],
    'CEMP-no':['k','v']}
 

saga_file='SAGA_MetalPoor_txt_table_recommended.tsv'
#saga_file='SAGA_MetalPoor_txt_table_all.tsv'


saga_file=cp.abundance_dir+saga_file

converters={'*/H*': np.float64, 'd(*':np.float64}
saga_table=Table.read(saga_file, delimiter='\t',format='csv',converters=converters,encoding='latin')
#saga_table=Table.read(saga_file, delimiter='\t',format='csv')

#saga_table=saga_table[:50]
#saga_table.pprint()

#saga_first_inds=np.where((saga_table['[Na/H]']-saga_table['[Fe/H]'])>0)
#saga_table=saga_table[saga_first_inds]


def get_el_vals(el_ratio, saga_subtable_main=saga_table):
        #need to handle Li though too...
        saga_subtable=saga_subtable_main.copy()
        def check_all_ions(el, saga_subtable):
            el_H=saga_subtable[el+'/H]']
            el_err=saga_subtable['d('+el.replace('[','')+')']
            missing_vals=np.where(el_H.mask==True)
            el_i_H=saga_subtable[missing_vals][el+' I/H]']
            saga_subtable[missing_vals][el+'/H]']=el_i_H
            return
        
        
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
                el1_err=saga_subtable['d('+el1.replace('[','')+')']
            except KeyError as error:
                print('el1',el1,'hopefully is Li')
                el1_H=saga_subtable['A('+el1.replace('[','')+')']-3.26#the solar A(Li)=3.26 from Asplund et al. (2009) cited by SAGA as its reference value to get [Li/H]
                el1_err=saga_subtable['d('+el1.replace('[','')+')']
            try:
                el2_H=saga_subtable[el2+'/H]']
                el2_err=saga_subtable['d('+el2.replace('[','')+')']
            except KeyError as error:
                print('el2',el2,'hopefully is Li')
                el2_H=saga_subtable['A('+el2.replace('[','')+')']-3.26#subtract the solar A(Li)=3.26 from Asplund et al. (2009) cited by SAGA as its reference value to get [Li/H]
                el2_err=saga_subtable['d('+el2.replace('[','')+')']
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

def get_subtable(MS_only=True,classifier='EMP',saga_table=saga_table):
    if MS_only:
        saga_subinds=np.where((saga_table['Teff']>RGB_max_teff) | (saga_table["log g"] > RGB_max_logg))
        saga_subtable=saga_table[saga_subinds].copy()
    else:
        saga_subtable=saga_table.copy()
    CFe,throwaway=get_el_vals('[C/Fe]',saga_subtable_main=saga_subtable)
    if classifier=='MP':
        saga_MP_inds=np.where((saga_subtable['[Fe/H]']>EMP_cut)& (CFe < Crich_cut))
        saga_output_table=saga_subtable[saga_MP_inds]
    elif classifier=='EMP':
        saga_EMP_inds=np.where((saga_subtable['[Fe/H]']<=EMP_cut) & (CFe < Crich_cut))
        saga_output_table=saga_subtable[saga_EMP_inds]
    elif classifier=='C-rich':
        saga_Crich_inds=np.where((saga_subtable['[Fe/H]']>EMP_cut) & (CFe >= Crich_cut))
        saga_output_table=saga_subtable[saga_Crich_inds]
    elif classifier[:4]=='CEMP':
        saga_CEMP_inds=np.where((saga_subtable['[Fe/H]']<=EMP_cut) & (CFe >= Crich_cut))
        saga_CEMP_table=saga_subtable[saga_CEMP_inds]
        BaFe,throway=get_el_vals('[Ba/Fe]',saga_subtable_main=saga_CEMP_table)
        if classifier=='CEMP-s':
            saga_CEMPs_inds=np.where(BaFe>=CEMPs_cut)
            saga_output_table=saga_CEMP_table[saga_CEMPs_inds]
        elif classifier=='CEMP-no':
            saga_CEMPno_inds=np.where(BaFe<CEMPs_cut)
            saga_output_table=saga_CEMP_table[saga_CEMPno_inds]
        elif classifier=='CEMP':
            saga_CEMP_subinds=np.where(BaFe.mask==True)
            saga_output_table=saga_CEMP_table[saga_CEMP_subinds]
        else:
            print('CEMP subclassifier:', classifier, 'not recognized')
    else:
        print('classifier:', classifier, 'unrecognized')
    
    
    
    return saga_output_table


def plot_el_abunds(x,y,markersize=4, alpha=1, color='k', errorbars=True, MS_only=True,require_uncertainties=True,saga_table=saga_table,marker='o'):
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
    
    
    #try:
        #print(saga_subtable[x][0])
        #x_vals=saga_subtable[x]
    #except KeyError as error:
        #print("KeyError:",error)
        #x_el1, x_el2=x.split('/')
        #x_el1_H=saga_subtable[x_el1+'/H]']
        #x_el2_H=saga_subtable[x_el2+'/H]']
        #x_vals=x_el1_H-x_el2_H
        
    x_vals,x_errs=get_el_vals(x,saga_subtable_main=saga_subtable)
    y_vals,y_errs=get_el_vals(y, saga_subtable_main=saga_subtable)
    
   
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
        plt.errorbar(x_vals,y_vals,xerr=x_errs, yerr=y_errs, linestyle='none',marker=marker,markersize=markersize,color=color,alpha=alpha,markeredgewidth=0)
    else:
        plt.errorbar(x_vals,y_vals,linestyle='none',marker=marker,markersize=markersize,color=color,alpha=alpha,markeredgewidth=0)
    plt.xlabel(x)
    plt.ylabel(y)
    #plt.show()
    
    return


def plot_class_el_abunds(x,y,markersize=4, alpha=1, color='k', errorbars=True, MS_only=True,require_uncertainties=True,saga_table=saga_table,use_class_color=True):
    """
    iterate through all of the classifiers and plot each of them for the given abundances
    
    """
    for classifier in classifier_dict:
        marker=classifier_dict[classifier][1]
        if use_class_color==True:
            color=classifier_dict[classifier][0]
        else:
            pass
        print('classifier:', classifier, 'color:', color)
        saga_subtable=get_subtable(MS_only=MS_only,classifier=classifier,saga_table=saga_table)
        plot_el_abunds(x,y,markersize=markersize, marker=marker,alpha=alpha, color=color, errorbars=errorbars, MS_only=MS_only,require_uncertainties=require_uncertainties,saga_table=saga_subtable)
        #plt.show()
    return 



def plot_el_vs_param(x,y,markersize=4, alpha=1, color='k', errorbars=True, MS_only=True,require_uncertainties=True):
    """
    Should be plotting [el/H] vs. Teff for example
    
    """
    print("MS_only:", MS_only)
    
    if MS_only:
        saga_subinds=np.where((saga_table['Teff']>RGB_max_teff) | (saga_table["log g"] > RGB_max_logg))
        saga_subtable=saga_table[saga_subinds]
    else:
        saga_subtable=saga_table
    
    
    #try:
        #print(saga_subtable[x][0])
        #x_vals=saga_subtable[x]
    #except KeyError as error:
        #print("KeyError:",error)
        #x_el1, x_el2=x.split('/')
        #x_el1_H=saga_subtable[x_el1+'/H]']
        #x_el2_H=saga_subtable[x_el2+'/H]']
        #x_vals=x_el1_H-x_el2_H
        
    x_vals=saga_subtable[x]
    y_vals,y_errs=get_el_vals(y,saga_subtable_main=saga_subtable)
    if (require_uncertainties or errorbars):
        #require only points that have uncertainties for all elements can be included, written to impose this condition if you're plotting errorbars so no points without errorbars on both axes can be plotted.
        y_vals.mask=y_errs.mask
    else:
        pass
    
    #print('max',x,np.nanmax(x_vals))
    #print(saga_subtable[np.argmax(x_vals)])
    #print('max',y,np.nanmax(y_vals))
    #print(saga_subtable[np.argmax(y_vals)])
    if errorbars:
        
        #plot errorbars now
        plt.errorbar(x_vals,y_vals, yerr=y_errs, linestyle='none',marker='o',markersize=markersize,color=color,alpha=alpha,markeredgewidth=0)
    else:
        plt.errorbar(x_vals,y_vals,linestyle='none',marker='o',markersize=markersize,color=color,alpha=alpha,markeredgewidth=0)
    plt.xlabel(x)
    plt.ylabel(y)
    
    
    return

if __name__=='__main__':
    plot_class_el_abunds('[Ca/H]','[Na/Ca]',markersize=8, alpha=1,errorbars=False, MS_only=True,require_uncertainties=True,saga_table=saga_table,use_class_color=True)
    #plot_el_abunds('[Ca/H]', '[Na/Ca', MS_only=True,errorbars=False,require_uncertainties=False, alpha=1,color='b')
    plt.show()
    
    plot_class_el_abunds('[Ca/H]','[Ca/Fe]',markersize=8, alpha=1,errorbars=False, MS_only=True,require_uncertainties=True,saga_table=saga_table,use_class_color=True)
    #plot_el_abunds('[Ca/H]', '[Na/Ca', MS_only=True,errorbars=False,require_uncertainties=False, alpha=1,color='b')
    plt.show()
    
    plot_class_el_abunds('[Fe/H]','[Ca/Fe]',markersize=8, alpha=1,errorbars=False, MS_only=True,require_uncertainties=True,saga_table=saga_table,use_class_color=True)
    #plot_el_abunds('[Ca/H]', '[Na/Ca', MS_only=True,errorbars=False,require_uncertainties=False, alpha=1,color='b')
    plt.show()

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














