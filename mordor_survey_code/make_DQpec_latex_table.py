"""
Created by Ben Kaiser (UNC-Chapel Hill) 2023-03-03.
Modified to be used for DQpecs on 2025-06-26

Read in the MORDOR survey CSV formatted astropy table and output the desired columns with names reformatted and column names (if possible) as a LaTeX table or at least in a format that will make it relatively easy to make a LaTeX table.



"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
from glob import glob
from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u
from astropy import constants as const
from astropy import convolution as conv
from astropy.table import Table, Column
#import scipy.interpolate as scinterp
import time
import csv


#plt.rc('font', size =18)
sys.path.append('../')
sys.path.append('/Users/BenKaiser/Desktop/radial_velocity_calculations/')


#print start
#import wdatmos
import spec_plot_tools as spt
import cal_params as cp

input_file='full_MORDOR_survey_1681141339_fullphotometry_fitparams_massadd.csv'


input_table = Table.read(input_file)


#selected_output_cols=['name','ra','dec','phot_g_mean_mag','g_rp','g_abs','sp_type','sed_sp_type','pwd']
selected_output_cols=['name','dr3_source_id','phot_g_mean_mag','teff','teff_err','logg','logg_err','m_wd','m_wd_err','c/he','c/he_err']

colname_replacements={'name':'Name',
                      'dr3_source_id':'\Gaia\ DR3 ID',
                      'phot_g_mean_mag':'$G$ (mag)',
                      'teff':'\Teff',
                      'logg':'$\log(g)$',
                      'm_wd':'\MWD (\Msol)'
    }

noerr_cols=[]

for thing in selected_output_cols:
    if "_err" in thing:
        pass
    else:
        noerr_cols.append(thing)

print('noerr_cols',noerr_cols)

desired_sptype='DQpec'


current_directory=os.getcwd()
output_name_base=input_file.split('/')[-1].split('.')[0]+'_'+current_directory.split('/')[-1]

object_name_base='MORDOR J'


    
def make_name(this_table):
    this_table['name']=this_table['name'].astype('S32')
    for row in this_table:
        input_coords = coord.SkyCoord(ra = row['ra'], dec =row['dec'], unit = (u.deg, u.deg), frame = 'icrs')
        string_coords= input_coords.to_string(style='hmsdms')
        replace_chars= ['d','h','m']
        ra_list=[]
        dec_list=[]
        epoch_list=[]
        mag_list=[]
        name_list=[]

        #print(name)
        #print(thing)
        #for char in replace_chars:
            #thing= thing.replace(char, ':')
        #thing=thing.replace('s', '')
        #split_string=thing.split(' ')
        #ra= split_string[0][:11] #limiting precision of the decimal
        #dec= split_string[1][:12] #limiting precision of the decimal, need the additional index for sign
        #small_ra= ra.replace(':', '')[:4]
        #small_dec=dec.replace(':','')[:5] #need additional index for + or -
        if 'GaiaJ' in row['name']:
            print("Gaia detected", row['name'])
            row['name']=row['name'].replace('GaiaJ',object_name_base)
            row['name']=row['name'].replace('-','–')
        elif 'SDSSJ' in row['name']:
            print("SDSS detected", row['name'])
            row['name']=row['name'].replace('SDSSJ','SDSS J')
            row['name']=row['name'].replace('-','–')
        elif 'WISEA0' in row['name']:
            print("WISEA detected", row['name'])
            row['name']=row['name'].replace('WISEA','WISEA J')
            row['name']=row['name'].replace('-','–')
        elif 'WISEAJ' in row['name']:
            print("WISEAJ detected", row['name'])
            row['name']=row['name'].replace('WISEAJ','WISEA J')
            row['name']=row['name'].replace('-','–')
        elif 'PSRJ' in row['name']:
            print("PSRJ detected", row['name'])
            row['name']=row['name'].replace('PSRJ','PSR J')
            row['name']=row['name'].replace('-','–')
        elif 'LPSMJ' in row['name']:
            print("LPSMJ detected", row['name'])
            row['name']=row['name'].replace('LPSMJ','LSPM J')
            row['name']=row['name'].replace('-','–')
        elif 'ULASJ' in row['name']:
            print("ULASJ detected", row['name'])
            row['name']=row['name'].replace('ULASJ','ULAS J')
            row['name']=row['name'].replace('-','–')
        elif 'LEHPM' in row['name']:
            print("LEHPM detected", row['name'])
            row['name']=row['name'].replace('LEHPM','LEHPM ')
        elif 'WDJ' in row['name']:
            row['name']=row['name'].replace('WDJ', 'WD J')
        elif 'LP' in row['name']:
            row['name']=row['name'].replace('LP', 'LP ')
        elif 'ULASJ' in row['name']:
            row['name']=row['name'].replace('ULASJ', 'ULAS J')
        row['name']=row['name'].replace('-','–')
        #if name=='.':
            #name=object_name_base+small_ra+small_dec
            #print('new name:', name)
        #elif name=='none':
            #name=object_name_base+small_ra+small_dec
            #print('new name:', name)
        #elif name == '0':
            #name=object_name_base+small_ra+small_dec
            #print('new name:', name)
        #print(thing)
    return


def format_coords(this_table):
    
    table_coordinates= coord.SkyCoord(ra=input_table['ra'], dec=input_table['dec'], unit = (u.deg, u.deg), frame='icrs')
    this_table['ra']=this_table['ra'].astype('S32')
    this_table['dec']=this_table['dec'].astype('S32')
    for row,coords in zip(this_table,table_coordinates):
        #coords= coord.SkyCoord(row['ra'], row['dec'], unit=(u.deg, u.deg))
        string_coords= coords.to_string(style='hmsdms')
        replace_chars= ['d','h','m']
        for char in replace_chars:
            string_coords= string_coords.replace(char, ':')
        string_coords=string_coords.replace('s', '')
        split_string=string_coords.split(' ')
        ra= split_string[0][:11] #limiting precision of the decimal
        dec= split_string[1][:12] #limiting precision of the decimal, need the additional index for sign
        row['ra']=ra
        row['dec']=dec
    return


def fix_sp_types(this_table):
    this_table['sp_type']=this_table['sp_type'].astype('S7')

    for row in this_table:
        if 'WDdM' in row['sp_type']:
            row['sp_type']='WD+dM'
        elif '??' in row['sp_type']:
            row['sp_type']='Unknown'
        else:
            pass
        
    
    return

def limit_sig_figs(this_table, mag_decimals=1, g_rp_decimals=2,pwd_decimals=2):
    #sig_fig_cols=['phot_g_mean_mag','g_rp','g_abs','pwd']
    sig_fig_dict={
        'phot_g_mean_mag':1,
        'teff':0,
        'teff_err':0,
        'logg':1,
        'logg_err':1,
        'm_wd':2,
        'm_wd_err':2,
        'c/he':1,
        'c/he_err':1
        
        }
    #for row in this_table:
        #row['phot_g_mean_mag']=np.round(row['phot_g_mean_mag'],decimals=mag_decimals)
    for row in this_table:
        for element in sig_fig_dict:
            row[element]=np.round(row[element],decimals=sig_fig_dict[element])
        #row['g_abs']=np.round(row['g_abs'],decimals=mag_decimals)
        #row['g_rp']=np.round(row['g_rp'],decimals=g_rp_decimals)
        #row['pwd']=np.round(row['pwd'],decimals=pwd_decimals)
    
    #this_table['pwd']=this_table['pwd'].astype('S4')
    #this_table['g_rp']=this_table['g_rp'].astype('S4')
    #for row in this_table:
    for column in this_table.colnames:
        this_table[column]=this_table[column].astype('S64')
        #if len(row['pwd'])==3:
            #row['pwd']=row['pwd']+'0'
        #if len(row['g_rp'])==3:
            #row['g_rp']=row['g_rp']+'0'
    #digits=sig_fig_dict.copy()
    for row in this_table:
        for colname in sig_fig_dict:
            try:
                print('\n',row[colname])
                row[colname]=row[colname]+("0"*(sig_fig_dict[colname]-len(row[colname].split('.')[1])))
                print(row[colname])
            except IndexError as error:
                print('IndexError:', error)
            #if len(row.split('.')[1])<sig_fig_dict[colname]:
                #row=row+("0"*sig_fig_dict[colname]-len(row.split('.')[1]))
            #else:
                #pass
    new_col_list=[]
    new_table=Table()
    newnewtable=Table()
    for name in this_table.colnames:
        if '_err' in name:
            print('\n',name)
            merged_col = [f"${v} \\pm {e}$" for v, e in zip(this_table[name.replace('_err','')], this_table[name])]
            # Add the merged column to the table
            print(name)
            new_table[name.replace('_err','')] = Column(merged_col)
    new_table.pprint()
    newnewtable=this_table[noerr_cols]
    newnewtable.pprint()
    for name in noerr_cols:
        try:
            newnewtable[name]=new_table[name]
        except KeyError:
            pass
    newnewtable.pprint()
    
    
    for name in colname_replacements:
        newnewtable.rename_column(name,colname_replacements[name])
        
    for name in newnewtable.colnames:
        if '/' in name:
            parts=name.split('/')
            parts[0]=parts[0].capitalize()
            parts[1]=parts[1].capitalize()
            newname='$\log(\\text{'+parts[0]+'}/\\text{'+parts[1]+'})$'
            newnewtable.rename_column(name,newname)
        else:
            pass
        
    return newnewtable


input_table.pprint()
subinds=np.where(input_table['sp_type']==desired_sptype)
subtable=input_table[subinds]
subtable=subtable[selected_output_cols]
#make_name(input_table)
#format_coords(input_table)
#fix_sp_types(input_table)
new_table=limit_sig_figs(subtable)

#input_table.pprint()
#output_table=Table(subtable[selected_output_cols])
output_table=new_table
print('output_table')
output_table.pprint()

output_name_base='DQpec_table'

final_name=output_name_base+'_'+spt.time_string()+'.tex'

output_dir='/Users/BenKaiser/Desktop/MORDOR_Survey_paper/tables/'

output_table.write(output_dir+final_name,format='ascii.latex')


#with open(output_dir+output_table,'r') as csvfile:
    #csvreader=csv.reader(csvfile)
    
    #skip the first row
    
    
    









