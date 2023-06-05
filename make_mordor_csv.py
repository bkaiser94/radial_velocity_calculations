"""
Created by Ben Kaiser (UNC-Chapel Hill) 2023-04-10.

Make the CSV to be posted on Zenodo of all of the MORDOR Survey objects with relevant Gaia 
DR2 data. Names should also be fixed to be what would be desired. I don't think I'm going to 
clean up the file otherwise though because I really can't be bothered at this point unfortunately.


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

import spec_plot_tools as spt
import cal_params as cp

input_file='20190516B_retargeted_purple_search_gaia_scbd_20230301_update.csv'

input_table = Table.read(input_file)


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
            #row['name']=row['name'].replace('-','–')
        elif 'SDSSJ' in row['name']:
            print("SDSS detected", row['name'])
            row['name']=row['name'].replace('SDSSJ','SDSS J')
            #row['name']=row['name'].replace('-','–')
        elif 'WISEA0' in row['name']:
            print("WISEA detected", row['name'])
            row['name']=row['name'].replace('WISEA','WISEA J')
            #row['name']=row['name'].replace('-','–')
        elif 'WISEAJ' in row['name']:
            print("WISEAJ detected", row['name'])
            row['name']=row['name'].replace('WISEAJ','WISEA J')
            #row['name']=row['name'].replace('-','–')
        elif 'PSRJ' in row['name']:
            print("PSRJ detected", row['name'])
            row['name']=row['name'].replace('PSRJ','PSR J')
            #row['name']=row['name'].replace('-','–')
        elif 'LPSMJ' in row['name']:
            print("LPSMJ detected", row['name'])
            row['name']=row['name'].replace('LPSMJ','LSPM J')
            #row['name']=row['name'].replace('-','–')
        elif 'ULASJ' in row['name']:
            print("ULASJ detected", row['name'])
            row['name']=row['name'].replace('ULASJ','ULAS J')
            #row['name']=row['name'].replace('-','–')
        elif 'LEHPM' in row['name']:
            print("LEHPM detected", row['name'])
            row['name']=row['name'].replace('LEHPM','LEHPM ')
        elif 'WDJ' in row['name']:
            row['name']=row['name'].replace('WDJ', 'WD J')
        elif 'LP' in row['name']:
            row['name']=row['name'].replace('LP', 'LP ')
        elif 'ULASJ' in row['name']:
            row['name']=row['name'].replace('ULASJ', 'ULAS J')
        #row['name']=row['name'].replace('-','–')
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


#################3
##################



make_name(input_table)
output_table=input_table


output_name_base='full_MORDOR_survey'

final_name=output_name_base+'_'+spt.time_string()+'.csv'

output_table.write(final_name,format='ascii.csv')



