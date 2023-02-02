"""
Created 2023-02-02 by Ben Kaiser (UNC-Chapel Hill).

@author: Ben Kaiser

This should take all the spectra in whatever directory and attempt to match them back to a master list of the objects in the survey that have spectra. It should then output a new list file with the names of each object, the coordinates, and the names of the spectral files that correspond to that object.

I can probably take the code from the cross matching into the Gentile Fusillo catalogue in gaia/ and then loosen the radius for searching because the RA and DEC in the Goodman headers does not correspond to the center of our ROI

***I should definitely load up a spectrum's header and search for that in the Gaia list because some objects won't have any spectra.

This raises another issue though, because the final list should almost certainly be sorted by RA




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


#plt.rc('font', size =18)

#print start
#import wdatmos
import spec_plot_tools as spt
import cal_params as cp


input_file= '../20190516B_retargeted_purple_search_gaia_scbd_20230131_update.csv'
input_table = Table.read(input_file)

current_directory=os.getcwd()
output_name_base=input_file.split('/')[1].split('.')[0]+'_'+current_directory.split('/')[-1]

search_radius = 5*60.#in arcseconds

spec_filenames=glob('*fits')
spec_filenames=sorted(spec_filenames)

new_cols=['400m1','400m2']

object_name_base='MORDOR J'

############################


ra_array = input_table['ra']
dec_array = input_table['dec']
name_array = input_table['name']
num_targets= ra_array.shape[0]

for new_col in new_cols:
    try:
        print('Test repeat column: ', input_table[new_col][0])
    except KeyError:
        print(new_col, " column doesn't exist, so we're making it.")
        new_length = len(input_table['dist'])
        if new_col=='400m1':
            false_array = np.full((new_length,), '', dtype=np.dtype('S64'))
            new_column = Column(false_array, name=new_col)
        elif new_col=='400m2':
            false_array = np.full((new_length,), '', dtype=np.dtype('S64'))
            new_column = Column(false_array, name=new_col)
        else:
            print("no matching names:", new_col)
        input_table.add_column(new_column)
        



input_table.pprint()
table_coordinates= coord.SkyCoord(ra=input_table['ra'], dec=input_table['dec'], unit = (u.deg, u.deg), frame='icrs')

print(table_coordinates)

def make_search_coords(ra,dec):
    """
    inputs:
    ra as string in decimal degrees or hour-angle
    dec as string in decimal degrees or degrees:arcmin:arcsec
    """    
    if ((":" in ra) or ( " " in ra) or ( "\t" in ra)):
        print("RA in hour angle")
        coordinate = coord.SkyCoord(ra = ra, dec =dec, unit = (u.hourangle, u.deg), frame = 'icrs')
    elif ('.' in ra):
        coordinate = coord.SkyCoord(ra = ra, dec =dec, unit = (u.deg, u.deg), frame = 'icrs')
    else:
        print("invalid RA present in dataset")
    #coordinate= coord.SkyCoord(ra=ra, dec=dec, unit = (u.hourangle, u.deg), frame='icrs')
    print('coordinate',coordinate)
    return coordinate


def check_coords(search_coord):
    inside_inds= np.where(table_coordinates.separation(search_coord) < (search_radius *u.arcsec))
    return inside_inds
    
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

table_row_inds=[]
for filename in spec_filenames:
    #output_table.pprint()
    header=fits.getheader(filename)
    search_coord=make_search_coords(header['ra'], header['dec'])
    inside_inds=check_coords(search_coord)
    print(filename,'inside_inds:',inside_inds)
    if '400m1' in filename:
        print('400m1 detected')
        input_table['400m1'][inside_inds]=filename
        #output_table.pprint()

    elif '400m2' in filename:
        print('400m2 detected')
        input_table['400m2'][inside_inds]=filename
        #output_table.pprint()

    else:
        print('No match found for coordinates:', header['ra'],header['dec'],filename)
    if inside_inds[0][0] in table_row_inds:
        pass
    else:
        table_row_inds.append(inside_inds[0][0])

print('\n\n',tuple(table_row_inds))

print(input_table[table_row_inds])

final_name=output_name_base+'_'+spt.time_string()+'.csv'
print('final_name:', final_name)
output_table=input_table[table_row_inds]
output_table.pprint()

#sort_order=sorted(output_table['ra'])
#print(sorted_order)
#output_table=output_table[sort_order]
output_table=output_table.group_by('ra')
print(output_table['name'])
make_name(output_table)
print(output_table['name'])
output_table.write(final_name,format='ascii.csv')


















