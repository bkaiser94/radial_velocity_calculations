"""
Created by Ben Kaiser (UNC-Chapel Hill) 2020-01-07

Basically this should take in some data file and then output it in a format that is far more friendly to some other file that could then use it.

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


#input_file='Lin2020_GALAH_thindisc.txt'
input_file='Lin2020_GALAH_thickdisc.txt'

input_file=cp.abundance_dir+input_file

output_filem=input_file.split('.')[0]+'_m'+'_clean.txt'
output_fileb=input_file.split('.')[0]+'_b'+'_clean.txt'

all_array=np.genfromtxt(input_file, dtype='str', delimiter='\t')
all_array=all_array.astype('<U20')
print(all_array.shape)
print(all_array.dtype)
print(all_array)

new_list=[]
#for thing in all_array[0]:
    #print(thing)
    #new_list.append(thing.decode)
    #thing=thing.replace(' ', '')
    #for other in thing:
        #print(other)
        #other=other.replace(' ', '')
        #print(other)

for i in range(0, all_array.shape[0]):
    print(all_array[i,0])
    all_array[i,0]=all_array[i,0].replace(' ', '')
    for j in range(0,all_array.shape[1]-1):
        all_array[i,j+1]=all_array[i,j+1].replace(' ', '')
        all_array[i,j+1]=all_array[i,j+1].replace('−', '-')
        if i==0:
            pass
        else:
            pass
            #all_array[i,j+1]='[' + all_array[i,j+1] +']'
            #all_array[i,j+1]=all_array[i,j+1].split(',')
print(all_array)
#all_array=all_array

#alltable=Table.read(input_file)
#alltable.pprint()

alltable=Table(all_array[1:,:], names=all_array[0])
copy_table=alltable.copy()
metal_bins=np.copy(all_array[0][1:])
alltable.pprint()
print(alltable.colnames)
alltable.add_index('Element')
print(alltable.loc['Na'])
print(metal_bins)

for row in alltable:
    for thing in metal_bins:
        #print(type(row[thing].split(',')))
        row[thing]=row[thing].split(',')[0]
        
for row in copy_table:
    for thing in metal_bins:
        #print(type(row[thing].split(',')))
        row[thing]=row[thing].split(',')[1]

#alltable[metal_bins[0]]=alltable[metal_bins[0]].astype(list)


#alltable.write(output_file, format='ascii', overwrite=True)
alltable.write(output_filem, format='ascii', overwrite=True)
copy_table.write(output_fileb, format='ascii', overwrite=True)

newtable=Table.read(output_filem, format='ascii')
newtable.pprint()
newtable.add_index('Element')
print(newtable['Element'])
#print(newtable.loc('Na')['[-0.5,-0.1]'])
print(newtable.loc['Na'])
print(type(newtable.loc['Na']['[-0.5,-0.1]']))
