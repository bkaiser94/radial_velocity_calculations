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


input_file='Lin2020_GALAH_thindisc.txt'

input_file=cp.abundance_dir+input_file

output_file=input_file.split('.')[0]+'_clean.txt'

all_array=np.genfromtxt(input_file, dtype='str', delimiter='\t')

print(all_array.shape)
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
print(all_array)
all_array=all_array

#alltable=Table.read(input_file)
#alltable.pprint()

alltable=Table(all_array[1:,:], names=all_array[0])
metal_bins=np.copy(all_array[0][1:])
alltable.pprint()
print(alltable.colnames)
alltable.add_index('Element')
print(alltable.loc['Na'])
print(metal_bins)

alltable.write(output_file, format='ascii')
