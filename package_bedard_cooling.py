"""
Created by Ben Kaiser (UNC-Chapel Hill) 2021-05-17

Take the white dwarf cooling model files from Bedard et al. 2020 (https://www.astro.umontreal.ca/
~bergeron/CoolingModels/), and compile them into single files that include all of the relevant information to 
be able to interpolate between them.

The files are structured in a strange way by default in which there are essentially lines with missing columns, 
namely the Mod column. Not to mention that it includes a #, so that will be shooting us in the foot when we 
try to read that in most likely. I'm not sure how to handle that one... I hope it doesn't comment out the whole 
header portion but I bet it will now...



"""

from __future__ import print_function
import numpy as np
import csv
#from astroquery.gaia import Gaia
#import astropy.units as u
#import astropy.coordinates as coord
from astropy.table import Table, QTable
#import matplotlib.pyplot as plt
#import scipy.stats as scistats
#import seaborn as sns
#import astropy

from glob import glob



input_file_string='*thin*'

input_file_list=sorted(glob(input_file_string))

input_file=input_file_list[0] #we're going to run the testing initially just using the first relevant file so we make sure we're at least looking in the right places and not screwing up formats... hopefully

cleaned_length=5



with open(input_file, 'r') as csvfile:
    reader=csv.reader(csvfile, delimiter=' ')
    index=0
    new_compiled_list=[]
    for row in reader:
        print('\n\n===',index,row,'\n\n====')
        copy_row=row[:]
        index_range=list(range(0,len(row)))
        index_range.sort(reverse=True)
        print(index_range, len(row))
        #print(index_range_reverse)
        #for number,entry in zip(index_range,row):
        for number in index_range:
            entry=row[number]
            
            #print(entry)
            if entry=='':
                #print('empty entry')
                del(row[number])
            elif entry==' ':
                pass
                #print('single space entry')
            else:
                pass
            #print(copy_row)
        if (len(row) > cleaned_length):
            print(len(row), 'is greater than', cleaned_length)
            print('so the first entry', row[0],"is presumably the row number in the original file...so we're going to delete it now")
            del(row[0])
        else:
            pass
        print('\n\n===',index,row,'\n\n====')
       
        #try:
                #if int(entry)==index:
                    #print(entry, 'appears to be', index, "so it's a row number and we don't want those... or it's a very unlucky coincidence... yikes!")
                    #del(row[number])
            #except ValueError:
                #pass
        index+=1
        new_compiled_list.append(row)
        
print(new_compiled_list)

#Ok now we've put together a full set of these things. We now need to wrap them around to make single rows for each actual row.

stacked_rows=3
stack_indices=range(0,stacked_rows)
for row in new_compiled_list:
    count=0
    if '==' in row[0]:
        print('equals signs row')
        pass
    else:
        new_row=row[:]
        if (count//stacked_rows==0):
            print('new row wrap')
        count+=1
        
    test_list= row
print(new_compiled_list[1])
print(new_compiled_list)
