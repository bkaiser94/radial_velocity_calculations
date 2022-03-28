"""
Created by Ben Kaiser (UNC-Chapel Hill) 2022-03-18

Read in the copy-pasted text from Moss et al. 2022 (the WD age-metallicity paper)


"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
from astropy.time import Time
from astropy.table import Table, QTable
from astropy.table import vstack as tablevstack

import os
import csv


import spec_plot_tools as spt
#from plot_spec import plot_telluric_spectrum, plot_spectrum

input_filename='Moss_2022_WDMS_age_metallicity_bad.txt'
output_filename='Moss_2022_WDMS_age_metallicity.csv'

#raw_table=np.genfromtxt(input_filename)

#print(raw_table)

names_list=['Gaia_eDR3_ID', 'parallax','parallax_err','ps_g','ps_g_err','ps_r','ps_r_err','ps_i','ps_i_err','ps_z','ps_z_err','ps_y','ps_y_err','fe/h','fe/h_err','total_age','total_age_err']

collected_rows=[]

with open(input_filename,'r') as csvfile:
    reader=csv.reader(csvfile, delimiter=' ')
    index=0
    for row in reader:
        if index==0:
            print(row)
            #individual_vals=row.split(' ')
            new_row=[]
            for val in row:
                individual_val_int=int(val)
                new_row.append(individual_val_int)
            collected_rows.append(new_row)
        else:
            #individual_vals=row.split(' ')
            sub_index=0
            val_row=[]
            err_row=[]
            for val in row:
                if sub_index%2==0:
                    print('even',val)
                    val_row.append(float(val))
                elif sub_index%2 ==1:
                    print('odd',val)
                    clean_val=val.replace('(','')
                    clean_val=clean_val.replace(')','')
                    err_row.append(float(clean_val))
                else:
                    print('somehow not even nor odd... bad', val)
                sub_index+=1
            collected_rows.append(val_row)
            collected_rows.append(err_row)
        index+=1
            
print(collected_rows)

print('\n\n')
for row in collected_rows:
    print(len(row))
    print(row)

print(len(collected_rows),len(names_list))
full_table=Table(collected_rows,names=names_list)
full_table.pprint()

full_table.write(output_filename,format='csv')
    
#full_array=np.array(collected_rows)
#print('\n\n',full_array)
