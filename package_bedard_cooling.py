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
from astropy.table import vstack as tablevstack
#import matplotlib.pyplot as plt
#import scipy.stats as scistats
#import seaborn as sns
#import astropy
import os

from glob import glob



#input_file_string='*thin*'
input_file_string='*thin*'
original_dir= os.getcwd()
input_file_list=sorted(glob(input_file_string))
output_dir='/Users/BenKaiser/Desktop/Bedard_WD_cooling_models/cleaned_individual_masses/'

super_output_dir='/Users/BenKaiser/Desktop/Bedard_WD_cooling_models/'

#input_file=input_file_list[0] #we're going to run the testing initially just using the first relevant file so we make sure we're at least looking in the right places and not screwing up formats... hopefully

cleaned_length=5
mass_index=1
decimal_index=1

save_singles= True
#######################
def clean_single_file(input_file, save_singles=save_singles):

    mass_digits=input_file.split('_')[mass_index]

    mass_number= '.'.join([mass_digits[:decimal_index],mass_digits[decimal_index:]])

    print('mass_number', mass_number)

    with open(input_file, 'r') as csvfile:
        reader=csv.reader(csvfile, delimiter=' ')
        index=0
        new_compiled_list=[]
        for row in reader:
            #print('\n\n===',index,row,'\n\n====')
            copy_row=row[:]
            index_range=list(range(0,len(row)))
            index_range.sort(reverse=True)
            #print(index_range, len(row))
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

    num_stack=3
    stack_indices=range(0,num_stack)
    stacked_rows=[]
    wrapped_rows=[]
        
    count=0

    for row in new_compiled_list:
        if '==' in row[0]:
            print('equals signs row')
            pass
        else:
            new_row=row[:]
            print('count', count, 'num_stack', num_stack, count % num_stack)
            if (count % num_stack==0):
                print('new row wrap')
                print('previous fully wrapped row:', wrapped_rows)
                stacked_rows.append(wrapped_rows)
                wrapped_rows=[]
                wrapped_rows.extend(new_row)
            else:
                wrapped_rows.extend(new_row)
            count+=1
        #test_list= row

    print('\n\n=====\n\n')
    print(stacked_rows)
    print('\n\n=====\n\n')
    #print(new_compiled_list[1])
    print(new_compiled_list)

    ready_list=stacked_rows[2:]
    names_list=stacked_rows[1]
    array_form=np.array(ready_list)

    mass_array=np.full(array_form.shape[0], mass_number)
    #mass_array[:]=mass_number
    print('mass_array',mass_array)
    mass_array=np.float_(mass_array)
    print(mass_array)

    #print(array_form)
    #print(array_form.T)
    #print(array_form.shape)

    single_file_table=Table(array_form, names=names_list)
    single_file_table['mass']=mass_array


    print(single_file_table.info)
    single_file_table.pprint()

    #now I want to re-order the table to have mass first.... hopefully this works...

    old_columns=single_file_table.columns
    print('old_columns', old_columns)
    print('old_columns', old_columns[0])
    #print('single_file_table.names', single_file_table.col_names)
    #print(single_file_table[old_columns[:-1]])

    hold_columns=old_columns[:-1]

    new_table_list=[single_file_table['mass']]
    for this_col in hold_columns:
        new_table_list.extend([single_file_table[this_col]])

    print(new_table_list)
    new_single_table=Table(new_table_list)
    #new_single_table=Table([single_file_table['mass'],hold_columns])
    #new_single_table.pprint()
    #num_cols=len(old_columns)
    #for name_index in range(0,num_cols-1):
        #print('old_columns[name_index]',old_columns[name_index])
        #print('trying to add new col')
        #new_single_table[old_columns[name_index]]=single_file_table[old_columns[name_index]]
        #print('new col added')
    print('new_single_table')
    new_single_table.pprint()
    
    if save_singles:
        output_filename=input_file.split('.')[0]+'_cleaned'+'.csv'
        print('output_filename:', output_filename)
        #single
        os.chdir(output_dir)
        print('saving', output_filename, 'in', os.getcwd())
        new_single_table.write(output_filename, format='ascii.csv')
        print('saved', output_filename)
        os.chdir(original_dir)
    else:
        pass
    
    
    return new_single_table

##########################

def make_super_table(save_super_table=False):
    super_table=clean_single_file(input_file_list[0],save_singles=False)
    
    for input_file in input_file_list[1:]:
        single_table=clean_single_file(input_file, save_singles=False)
        super_table=tablevstack([super_table,single_table], join_type='exact')
    
    super_table.pprint()
    
    if save_super_table:
        single_input_file=input_file_list[0]
        first_strings=single_input_file.split('_')[:mass_index]
        last_strings=single_input_file.split('_')[mass_index+1:]
        print('first_strings', first_strings)
        print('last_strings', last_strings)
        #name_strings=first_strings.extend(last_strings)
        #print('name_strings', name_strings)
        new_name='_'.join([first_strings[0],last_strings[0]])
        output_filename='bedard2020_'+new_name.split('.')[0]+'_cleaned'+'.csv'
        print('output_filename:', output_filename)
        #single
        os.chdir(super_output_dir)
        print('saving', output_filename, 'in', os.getcwd())
        super_table.write(output_filename, format='ascii.csv')
        print('saved', output_filename)
        os.chdir(original_dir)
    else:
        pass
    return super_table




############

#for input_file in input_file_list:
    #clean_single_file(input_file)
    
make_super_table(save_super_table=False)

