"""
Created by Ben Kaiser (UNC-Chapel Hill) 2021-08-11

This is intended to be run in Python 3. Allegedly the "filter" function for lists behaves 
differently in the two pythons, so running in python 2 could mess up the different levels of the 
indexing. And realistically writing for python 2 now is an act of self-hatred pretty much... 
*shrug*. 

Take the diffusion timescale tables from http://www1.astrophysik.uni-kiel.de/~koester/
astrophysics/astrophysics.html that were sort of described in Koester's 2020 paper (https://
ui.adsabs.harvard.edu/abs/2020A&A...635A.103K/abstract), and reformat those files to nicely 
import as astropy tables in other scripts that actually use them.

This is a duplicated "package_bedard_cooling.py" so there could be artifact lines and 
comments that don't apply to this specific application. That might have been a mistake in 
 hindsight...



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
import sys

from glob import glob



#input_file_string='DA_diff_ov0.0.txt'
#input_file_string='DB_diff_ov0.0.txt'
#input_file_string='DA_diff_ov1.0.txt'
input_file_string='DB_diff_ov1.0.txt'

#input_file_string='*thin*'
original_dir= os.getcwd()
#input_file_list=sorted(glob(input_file_string))
output_dir='/Users/BenKaiser/Desktop/Koester_2020_diffusion_timescales/cleaned_files/'
 
super_output_dir='/Users/BenKaiser/Desktop/Koester_2020_diffusion_timescales/'

#input_file=input_file_list[0] #we're going to run the testing initially just using the first relevant file so we make sure we're at least looking in the right places and not screwing up formats... hopefully


#########################

if input_file_string[1]=='A':
    el_num_list=list(range(2,31))
elif input_file_string[1]=='B':
    el_num_list=list(range(3,31))
else:
    pass
print('el_num_list:',el_num_list)
col_names=['logg', 'teff','qcvz']
col_names.extend(el_num_list)
print('col_names:',col_names)

front_end_chaff= 30

block_skip_num=2
logg_row_num=2
diff_times=[]
diff_velocities=[]
time_chunk=True

#row types to choose between:
#garbage
#teff
#data
#diffusion_name
#logg

##########################


def pack_tau_group():
    
    
    return



############################


with open(input_file_string, 'r') as csvfile:
    reader=csv.reader(csvfile, delimiter=' ')
    index=0
    block_index=0
    lists_added=0
    new_compiled_list=[]
    new_compiled_array=np.array([[]])
    row_type='garbage'
    #el_num_list=[]
    chunk_list=[]
    for row in reader:
        print('++++++')
        print('\n\nindex', index,'\n\n')
        #print(row)
        if index < front_end_chaff:
            #print('skipping because front_end_chaff')
            row_type='garbage'
            #print(row)
            print('++++++')
        elif block_index < (block_skip_num):
            print('garbage row before table data')
            row_type='garbage'
            print(row)
        #elif block_index==logg_row_num:
            #print(index, block_index, row)
            #row_type='logg'
            #equal_index=row.index('=') #index in row of the equals sign (the logg value should be right after it
            #block_logg=row[equal_index+1] #should be the logg value of the overall chunk of data points. Will need to append it into rows at some point.
            #print('Found log(g) (I think...):',block_logg)
        #elif index > 200:
            #print('exiting to keep from running whole file')
            #sys.exit()
        else:
            print('data row cleaning')
            row_type='data'
            #print(row)
            row=list(filter(None, row)) #remove the empty entries retaining only the important numerical entries
            #print(row)
            #try:
                #row=row.remove('')
                #print(row)
            #except ValueError as error:
                #print(error)
            
        index += 1
        block_index +=1
        if row==[]:
            print('empty list', row)
            row_type='garbage'
        elif '-----' in row[0]:
            print('dashes row it seems')
            #print(row)
            #row_type='garbage'
            row_type='hyphen'
        elif '===' in row[0]:
            print('equals row it seems')
            #print(row)
            row_type='garbage'
        elif row[0]=='EL':
            #print(index, block_index, row)
            row_type='logg'
            print('logg row')
            
        elif 'Z' in row[0]:
            print('teff row')
            row_type='teff'
            
            
        elif row[0]=='Diffusion':
            print(row[0],'indicating beginning of table chunk')
            row_type='diffusion_name'
            print('second entry in row', row[1])
            if row[1]=='time':
                print('should be diffusion timescales in this part')
                time_chunk=True
            elif row[1]=='velocities':
                time_chunk=False
            block_index = 0
        elif row[0] == 'qcvz':
            row_type='qcvz'
            
        elif row_type=='garbage':
            pass
        else:
            print('should be actual data to include in the table')
        
        #index += 1
        #block_index +=1
        
        ####### Now the row_type should have been designated and we can proceed without embedding increasingly complicated stuff.
        if time_chunk==False:
            #This should exclude us in the event that we're in the diffusion velocities portion of the data tables, so it doesn't compile it into whatever is going on.
            print('not in a diffusion timescale chunk of the data tables')
        elif row_type=='data':
            print("row_type=='data'", row_type)
            print(row)
            #el_num_list.append([row[0]])
            add_row=row[1:]
            chunk_list.append(add_row)
            print('chunk_list appended', chunk_list)
        elif row_type=='garbage':
            print("row_type=='garbage'", row_type)
        elif row_type=='teff':
            print("row_type=='teff'", row_type)
            print(row)
            add_row=row[1:]
            chunk_list.append(add_row)
            print('chunk_list appended', chunk_list)
        elif row_type=='logg':
            print("row_type=='logg'", row_type)
            
            ####stuff to do before moving to the next thing
            if len(chunk_list)>0:
                print('\nlogg row and chunk_list is empty, so we need to output before moving on.\n')
                print('time_chunk:', time_chunk)
                row_length=len(chunk_list[0])
                print('row_length',row_length)
                logg_list=[]
                for thing in range(0,row_length):
                    logg_list.append(block_logg)
                #print('chunk_list', chunk_list)
                #print('logg_list', logg_list)
                logg_list=[logg_list]
                #print('logg_list', logg_list)
                logg_list.extend(chunk_list)
                #print('chunk_list', chunk_list)
                #print('logg_list', logg_list)
                chunk_list=logg_list
                #print('chunk_list', chunk_list)
                #print('chunk_list[0]',chunk_list[0])
                #print('\n\nchunk_list[1]',chunk_list[1])
                
                chunk_array=np.array(chunk_list)
                #print('chunk_array', chunk_array)
                print('chunk_array.shape',chunk_array.shape)
                #print('transpose:',chunk_array.T)
                #for thing in chunk_list:
                    #print(len(thing))
                print('chunk_list resetting')
                chunk_list=[]
                if lists_added==0:
                    new_compiled_array=np.copy(chunk_array.T)
                else:
                    new_compiled_array=np.append(new_compiled_array, np.copy(chunk_array.T),axis=0)
                lists_added+=1
                #print('trying to create an astropy table')
                #trial_table=Table(chunk_array.T,names=col_names)
                #trial_table.pprint()
            
            print(row)
            equal_index=row.index('=') #index in row of the equals sign (the logg value should be right after it
            block_logg=row[equal_index+1] #should be the logg value of the overall chunk of data points. Will need to append it into rows at some point.
            print('Found log(g) (I think...):',block_logg)
        elif row_type=='qcvz':
            print("row_type=='qcvz'", row_type)
            print(row)
            add_row=row[1:]
            chunk_list.append(add_row)
            print('chunk_list appended', chunk_list)
        #elif ((row_type=='hyphen') and (len(chunk_list)>0)):
            #print("\nhit a hyphen row in what should be a time chunk, so we're going to transpose and tack it on.\n")
            #print('time_chunk:', time_chunk)
            #row_length=len(chunk_list[0])
            #print('row_length',row_length)
            #logg_list=[]
            #for thing in range(0,row_length):
                #logg_list.append(block_logg)
            #print('chunk_list', chunk_list)
            #print('logg_list', logg_list)
            #logg_list=[logg_list]
            #print('logg_list', logg_list)
            #logg_list.extend(chunk_list)
            ##print('chunk_list', chunk_list)
            #print('logg_list', logg_list)
            #chunk_list=logg_list
            #print('chunk_list', chunk_list)
            #print('chunk_list[0]',chunk_list[0])
            #print('\n\nchunk_list[1]',chunk_list[1])
            
            #chunk_array=np.array(chunk_list)
            #print('chunk_array', chunk_array)
            #print('chunk_array.shape',chunk_array.shape)
            #print('transpose:',chunk_array.T)
            #for thing in chunk_list:
                #print(len(thing))
            #print('chunk_list resetting')
            #chunk_list=[]
        if row_type=='diffusion_name':
            print("row_type=='diffusion_name'", row_type)
            print(row)
            if time_chunk == False:
                #this means we just finished a time_chunk section because we have now encountered a diffusion timescale, so we need to add whatever actions are needed related to compiling the data.
                row_length=len(chunk_list[0])
                print('row_length',row_length)
                logg_list=[]
                for thing in range(0,row_length):
                    logg_list.append(block_logg)
                #print('chunk_list', chunk_list)
                #print('logg_list', logg_list)
                logg_list=[logg_list]
                #print('logg_list', logg_list)
                logg_list.extend(chunk_list)
                #print('chunk_list', chunk_list)
                #print('logg_list', logg_list)
                chunk_list=logg_list
                #print('chunk_list', chunk_list)
                #print('chunk_list[0]',chunk_list[0])
                #print('\n\nchunk_list[1]',chunk_list[1])
                
                chunk_array=np.array(chunk_list)
                #print('chunk_array', chunk_array)
                print('chunk_array.shape',chunk_array.shape)
                
                #for thing in chunk_list:
                    #print(len(thing))
                #print('transpose:',chunk_array.T)
                #for thing in chunk_list:
                    #print(len(thing))
                print('chunk_list resetting')
                chunk_list=[]
                new_compiled_array=np.append(new_compiled_array, np.copy(chunk_array.T),axis=0)
                #print('trying to create an astropy table')
                #trial_table=Table(new_compiled_array,names=col_names)
                #trial_table.pprint()
            elif time_chunk==True:
                pass
            else:
                pass
        else:
            #print('row_type not recognized...', row_type)
            #print(row)
            pass
        
        print('++++++')



print('trying to create an astropy table')
output_table=Table(new_compiled_array,names=col_names)
output_table.pprint()

core_name=''.join(input_file_string.split('.')[0:-1])
output_filename="Koester2020_"+core_name+'_diffusion_timescales'+'.csv'


os.chdir(output_dir)
print('saving', output_filename,'in',os.getcwd())
output_table.write(output_filename,format='ascii.csv')
print('saved', output_filename)























