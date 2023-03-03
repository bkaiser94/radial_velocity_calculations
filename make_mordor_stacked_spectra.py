"""
Created 2023-02-02 by Ben Kaiser (UNC-Chapel Hill)

@author: Ben Kaiser


This should make the stacked plots of the spectra for the MORDOR Survey chapters in my dissertation. I'll probably have to make the vertical spacing flexible, but the vertical extent of the plots will probably have to be rigidly set.

So with that in mind, the margins are supposed to be 1" on each side, so an 11" tall paper would have 9" to be occupied by the figure and caption... so I suppose that means I have 9"-caption height worth of inches to alot to the figure itself.

The caption for the J1644 spectrum from the Science paper is taking up 1.6" vertically in the dissertation.


"""
from __future__ import print_function


import matplotlib

matplotlib.use('pdf')
savefig=True


import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
from glob import glob
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const
from astropy import convolution as conv
from astropy.table import Table, Column
import scipy.interpolate as scinterp
import time
start = time.time()


#plt.rc('font', size =18)

#print start
#import wdatmos
import spec_plot_tools as spt
import cal_params as cp
import plot_spec as ps



#input_file='20190516B_retargeted_purple_search_gaia_scbd_20230131_update_DQpecs_1675365552.csv'
#input_file='20190516B_retargeted_purple_search_gaia_scbd_20230301_update_DQpecs_1677714150.csv'
input_file='20190516B_retargeted_purple_search_gaia_scbd_20230131_update_DZs_1675367585.csv'
#input_file='20190516B_retargeted_purple_search_gaia_scbd_20230131_update_WDdM_binaries_1675621755.csv'
#input_file='20190516B_retargeted_purple_search_gaia_scbd_20230131_update_DCs_1675625828.csv'
#input_file='20190516B_retargeted_purple_search_gaia_scbd_20230131_update_DCs_1676221722.csv'
#input_file='20190516B_retargeted_purple_search_gaia_scbd_20230301_update_DCs_1677713100.csv'
#input_file='20190516B_retargeted_purple_search_gaia_scbd_20230301_update_SDSS_objects_1677771324.csv'
#input_file='20190516B_retargeted_purple_search_gaia_scbd_20230131_update_polluted_WD_1677081843.csv' #WDdM binaries with polluted white dwarfs
input_table=Table.read(input_file)
current_directory=os.getcwd()
output_name_base=current_directory.split('/')[-1]


##DQpec 400M2
#smooth_size=5. #pixels
#offset_scale=1.
#max_spec_per_frame=6.
#text_y_offset=offset_scale/5.
##text_x_position=5400.
#text_x_position=8000.
#norm_range=[6630,6690]#wider double norm range
##norm_range=[7440, 7550]
#default_height=7.
#default_width=6.
#setup='400m2'
#line_id=''

#DZ 400M2 adjusted settings
smooth_size=3. #pixels
offset_scale=1.
max_spec_per_frame=6.
text_y_offset=-1*offset_scale/7.
#text_y_offset=0.
text_x_position=5210.
#norm_range=[6630,6690]#wider double norm range
norm_range=[7440, 7550]
default_height=7.
default_width=6.
setup='400m2'
line_id='cool_wd'


##WD+dM 400M2 adjusted settings
#smooth_size=3. #pixels
#offset_scale=1.
#max_spec_per_frame=6.
#text_y_offset=1*offset_scale/7.
##text_y_offset=0.
#text_x_position=6700.
##norm_range=[6630,6690]#wider double norm range
#norm_range=[7440, 7550]
#default_height=7.
#default_width=6.
#setup='400m2'
#line_id=''

##WD+dM 400M1 adjusted settings
#smooth_size=3. #pixels
#offset_scale=1.
#max_spec_per_frame=6.
#text_y_offset=3*offset_scale/7.
##text_y_offset=0.
#text_x_position=6000.
#norm_range=[6630,6690]#wider double norm range
##norm_range=[7440, 7550]
#default_height=7.
#default_width=6.*3350./4030.
#setup='400m1'
#line_id=''

###DC 400M2 adjusted settings
#smooth_size=5. #pixels
#offset_scale=1.
#max_spec_per_frame=6.
#text_y_offset=2*offset_scale/7.
##text_y_offset=0.
#text_x_position=5950.
#norm_range=[6630,6690]#wider double norm range
##norm_range=[7440, 7550]
#default_height=7.
#default_width=6.
#setup='400m2'
#line_id=''

####DC 400M1 adjusted settings
#smooth_size=5. #pixels
#offset_scale=1.
#max_spec_per_frame=6.
#text_y_offset=2*offset_scale/7.
##text_y_offset=0.
#text_x_position=6000.
#norm_range=[6630,6690]#wider double norm range
##norm_range=[7440, 7550]
#default_height=7.
#default_width=6.*3350./4030.
#setup='400m1'
#line_id=''


#DZ 400M1
smooth_size_400m1=smooth_size #pixels
offset_scale_400m1=offset_scale
text_y_offset_400m1=text_y_offset
#text_y_offset=0.
text_x_position_400m1=text_x_position
norm_range_400m1=norm_range#wider double norm range
#norm_range=[7440, 7550]
default_width_400m1=6.*3350./4030.





counter=0

spt.initiate_science_plot()
plt.rc('lines',linewidth=0.5)
tally=0
for row in input_table:
    if row['400m1'] != '':
        tally+=1


num_spec=[tally,max_spec_per_frame][np.greater(tally, max_spec_per_frame)]
print('num_spec',num_spec)

#plt.figure(figsize=(default_width_400m1,default_height*num_spec/max_spec_per_frame),constrained_layout=True)


#for row in input_table:
    #if row['400m1'] != '':
    ##row=input_table[index]
        #filename=row['400m1']
        #target_spec, header, target_noise= spt.retrieve_spec(filename)
        #hdu= fits.open(filename)
        #print('angstrom range',np.nanmax(target_spec[0])-np.nanmin(target_spec[0]))
        #ps.plot_spectrum(target_spec,'',header, smooth=True,norm=True,kernel_type='box', pix_width=smooth_size_400m1, offset=num_spec-1-(counter%num_spec)*offset_scale_400m1, norm_range=norm_range_400m1,color='k')
        #plt.text(text_x_position_400m1, num_spec-(counter%num_spec)*offset_scale_400m1+text_y_offset_400m1, row['name'])
        #counter+=1

#plt.title('')
#plt.ylabel(r'$\mathrm{F}_{\lambda}$ (Arbitrary Units)')
#plt.xlabel(r'Wavelength $(\mathrm{\AA})$')
#final_name=output_name_base+'_400m1_'+spt.time_string()+'.pdf'
##spt.show_plot(show_legend=False, actually_show=False)
#spt.show_plot(show_legend=False, line_id='cool_wd',show_label=False,actually_show=False,convert_to_air=True)
##plt.savefig(final_name)
#spt.show_plot(show_legend=False)



spt.initiate_science_plot()
plt.rc('lines',linewidth=0.5)
tally=0
#for row in input_table:
    #if row['400m2'] != '':
        #tally+=1
        
for row in input_table:
    if row[setup] != '':
        tally+=1
    #row=input_table[index]
print('total spectra that should be done', tally)    

num_frames=tally/int(max_spec_per_frame)

#list_spec_in_frame=[]
#for entry in range(0,num_frames):
    #list_spec_in_frame.append(max_spec_per_frame)
    
#list_spec_in_frame.append(tally%max_spec_per_frame)

#print('total number of spectra to be plotted and number alotted in frames summed', tally, np.sum(list_spec_in_frame))

num_spec=[tally,max_spec_per_frame][np.greater(tally, max_spec_per_frame)]
print('num_spec',num_spec)

plt.figure(figsize=(default_width,default_height*num_spec/max_spec_per_frame),constrained_layout=True)

current_frame=0
counter=0
for row in input_table:
    if row[setup] != '':
    #row=input_table[index]
        filename=row[setup]
        target_spec, header, target_noise= spt.retrieve_spec(filename)
        hdu= fits.open(filename)
        #print('angstrom range',np.nanmax(target_spec[0])-np.nanmin(target_spec[0]))
        ps.plot_spectrum(target_spec,'',header, smooth=True,norm=True,kernel_type='box', pix_width=smooth_size, offset=num_spec-1-(counter%num_spec)*offset_scale, norm_range=norm_range,color='k')
        plt.text(text_x_position, num_spec-(counter%num_spec)*offset_scale+text_y_offset, row['name'])
        counter+=1
        tally-=1
        if counter/num_spec==1:
            counter=0
            plt.title('')
            plt.ylabel(r'$\mathrm{F}_{\lambda}$ (Arbitrary Units)')
            plt.xlabel(r'Wavelength $(\mathrm{\AA})$')
            final_name=output_name_base+'_'+ setup+'_'+'frame'+str(current_frame)+'_'+spt.time_string()+'.pdf'
            #spt.show_plot(show_legend=False, actually_show=False)
            print('show plot right after this')
            #spt.show_plot(show_legend=False, line_id='cool_wd',show_label=False,actually_show=False,convert_to_air=True)
            #print('show plot just happened successfully')
            spt.show_plot(show_legend=False, line_id=line_id,show_label=False,actually_show=False,convert_to_air=True)
            print('show plot just happened successfully')
            if tally==0:
                #plt.ylim(-2,3) #DC 400M1 last plot. 
                #plt.ylim(-1,4) #DC 400M1 last plot. with J2320 correctly included.
                #plt.ylim(0,4.5) #DC 400M2 last plot. now that the 2 unknowns have been removed
                pass
            plt.savefig(final_name)
            plt.show()
            
            num_spec=[tally,max_spec_per_frame][np.greater(tally, max_spec_per_frame)]
            print('num_spec',num_spec)
            plt.figure(figsize=(default_width,default_height*num_spec/max_spec_per_frame),constrained_layout=True)
            current_frame+=1
        else:
            pass

plt.title('')
plt.ylabel(r'$\mathrm{F}_{\lambda}$ (Arbitrary Units)')
plt.xlabel(r'Wavelength $(\mathrm{\AA})$')
print('y lim set')
plt.ylim(-2,3)
print('current frame',current_frame)
final_name=output_name_base+'_'+setup+'_'+'frame'+str(current_frame)+'_'+spt.time_string()+'.pdf'
#spt.show_plot(show_legend=False, actually_show=False)
spt.show_plot(show_legend=False, line_id='cool_wd',show_label=False,actually_show=False,convert_to_air=True)
plt.savefig(final_name)


