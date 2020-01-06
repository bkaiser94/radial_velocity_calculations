"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-12-19

Plot SED's with correctly indicated limiting error bars

I tried to do this in Jupyter notebooks, but it was refusing to plot the error bars for whatever reason in the correct way that I wanted them. I also couldn't really resize the windows which is pretty infuriating.

"""

from __future__ import print_function
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

import spec_plot_tools as spt
import cal_params as cp
#import plot_spec as ps
print(os.getcwd())

#target_phot_file=''
#target_phot_file=''

#J1644
wd_name='GaiaJ1644-0449'
target_phot_file='J1644_photometry.vot'

#J2356
#wd_name='WDJ2356-209'
#target_phot_file='J2356_photometry.vot'


##J1330
#wd_name='SDSSJ1330+6435'
#target_phot_file='J1330_photometry.vot'

def resave_vot_as_csv():
    target_table=Table.read(target_phot_file)
    new_name=target_phot_file.split('.')
    new_name=new_name[0]+'.csv'
    target_waves=target_table['sed_freq'].to(u.micron, equivalencies=u.spectral())
    target_table.add_column(Column(target_waves, name=('sed_wave')))
    target_table.write(new_name)
    return

#resave_vot_as_csv()

fill_val=5.e-5

target_table=Table.read(target_phot_file)
target_waves=target_table['sed_freq'].to(u.micron, equivalencies=u.spectral())
target_table.add_column(Column(target_waves, name=('sed_wave')))
#plt.figure(figsize=(15,15))
plt.scatter(target_table['_RAJ2000'], target_table['_DEJ2000'])
for row in target_table:
    plt.annotate(str(row['sed_filter']),xy=(row['_RAJ2000'], row['_DEJ2000']), xycoords='data', xytext=(row['_RAJ2000'], row['_DEJ2000']), textcoords= 'data' , fontsize=12)

#plt.xlim(202.498,202.508)
#plt.ylim(64.5885, 64.5905)
plt.show()

#plt.figure(figsize=(15,15))
measured=np.where(target_table['sed_eflux'].mask== False)
limits= np.where(target_table['sed_eflux'].mask!= False)
print(target_table['sed_eflux'][measured])
#for wave, row in zip(target_waves, target_table):
    ##plt.scatter(wave, row['sed_flux'], marker='o', color='b')
    #try:
        #print('trying')
        #if row['sed_eflux'].mask==True:
            #print('masked')
            ##print(type(row['sed_eflux']).filled())
            #plt.errorbar(wave.value, row['sed_flux'].quantity.value, yerr=row['sed_eflux'].filled(5e-5).quantity.value, marker='o',uplims=True)
            ##print(row['sed_eflux'].filled(5).quantity.value)
    #except AttributeError:
        #pass
    ##else:
        #plt.errorbar(wave.value, row['sed_flux'], yerr=row['sed_eflux'], linestyle='none', marker='o')


#    try:
#        plt.errorbar(wave, row['sed_flux'],yerr=row['sed_eflux'], linestyle='none', marker='o', color='b')
#    except ValueError:
#        plt.errorbar(wave, row['sed_flux'], linestyle='none', marker='o')
measured_table=target_table[measured]
limits_table=target_table[limits]
limits_table=limits_table.filled(limits_table['sed_flux']*0.5)
measured_table.pprint()
#plt.errorbar(target_waves[measured], target_table['sed_flux'][measured],yerr=target_table['sed_eflux'][measured], linestyle='none', marker='o')
#plt.errorbar(target_waves[limits], target_table['sed_flux'][limits], yerr=target_table['sed_eflux'][limits].filled(fill_val).quantity.value,linestyle='none', marker='o', color='r')
plt.errorbar(measured_table['sed_wave'], measured_table['sed_flux'],yerr=measured_table['sed_eflux'], linestyle='none', marker='o')
plt.errorbar(limits_table['sed_wave'], limits_table['sed_flux'],yerr=limits_table['sed_eflux'], linestyle='none', marker='o', color='r', uplims=True)
#plt.scatter(target_waves, target_table['sed_flux'], marker='s', color='b')
#plt.errorbar(target_waves.value, J1330_table['sed_flux'].quantity.value, yerr=J1330_table['sed_eflux'].filled(0).quantity.value, linestyle='none', marker='o',solid_capstyle='projecting')
plt.xscale('log')
plt.ylim(np.min(target_table['sed_flux'])*0.9, np.max(target_table['sed_flux'])*1.1)
plt.yscale('log')
plt.xlabel(r'Wavelength ($\mu$m)')
plt.ylabel(r'$f_{\nu}$ (Jy)')
plt.axvline(x=1.8,linestyle='--', color='k')
plt.title(wd_name)
plt.show()





