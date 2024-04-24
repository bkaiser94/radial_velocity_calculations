"""
Created by Ben Kaiser (UNC-Chapel Hill) 2024-04-24.

This should handle ZTF forced-photometry light curves (or make them I suppose) from the originally output ASCII table.



"""

from __future__ import print_function


#import matplotlib
#matplotlib.use('pdf')


import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const
from astropy import convolution as conv
from astropy.table import Table, Column
import scipy.interpolate as scinterp

import astropy.timeseries as ts

import time
start = time.time()

plt.rc('lines',markersize=4)

input_file='SDSSJ1312m0229_photometry_fromstart_thru_20240423.txt'


input_table=Table.read(input_file,format='ascii.basic') #using ascii.basic requires you to remove the commas from the header row.

input_table.pprint()
print(input_table.__doc__)
print(input_table.colnames)

field_array=np.unique(input_table['field'])
print('fields', field_array)
for field in field_array:
    field_inds=np.where(input_table['field']==field)
    print('field:', field, 'num obs:', len(input_table[field_inds]))
    cadence=np.diff(input_table['jd'][field_inds])
    print('min cadence for field '+str(field)+ ':', np.nanmin(cadence), ' days')
    cadence_hours=(cadence*u.day).to(u.hour)
    print('min cadence for field '+str(field)+ ':', np.nanmin(cadence_hours), ' hours')
    cadence_minutes=cadence_hours.to(u.minute)
    print('min cadence for field '+str(field)+ ':', np.nanmin(cadence_minutes), ' minutes',np.nanmax(cadence_minutes))
    print('time spanned', np.nanmax(input_table['jd'])-np.nanmin(input_table['jd']), 'days')
    
#plt.hist(cadence_minutes.value,bins=np.arange(0,223614,3.))
#plt.xlabel('cadence in minutes')
#plt.title('frequency of cadences')
#plt.show()

def get_cadence(table):
    cadence=np.diff(table['jd'])
    table['jd'].pprint()
    cadence_minutes=(cadence*u.day).to(u.minute)
    print('min cadence', np.min(cadence_minutes))
    print('max cadence', np.max(cadence_minutes))
    print('len(cadence)', len(cadence))
    return
def make_lightcurve(table, band='r',field=424,marker='o'):
    print('making light curve using band '+band+' and field ' +str(field))
    if len(band)==1:
        band='ZTF_'+band
        print('band name updated to ' + band)
    else:
        pass
    sub_inds=np.where(table['field']==field)
    sub_table=table[sub_inds].copy()
    subsub_inds=np.where(sub_table['filter'])
    subsub_table=sub_table[subsub_inds]
    get_cadence(subsub_table)
    print('len(subsub_table)', len(subsub_table))
    def fix_nulls(colname):
        null_inds=np.where(subsub_table[colname]=='null')
        subsub_table[colname][null_inds]=np.nan
        subsub_table[colname]=np.float_(subsub_table[colname])
        return
    fix_nulls('forcediffimflux')
    fix_nulls('forcediffimfluxunc')
    #null_inds=np.where(subsub_table['forcediffimflux']=='null')
    #subsub_table['forcediffimflux'][null_inds]=np.nan
    #subsub_table['forcediffimflux']=np.float_(subsub_table['forcediffimflux'])
    if band[-1]!='i':
        color=band[-1]
    else:
        color='magenta'
    plt.errorbar(subsub_table['jd'],subsub_table['forcediffimflux'],yerr=subsub_table['forcediffimfluxunc'],linestyle='none',color=color,label=band,marker=marker,alpha=0.4)
    print(np.nanmin(subsub_table['jd']), np.nanmax(subsub_table['jd']))
    print(np.nanmin(subsub_table['forcediffimflux']), np.nanmax(subsub_table['forcediffimflux']))
    #plt.scatter(subsub_table['jd'],np.float_(subsub_table['forcediffimflux']),label=band)

    
    return subsub_table

plt.show()
r_table=make_lightcurve(input_table,band='r',marker='o')
g_table=make_lightcurve(input_table,band='g',marker='s')
i_table=make_lightcurve(input_table,band='i',marker='^')
plt.xlabel('JD')
plt.ylabel('Flux difference')
plt.legend()
plt.show()


r_frequency, r_power = ts.LombScargle((r_table['jd']*u.day).to(u.hour),r_table['forcediffimflux'],r_table['forcediffimfluxunc']).autopower(minimum_frequency=1/(0.1*u.hour), maximum_frequency=1/(60.*u.hour))
#r_frequency, r_power = ts.LombScargle((r_table*u.day).to(u.hour),r_table['forcediffimflux']).autopower()
print(1/r_frequency)
print(r_power)
plt.plot(r_frequency,r_power)
plt.xlabel('frequency (1/hours)')
plt.show()

plt.plot(1./r_frequency,r_power)
plt.plot('period (hours)')
plt.show()
    
    
