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
import csv

import astropy.timeseries as ts

import time
start = time.time()

plt.rc('lines',markersize=4)

#input_file='SDSSJ1312m0229_photometry_fromstart_thru_20240423.txt'
#input_file='elBadry_P2p00_photometry_from_start_thru_20240512.txt'
#input_file='J1312_nearby_G12_photometry_from_start_thru_2040425.txt'
#input_file='J1312_nearby_G19_photometry_from_start_thru_20240424.txt'
input_file_core='*photometry_from_start_thru*'

command_arg=sys.argv[1]

input_file=command_arg+input_file_core
print(input_file)
print(type(input_file))
input_file=glob(input_file)[0]

print(input_file)
print(type(input_file))

snr_threshold=3. #value per the documentation
snr_uplim=5. #5-sigma upper limit as the value; also per the documentation
min_period=0.1 #minimum period to test in hours
max_period=24. #maximum period to test in hours.

crude_file=open(input_file,'r')
file_lines=crude_file.readlines()
ra_found=False
dec_found=False
ra=0. #is overwritten later
dec=0. #is overwritten later
alpha=0.4
field_num=484


def fix_band_name(band):
    if len(band)==1:
        band='ZTF_'+band
        print('band name updated to ' + band)
    else:
        pass
    
    return band

#get target coordinates from the forced photometry file
for line in file_lines:
    if (ra_found and dec_found):
        print('both RA and Dec found:', ra, ',', dec)
        break
    else:
        if 'Requested input R.A.' in line:
            ra_string=line.split('=')[1]
            ra_substring=ra_string.split(' ')[1]
            ra=float(ra_substring)
            ra_found=True
            #print('ra_substring', ra_substring)
        elif 'Requested input Dec.' in line:
            dec_string=line.split('=')[1]
            dec_substring=dec_string.split(' ')[1]
            print('dec_substring',dec_substring)
            dec=float(dec_substring)
            dec_found=True
crude_file.close()




input_table=Table.read(input_file,format='ascii.basic') #using ascii.basic requires you to remove the commas from the header row in the ZTF tables as output based on my attempts. Maybe there's a better format to choose, but this has worked for me and wasn't that hard to do.

#input_table.pprint()
#print(input_table.__doc__)
print(input_table.colnames)

print('first JD', input_table['jd'][0])

proc_inds=np.where(input_table['procstatus']==0)
input_table=input_table[proc_inds]
input_table['forcediffimflux']=Column(input_table['forcediffimflux'], dtype=np.float64) #fixing this to be a datatype that will play nicely with math later.
input_table['forcediffimfluxunc']=Column(input_table['forcediffimfluxunc'],dtype=np.float64)
print(input_table['procstatus'])
##### add a column of the bmjd_TDB times from the JD times in the file. ###
#obs_time = Time(input_table['jd'], format = 'jd', location =coords.EarthLocation.of_site( 'Palomar'))
obs_time = Time(input_table['jd'], format = 'jd', location =coords.EarthLocation.of_site( 'Palomar'))

target_coord = coords.SkyCoord(ra, dec, frame = 'icrs', unit= (u.hourangle, u.deg), )
bary_corr =obs_time.tdb.light_travel_time(target_coord)
bmjd_tdb_val = (obs_time.tdb+ bary_corr.tdb).mjd
print('BMJD_TDB', bmjd_tdb_val, type(bmjd_tdb_val))
input_table.add_column(bmjd_tdb_val, name='bmjd_tdb')

################

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
    
    
print('\n\nresetting target field_num',field_num, 'to be first field number in the lightcurve')
field_num=field_array[0]
print('New field_num:',field_num)
#plt.hist(cadence_minutes.value,bins=np.arange(0,223614,3.))
#plt.xlabel('cadence in minutes')
#plt.title('frequency of cadences')
#plt.show()

def get_cadence(table):
    cadence=np.diff(table['jd'])
    #table['jd'].pprint()
    cadence_minutes=(cadence*u.day).to(u.minute)
    min_cadence_ind=np.argmin(cadence)
    print('min cadence JD', table['jd'][min_cadence_ind])
    print('min cadence', np.min(cadence_minutes))
    print('max cadence', np.max(cadence_minutes))
    print('len(cadence)', len(cadence))
    return


def make_lightcurve(table, band='r',field=field_num,marker='o'):
    print('making light curve using band '+band+' and field ' +str(field))
    band=fix_band_name(band)
    sub_inds=np.where(table['field']==field)
    sub_table=table[sub_inds].copy()
    subsub_inds=np.where(sub_table['filter'])
    subsub_table=sub_table[subsub_inds]
    get_cadence(subsub_table)
    #print('len(subsub_table)', len(subsub_table))
    def fix_nulls(colname):
        null_inds=np.where(subsub_table[colname]=='null')
        print('null_inds', null_inds)
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

    #print(np.nanmin(subsub_table['jd']), np.nanmax(subsub_table['jd']))
    #print(np.nanmin(subsub_table['forcediffimflux']), np.nanmax(subsub_table['forcediffimflux']))
    #plt.scatter(subsub_table['jd'],np.float_(subsub_table['forcediffimflux']),label=band)

    
    return subsub_table


def get_abs_flux(band_table):
    """
    inputs:
        band_table: astropy table that has already been limited to a single photometric band... I suppose this isn't necessarily required in reality though because it goes row by row in its evaluation, so different bands will have different ref mags most likely. Whatever.
    This is all based on section 13 of the forced photometry documentation.
    """
    nearest_ref_flux=10.**(0.4*(band_table['zpdiff']-band_table['nearestrefmag']))
    nearest_ref_flux_unc=band_table['nearestrefmagunc']*nearest_ref_flux/1.0857 
    print('inside absflux', band_table['forcediffimflux'],nearest_ref_flux)
    print(type(band_table['forcediffimflux']),type(nearest_ref_flux))
    print(band_table['forcediffimflux'].dtype)
    print(nearest_ref_flux.dtype)
    #print(band_table['forcediffimflux'].astype(np.float64))
    #total_flux=band_table['forcediffimflux'].astype(np.float64)+nearest_ref_flux
    total_flux=band_table['forcediffimflux']+nearest_ref_flux
    total_flux_unc=np.sqrt(band_table['forcediffimfluxunc']**2.-nearest_ref_flux_unc**2.)
    total_snr=total_flux/total_flux_unc
    band_table.add_columns([total_flux, total_flux_unc, total_snr], names=['totalflux','totalfluxunc', 'totalfluxsnr'])
    return

def get_mags(table, band='r'):
    #good_snr_inds=np.where(table['forcediffimflux']/table['forcediffimfluxunc']>snr_threshold)
    #bad_snr_inds=np.where(table['forcediffimflux']/table['forcediffimfluxunc']<=snr_threshold)
    #band=fix_band_name(band)
    #mags=table['zpdiff']-2.5*np.log10(table['forcediffimflux'])
    mags=table['zpdiff']-2.5*np.log10(table['totalflux'])
    print(mags)
    return mags

def plot_table(table, marker='o', color='r', label='r'):
    
    
    
    return


get_abs_flux(input_table)
all_mags=get_mags(input_table)

plt.errorbar(input_table['bmjd_tdb'],all_mags,marker='o',linestyle='none')
plt.xlabel('BMJD_TDB')
plt.ylabel('All mags')
plt.title(input_file)
plt.show()

plt.errorbar(input_table['bmjd_tdb'], input_table['totalflux'], yerr=input_table['totalfluxunc'], linestyle='none', marker='o')
plt.xlabel('bmjd_tdb')
plt.ylabel('total flux')
plt.show()

r_table=make_lightcurve(input_table,band='r',marker='o')
g_table=make_lightcurve(input_table,band='g',marker='s')
i_table=make_lightcurve(input_table,band='i',marker='^')
print('r_table cadence')
get_cadence(r_table)
plt.xlabel('JD')
plt.ylabel('Flux difference')
plt.legend()
plt.show()

print('getting r mags')
r_mags=get_mags(r_table)
print('getting g mags')
g_mags=get_mags(g_table)
print(' getting i mags')
i_mags=get_mags(i_table)

#nan_inds=np.where(np.isnan(r_mags))
#r_mags[nan_inds]=0.
#print('nan_inds',nan_inds)
#print(r_mags[0])
#print(r_mags[0]==np.nan)
#print(type(r_mags[0]), r_mags[0]*10., np.nan*10.)

plt.errorbar(r_table['bmjd_tdb'],r_table['totalflux'],yerr=r_table['totalfluxunc'], linestyle='none',marker='o',color='r',label='r',alpha=alpha)
plt.errorbar(g_table['bmjd_tdb'],g_table['totalflux'],yerr=g_table['totalfluxunc'], linestyle='none',marker='s',color='g',label='g',alpha=alpha)
plt.errorbar(i_table['bmjd_tdb'],i_table['totalflux'],yerr=i_table['totalfluxunc'], linestyle='none',marker='^',color='magenta',label='i',alpha=alpha)
plt.xlabel('BMJD_TDB')
plt.ylabel('Absolute Flux')
plt.legend()
plt.show()


plt.errorbar(r_table['bmjd_tdb'],r_mags,color='r', label='r',linestyle='none', marker='o', alpha=0.4)
plt.errorbar(r_table['bmjd_tdb'],g_mags,color='g', label='g',linestyle='none', marker='s', alpha=0.4)
plt.errorbar(r_table['bmjd_tdb'],i_mags,color='magenta', label='i',linestyle='none',marker='^',alpha=0.4)
plt.legend()
plt.ylabel('Mag')
plt.xlabel('BMJD_TDB')
plt.show()

plt.errorbar(r_table['forcediffimflux'],r_mags,color='r', label='r',linestyle='none', marker='o', alpha=0.4)
plt.ylabel('Mag (0 is used for nans)')
plt.xlabel('forcediffimflux')
plt.show()


elapsed_hours_col=Column(((r_table['bmjd_tdb']-r_table['bmjd_tdb'][0])*u.day).to(u.hour), name='time_hours')
##now outputing
#print('saving file')
##output_table=Table([r_table['bmjd_tdb'],r_table['forcediffimflux'],r_table['forcediffimfluxunc']])
#output_table=Table([elapsed_hours_col,r_table['forcediffimflux'],r_table['forcediffimfluxunc']])
######output_table=Table([r_table['jd'],r_table['forcediffimflux'],r_table['forcediffimfluxunc']])
##output_table.write('SDSSJ1312m0229_rband_photometry_'+str(time.time()).split('.')[0]+'.txt',format='csv', delimiter='\t')
#output_table.write('SDSSJ1312m0229_gband_photometry_'+str(time.time()).split('.')[0]+'.txt',format='csv', delimiter='\t')

#print('file saved')

r_frequency, r_power = ts.LombScargle((r_table['bmjd_tdb']*u.day).to(u.hour),r_table['forcediffimflux'],r_table['forcediffimfluxunc']).autopower(minimum_frequency=1/(max_period*u.hour), maximum_frequency=1/(min_period*u.hour))
#r_frequency, r_power = ts.LombScargle((r_table*u.day).to(u.hour),r_table['forcediffimflux']).autopower()
print(1/r_frequency)
print(r_power)
top_args=np.argsort(r_power)
print('top powers', np.flip(r_power[top_args]))
top6_freqs=np.flip(r_frequency[top_args])[:6]
print('top 6 frequencies:', top6_freqs)
top6_periods=1./top6_freqs
print('top 6 periods:', top6_periods, 'hrs')

plt.plot(r_frequency,r_power)
#for x in top6_freqs:
    #plt.axvline(x=x, linestyle=':', color='r')

plt.xlabel('frequency (1/hours)')
plt.show()

plt.plot(1./r_frequency,r_power)
plt.plot('period (hours)')
plt.show()

for period in top6_periods:
    folded_times=(r_table['bmjd_tdb']*u.day).to(u.hour)%period
    plt.errorbar(folded_times.value, r_table['forcediffimflux'], yerr=r_table['forcediffimfluxunc'], linestyle='none',marker='o')
    plt.xlabel('Time (hours)')
    plt.title('Period: '+str(period) + ' ZTF r-band')
    plt.ylabel('Force Difference Image flux')
    plt.show()

total_time=12.*60.#number of minutes in a long-period ELM approximately

spacings=np.array([0.,1.,3.,5.,7.,11.,13.,17.])
phases=np.random.rand(spacings.shape[0])
sub_periods=np.random.rand(spacings.shape[0])*(total_time-30.)+30.
#sub_periods=total_time/spacings
all_times=np.linspace(0.,total_time,10000)
def make_sine_function(x_vals,period,amp=4.,phase=0.):
    interior=2*np.pi*x_vals/period+(phase*2.*np.pi)
    
    return amp*np.sin(interior)


for period,phase in zip(sub_periods,phases):
    plt.plot(all_times,make_sine_function(all_times,period,phase=phase), label='period: '+str(period/60.)[:4]+' hr, phase: '+str(phase)[:4])
    plt.plot(spacings*10.,make_sine_function(spacings*10.,period,phase=phase), linestyle='none',marker='o')
plt.legend()
for pos in spacings*10.:
    plt.axvline(x=pos, linestyle=':', color='k')

plt.xlabel('Time (minutes)')
plt.show()


    
    
