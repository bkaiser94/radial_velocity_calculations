"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-06-06, D-Day

This should be a test-bed for the more complicated atmospheric correction processes that I'll be undertaking
I'm not sure if this is going to evolve into a separate script for atmospheric corrections or if it will be copied and 
pasted into flux_calibration.py and calibrate_flux.py; I just don't want to clutter those even further with
unnecessary plotting as I'm about to do here.

"""


from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
from astropy.time import Time
from astropy import coordinates as coords
from astropy.table import Table, hstack
from astropy import units as u
from astropy import constants as const
from astropy import convolution as conv
import scipy.interpolate as scinterp
import csv
import time
start = time.time()

import spec_plot_tools as spt
import cal_params as cp
import get_cal_params as gcp
import plot_spec as ps

#sens_names = glob('E*sensitivity_curve.txt')+glob('G*sensitivity_curve.txt')+glob('e*sens*curve.txt')
sens_names=glob('sens_curv*')
#sens_names=glob('sens*')
#sens_names = glob('*sensitivity*.txt')

resid_names= glob('resid*.txt')
#resid_names=glob('tell*')
tell_names= glob('tell*')
do_fnu=False

sens_names=sorted(sens_names)
resid_names=sorted(resid_names)
tell_names=sorted(tell_names)

#wavelengths = np.linspace(4940,8980, 8080) #400M2 approximately
#wavelengths = np.linspace(3800,7200, 8080) #400M1 approximately
wavelengths= np.linspace(3800, 8980, 10000)


def extract_AM_MJD(sens_curve_file):
    """
    INPUT: filename string for one of the sensitivity curve files
    
    OUTPUT: tuple of the airmass and MJD value for the sensitivity curve file from its 'header'
    
    """
    with open(sens_curve_file, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        #print('reader[0]',reader[0])
        index=0
        header_dict={}
        for row in reader:
            #print(row)
            if index==0:
                ##print(row[0])
                ##airmass_string=row[0][-5:]
                #airmass_string= row[0].split(':')[1]
                ##print('airmass:',airmass_string)
                ##airmass_string=airmass_string.replace(' ', '')
                #airmass= float(airmass_string)
                #mjd_string= row[1].split(':')[1]
                #mjd=float(mjd_string)
                ##print('mjd:', mjd_string)
                ##print('mjd:', mjd)
                for entry in row:
                    #print('entry', entry)
                    entry=entry.replace('#', '')
                    parts = entry.split(':')
                    parts[0]=parts[0].replace(' ', '')
                    #print('parts', parts)
                    try:
                        header_dict[parts[0]]=float(parts[1])
                    except ValueError as error:
                        print("ValueError:", error)
                        print("That was presumably the newly added name of the standard that was used as the reference spectrum, which BETTER match the name of the sens_curv file! I'll just ignore this for now since by default things were trying to cast to floats for the header.")
            else:
                pass
            index+=1
    
    return header_dict

def get_absorption(tell_spec, wave_range):
    absorptions= 1.-tell_spec[1]
    inbounds= np.where((tell_spec[0]>wave_range[0]) & (tell_spec[0] < wave_range[1]))
    total_abs= np.sum(absorptions[inbounds])
    total_bins= tell_spec[1][inbounds].shape[0]
    dlambdas=np.roll(tell_spec[0], -1)-tell_spec[0]
    dlambda_in=dlambdas[inbounds]
    dlambda_in[np.where(dlambda_in <0)] = np.nan
    ew= np.nansum(absorptions[inbounds]*dlambda_in)/1.
    print('ew', ew)
    avg_abs= total_abs/total_bins
    return ew

def generate_abs_row(tell_spec):
    absorp_vals= []
    for tell_range in cp.telluric_lines:
        this_absorp = get_absorption(tell_spec, tell_range)
        absorp_vals.append(this_absorp)
        #print(tell_range, this_absorp)
    return absorp_vals

coll_max_thrus= []

for sens_name in sens_names:
    #airmass, mjd= extract_AM_MJD(sens_name)
    header_dict= extract_AM_MJD(sens_name)
    airmass= header_dict['Airmass']
    mjd= header_dict['MJD']
    sens_curve_coeffs = np.genfromtxt(sens_name)
    sens_curve = np.polyval(sens_curve_coeffs,wavelengths)
    max_index= np.argmax(sens_curve)
    label= ','.join([sens_name, str(airmass), str(mjd)])
    coll_max_thrus.append([mjd, airmass, wavelengths[max_index], sens_curve[max_index]])
    plt.plot(wavelengths, sens_curve, label=label)
    plt.scatter(wavelengths[max_index], sens_curve[max_index], marker='*')
    plt.text(wavelengths[max_index],sens_curve[max_index],sens_name)
    print("\n=============")
    print(sens_name)
    print('airmass:', airmass, 'mjd:', mjd)
    print('Peak flux at', wavelengths[max_index], 'angstroms')
plt.xlabel(r'wavelength ($\AA$)')
plt.xlim(np.nanmin(wavelengths), np.nanmax(wavelengths))
plt.ylim(0,1)
spt.show_plot(show_legend=False)

coll_max_thrus=np.array(coll_max_thrus).T
plt.plot(coll_max_thrus[1], coll_max_thrus[3], linestyle='none', marker='o')
plt.xlabel('airmass')
plt.ylabel('max throughput of sensitivity curve')
plt.show()

for resid_name in resid_names:
    header_dict= extract_AM_MJD(resid_name)
    airmass= header_dict['Airmass']
    mjd=header_dict['MJD']
    resid_array = np.genfromtxt(resid_name)
    max_index= np.argmax(sens_curve)
    label= ','.join([resid_name, str(airmass), str(mjd)])
    plt.plot(resid_array, label=label)
    print("\n=============")
    print(sens_name)
    print('airmass:', airmass, 'mjd:', mjd)
plt.xlabel(r'Pixel')
spt.show_plot()

tell1_array=  np.genfromtxt(tell_names[0], skip_header=1).T
tell1_waves= tell1_array[0]
counter=0


coll_absorps= []
header_table=[]
header_names=[]

header_dict= extract_AM_MJD(tell_names[0])
for thing in header_dict:
    header_names.append(thing)


print('header_names',type(header_names), header_names)


#ps.plot_telluric_spectrum([3700,9000], smooth=False, pix_width=30, tell_filename='LBL_A30_s0_w015_R0060000_T.fits')
#ps.plot_telluric_spectrum([3700,9000], smooth=True, pix_width=30, tell_filename='LBL_A30_s0_w015_R0060000_T.fits')
pix_width_tell=30
water='050'
ps.plot_telluric_spectrum([3700,9000], smooth=True, pix_width=pix_width_tell, tell_filename='LBL_A10_s0_w'+water+'_R0060000_T.fits',show_filename=True)
ps.plot_telluric_spectrum([3700,9000], smooth=True, pix_width=pix_width_tell, tell_filename='LBL_A15_s0_w'+water+'_R0060000_T.fits',show_filename=True)
ps.plot_telluric_spectrum([3700,9000], smooth=True, pix_width=pix_width_tell, tell_filename='LBL_A20_s0_w'+water+'_R0060000_T.fits',show_filename=True)
#ps.plot_telluric_spectrum([3700,9000], smooth=True, pix_width=pix_width_tell, tell_filename='LBL_A25_s0_w'+water+'_R0060000_T.fits',show_filename=True)
#ps.plot_telluric_spectrum([3700,9000], smooth=True, pix_width=pix_width_tell, tell_filename='LBL_A30_s0_w'+water+'_R0060000_T.fits', show_filename=True)


for tell_name in tell_names:
    #airmass, mjd= extract_AM_MJD(tell_name)
    header_dict= extract_AM_MJD(tell_name)
    row_add= []
    for thing in header_dict:
        row_add.append(header_dict[thing])
    header_table.append(row_add)
    airmass= header_dict['Airmass']
    mjd= header_dict['MJD']
    tell_array = np.genfromtxt(tell_name, skip_header=1).T
    absorp_row= generate_abs_row(tell_array)
    coll_row= np.hstack([mjd,airmass, absorp_row])
    coll_absorps.append(coll_row)
    #tell_array[1]=tell_array[1]+counter
    subname= tell_name.split('_')[4:6]
    subname='_'.join(subname)
    label= ','.join([subname, str(airmass), str(mjd)])
    plt.plot(tell_array[0], tell_array[1], label=label)
    #plt.plot(tell_array[0]+2.7, tell_array[1], label=label)
    #plt.axhline(y=1, color='k', linestyle='--')
    #plt.xlabel(r'Wavelength ($\AA$)')
    #plt.ylabel('Transmission')
    #plt.title('Telluric factors')
    #spt.show_plot(show_legend=False)
    
    #plt.plot(tell_array[1], label=label)
    #plt.plot(tell1_waves, np.interp(tell1_waves, tell_array[0], tell_array[1]), label='interp_'+label)
    #plt.text(np.nanmin(tell_array[0]), counter+1, label, color='k')
    #plt.text(0, 1, label, color='k')
    print("\n=============")
    print(tell_name)
    print('airmass:', airmass, 'mjd:', mjd)
    counter+=1
header_table=np.array(header_table)
print('header_table.shape', header_table.shape)
print('len(header_names)', len(header_names))
header_table= Table(header_table, names= header_names)
header_table.pprint()

coll_absorps= np.array(coll_absorps).T

plt.axvline(x=7599, color='k', linestyle='--')
plt.axvline(x=7623, color='k', linestyle='--')
plt.axvline(x=6866, color='k', linestyle='--')
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Transmission')
spt.show_plot(show_legend=True)
#plt.show()

index= 0

#for tell_region in coll_absorps[2:]:
    #plt.scatter(coll_absorps[1], tell_region, label='region '+ str(cp.telluric_lines[index]))
    #index+=1
for tell_region in coll_absorps[2:]:
    plt.plot(coll_absorps[1], tell_region, label='region '+ str(cp.telluric_lines[index]), linestyle='None', marker='o')
    index+=1
plt.legend(loc='best')
plt.xlabel('Airmass')
plt.ylabel(r'Equivalent width of absorption ($\AA$)')
plt.title("Telluric absorption EW's in 400M2")
plt.grid()
plt.show()

plt.plot(coll_absorps[1], coll_absorps[0], linestyle='None', marker='o')
plt.ylabel('MJD')
plt.xlabel('airmass')
plt.grid()
plt.show()

#plt.plot(header_table['Ext_width']/header_table['See_FWHM'], coll_max_thrus[3], linestyle='none', marker='o')
plt.plot(header_table['Ext_width'], coll_max_thrus[3], linestyle='none', marker='o')
plt.xlabel('Extraction_width')
plt.ylabel('Maximum throughput')
plt.show()

plt.plot(header_table['See_FWHM'], coll_max_thrus[3], linestyle='none', marker='o')
plt.xlabel('FWHM')
plt.ylabel('Maximum throughput')
plt.grid()
plt.show()



coll_absorps_names=[]
coll_absorps= coll_absorps[2:] #remove MJD and airmass
for num, col in enumerate(coll_absorps):
    coll_absorps_names.append('EW_'+str(num))
print('coll_absorps_names', coll_absorps_names, len(coll_absorps_names))
print('coll_absorps')
print(coll_absorps.shape)

coll_absorps= coll_absorps.T
print('coll_absorps', coll_absorps.shape)
print('coll_absorps_names', len(coll_absorps_names))
coll_absorps_table= Table(coll_absorps, names= coll_absorps_names)
output_table= hstack([header_table, coll_absorps_table])

#output_table.write('collected_EWs.csv', format='ascii.csv', overwrite=True)





