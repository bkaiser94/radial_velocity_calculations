"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-07-22


This is supposed to take a line list (I'm using one from NIST at the moment, but there's no reason it has to be
from there other than consistency of headers), and it uses that line list on a lamp image that is wavelength 
calibrated already using existing line lists and a different lamp image (one that doesn't have the new lines to be 
identified or actually might as a check).

This doesn't do any actual polynomial fitting. It should literally produce an intensity plot, and that's it. Well, it 
should also plot the line labels over top of it...





"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import balmer_line_ranges as blr
from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column
import scipy.interpolate as scinterp
import scipy.optimize as sciop

import cal_params as cp
import spec_plot_tools as spt


#output_filename= '400m2_HgArNe_calc.txt'
#linelist_file= '400M2_HgAr.txt'
#nelist_file='NIST_NeI_linelist_copy.csv'
#nelist_file='NIST_NeI_linelist_copy.csv'
#nelist_file='NIST_HgIArINeI_480to920_vac.csv'
#nelist_file='NIST_HgIArINeI_480to920_air_Cup.csv'
#nelist_file='NIST_HgIArINeI_350to750_air_Cup.csv'
#wave_soln_file= 'wsoln_0114_SDSSJ1150p2403_400m2.txt'
#fits_file= 'ctb.0059_HgAr_400m2_simple_fe.fits'
#fits_file2= 'ctb.0056_HgArNe_400m2_simple_fe.fits'


output_filename= '400m1_HgArNe_calc.txt'
linelist_file= '400M1_HgAr_redcam.txt'
nelist_file='NIST_HgIArINeI_350to750_air_Cup.csv'
wave_soln_file='wsoln_0213_SDSSJ1150p2403_400m1.txt'
fits_file='ctb.0052_HgAr_400m1_simple_fe.fits'
fits_file2='ctb.0049_HgArNe_400m1_simple_fe.fits'

wavelength_key= 'obs_wl_air(A)'
#wavelength_key= 'obs_wl_vac(A)'

wave_sol_binning= 2 #binning of the wavelength solution

trace_width=10
trace_mid=100

intensity_threshold= np.float_(499. )#minimum "intens" value that the NIST lines can have to make it into the plotting
element_sel= 'Ar'

wave_soln_file= cp.wave_sol_dir+wave_soln_file
linelist_file= cp.line_list_dir+linelist_file
nelist_file= cp.line_list_dir+nelist_file

color_dict={
    'Ar':'b',
    'Hg':'r',
    'Ne':'g'
    }
############################

def get_binning(header):
    """
    get the dispersion direction binning of pixels to figure out how to map a polynomial from the wavelength calibrations.
    """
    binning= header['CCDSUM'].split(' ')
    xbinning= binning[0]
    return xbinning



def wavelength_to_pixel(lambda_val, in_wave_coeffs, lamp_poly_degree=5, bounds=[0,2100]):
    """
    input wave_coeffs should already have an offset subtracted from the x-values everywhere.... you can't really 
    do that...
    
    bounds: the pixel boundaries at which the lambda_val might be located (it should be 
    """
    wave_coeffs= np.copy(in_wave_coeffs)
    wave_coeffs[-1]= wave_coeffs[-1]-lambda_val
    
    def func_to_solve(x):
        if lamp_poly_degree==5:
            return wave_coeffs[0]*x**5+ wave_coeffs[1]*x**4 +wave_coeffs[2]*x**3+ wave_coeffs[3]*x**2+wave_coeffs[4]*x + wave_coeffs[5]
        else:
            print("don't have function to solve for inversion of wavelengths for that lamp_poly_degree:", lamp_poly_degree)
            return np.polyval(wave_coeffs, x)
    #plt.plot(np.polyval(wave_coeffs, np.linspace(0,2000,2000)), label='changed wave_coeffs')
    #plt.plot(np.polyval(in_wave_coeffs, np.linspace(0,2000,2000)), label='og wave_coeffs')
    #plt.plot(func_to_solve(np.linspace(0,2000,2000)),label='func_to_solve')
    #plt.legend(loc='best')
    #plt.show()
    pixel= sciop.brentq(func_to_solve, bounds[0],bounds[1])
    return pixel


def retrieve_ben_list():
    fear_array= np.genfromtxt(linelist_file, names = True)
    line_x_checks = np.copy(fear_array['Pixel'])
    print "line_x_checks should have just been created"
    print line_x_checks
    lamp_lines = np.copy(fear_array['User'])
    #line_sides = np.ones(line_x_checks.shape[0])*line_search_width
    names= np.str_(lamp_lines)
    #try:
        #good_inds= np.where(fear_array['use']>0)
        #line_x_checks=line_x_checks[good_inds]
        #lamp_lines= lamp_lines[good_inds]
    #except KeyError:
        #pass
    return line_x_checks, lamp_lines


def remove_stars(ne_table):
    for row in ne_table:
        print(row['intens'])
        if '*' in row['intens']:
            print('* in intens')
            row['intens']='0'
        elif ''==row['intens']:
            print('blank row')
            row['intens']='0'
        else:
            pass
    return ne_table

def retrieve_nist_list():
    ne_table=Table.read(nelist_file, format='ascii.csv')
    #print(ne_table['intens'])
    #ne_table=remove_stars(ne_table)
    #print(ne_table['intens'])
    #try:
        #good_inds= np.where(ne_table['use']>0)
        #ne_table=ne_table[good_inds]
    #except KeyError:
        #pass
    #for count, row in enumerate(ne_table['obs_wl_air(A)']):
        #print(count, row)
    for count, row in enumerate(ne_table):
        print(count, row['intens'])
        try:
            row['intens']=int(row['intens'])
            ne_table[count]['intens']=float(row['intens'])
            print('converted to int!', row['intens'])
        except ValueError:
            print('ValueError ^^^^')
    #ne_table['intens']=ne_table['intens'].astype(float) #changing this column that gets read as strings for whatever reason
    ne_table.pprint()
    return ne_table

def inbounds_table(ne_table, im_waves):
    lower_indices = np.where(ne_table[wavelength_key] < np.nanmax(im_waves))
    print(np.nanmax(im_waves).shape)
    ne_table=ne_table[lower_indices]
    print('less than max', np.nanmax(im_waves))
    ne_table.pprint()
    upper_indices= np.where(ne_table[wavelength_key] > np.nanmin(im_waves))
    ne_table=ne_table[upper_indices]
    print('greater than min', np.nanmin(im_waves))
    ne_table.pprint()
    print(ne_table['intens'].dtype)
    intense_indices = np.where(ne_table['intens'] > intensity_threshold)
    ne_table=ne_table[intense_indices]
    print('greater than threshold', intensity_threshold)
    ne_table.pprint()
    return ne_table

def element_table(input_table, element_string):
    element_inds= np.where(input_table['element']==element_string)
    output_table= input_table[element_inds]
    
    return output_table

def plot_ne_lines(ne_table, counts):
    for row in ne_table:
        #print(name+name2, type(name))
        try:
            air_name=row['element']+str(row['sp_num'])+'-'+str(row[wavelength_key])[:4]
        except KeyError as error:
            print('KeyError:', error)
            air_name='NeI '+str(row[wavelength_key])[:4]
        plt.axvline(x=row[wavelength_key], linestyle='--', color=color_dict[row['element']])
        #plt.text(row['obs_wl_air(A)'], np.nanmax(counts), air_name, color='g', rotation=90)
        plt.text(row[wavelength_key], 1000, air_name, color=color_dict[row['element']], rotation=90)
    return


def plot_ne_pixels(ne_table, ne_pixels, counts):
    for row,pixel in zip(ne_table, ne_pixels):
        #print(name+name2, type(name))
        air_name='NeI '+str(row[wavelength_key])[:4]
        plt.axvline(x=pixel, linestyle='--', color='g')
        #plt.text(row['obs_wl_air(A)'], np.nanmax(counts), air_name, color='g', rotation=90)
        plt.text(pixel, 1000, air_name, color='g', rotation=90)
    return

######################


hdu = fits.open(fits_file)
header= fits.getheader(fits_file)
image_data = hdu[0].data

hdu2=fits.open(fits_file2)
header2=fits.getheader(fits_file2)
image_data2=hdu2[0].data


roi= image_data[trace_mid-trace_width/2: trace_mid+trace_width/2,:]
counts= np.nanmean(roi,axis=0)


roi2= image_data2[trace_mid-trace_width/2: trace_mid+trace_width/2,:]
counts2= np.nanmean(roi2,axis=0)

band_inds= np.indices(counts.shape)[0]
print(band_inds.shape)


wave_poly_coeffs= np.genfromtxt(wave_soln_file)

wave_vals= np.polyval(wave_poly_coeffs, band_inds)


hgar_pixel_vals, hgar_wave_vals= retrieve_ben_list()
ne_table= retrieve_nist_list()
ne_table= inbounds_table(ne_table, wave_vals)
hg_table= element_table(ne_table, 'Hg')
ar_table= element_table(ne_table, 'Ar')
ne_table= element_table(ne_table,'Ne')
#ne_table= element_table(ne_table, element_sel)

ne_table.pprint()
for x_spot in hgar_pixel_vals:
    plt.axvline( x= x_spot, color = 'k',linestyle = '--')
plt.plot(counts)
plt.ylabel('intensity')
plt.xlabel('pixel')
plt.show()


for x_spot in hgar_wave_vals:
    plt.axvline( x= x_spot, color = 'k',linestyle = '--')
plot_ne_lines(ne_table, counts)
plot_ne_lines(ar_table, counts)
plot_ne_lines(hg_table, counts)

plt.plot(wave_vals,counts2, color='gray')
plt.plot(wave_vals,counts, color='k')
plt.ylabel('intensity')
#plt.yscale('log')
plt.xlabel('wavelength (angstroms)')
plt.show()


ne_pixels=[]
for row in ne_table:
    ne_wave= float(row[wavelength_key])
    ne_pixel= wavelength_to_pixel(ne_wave, wave_poly_coeffs)
    ne_pixels.append(ne_pixel)

hgar_calc_pixels= []
for hgar_wave in hgar_wave_vals:
    hgar_pixel= wavelength_to_pixel(hgar_wave, wave_poly_coeffs)
    hgar_calc_pixels.append(hgar_pixel)
    

plot_ne_pixels(ne_table, ne_pixels, counts)
for x_spot in hgar_calc_pixels:
    plt.axvline( x= x_spot, color = 'k',linestyle = '--')
plt.plot(counts2, color='r')
plt.plot(counts)
plt.xlabel('pixel')
plt.show()


ne_waves= np.copy(ne_table[wavelength_key])
ne_pixels= ne_pixels
ne_use= np.copy(ne_table['use'])
ne_relheight= np.copy(ne_table['intens'])
ne_names= []
ne_name2= []
ne_other= []
ne_residuals= np.zeros(ne_use.shape)
for row in ne_table:
    ne_name= 'NeI-'+str(row[wavelength_key])[:4]
    ne_names.append(ne_name)
    ne_name2.append(0)
    ne_other.append(0)
    
hgar_array= np.genfromtxt(linelist_file, names = True)
    
hgar_waves= np.copy(hgar_wave_vals)
hgar_pixels= np.copy(hgar_calc_pixels)
hgar_use= np.copy(hgar_array['use'])
hgar_residuals= np.copy(hgar_array['Residual'])

hgar_name= []
for row in hgar_waves:
    hgar_name.append('HgAr-'+str(row)[:4])

hgar_name2= hgar_array['Name2']
hgar_relheight= hgar_array['Rel_Height']
hgar_other= hgar_array['other']



######################3
#now put those arrays together and make some stuff

all_waves= np.append(hgar_waves, ne_waves)
print(all_waves.shape)
print(all_waves)
all_pixels= np.append(hgar_pixels, ne_pixels)
all_residual=np.append(hgar_residuals, ne_residuals)
all_name= np.append(hgar_name, ne_names)
all_name2= np.append(hgar_name2, ne_name2)
all_relheight= np.append(hgar_relheight, ne_relheight)
all_other= np.append(hgar_other, ne_other)
all_use= np.append(hgar_use, ne_use)

####################
#make the Columns and table

all_table=Table()
all_table['User']= all_waves
all_table['Pixel']=all_pixels
all_table['Residual']=all_residual
all_table['Name']=all_name
all_table['Name2']=all_name2
all_table['Rel_Height']=all_relheight
all_table['other']=all_other
all_table['use']=all_use

all_table.pprint()

print('saving ', output_filename)
#all_table.write(output_filename, format='ascii', delimiter='\t')


#ne_table.pprint()

#plt.plot(ne_table['obs_wl_air(A)'], ne_table['intens'])
#plt.show()


