"""


"""


import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u
from astropy import constants as const
from astropy import convolution as conv
from astropy.table import Table, Column
import scipy.interpolate as scinterp
import scipy.stats as scistats
import time
import astroplan as ap

from astropy.visualization import time_support
time_support()


start = time.time()

plt.rc('lines',linewidth=0.5)
#plt.rc('font', size =18)

#print start
#import wdatmos
import spec_plot_tools as spt
import cal_params as cp
import plot_spec as ps

cerro_pachon_location = coord.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)

speclistname = 'listSpec'
speclist= np.genfromtxt(speclistname, dtype = 'str')
logg_file='c1_logg_differences.csv'
teff_file='c1_teff_differences.csv'

fitsfile='ravg_fwctb.EC01578m1743_930_blue_045asec.fits'

textfile=''

clist=['','r','purple','yellow','g','cyan','blue','orange']
name_string=''
color_index=0

teff_table=Table.read(teff_file)
logg_table=Table.read(logg_file)

time_list=[]
PA_list=[]
parallactic_list=[]
rotator_list=[]
#print(sorted(coord.EarthLocation.get_site_names()))
soar=ap.Observer.at_site('gems',timezone='utc')
#obs=ap.Observer(location
point_colors=[]
exptimes=[]
bmjd_list=[]

def get_parallactic(header):
    obs_time=header['DATE-OBS']
    target_coords=coord.SkyCoord(ra=header['ra'],dec=header['dec'], frame='icrs', unit=(u.hourangle,u.deg))
    parallactic=soar.parallactic_angle(obs_time,target_coords)
    #probably need to convert the parallactic angle from radians to degrees...
    parallactic=parallactic*180./np.pi
    #time_list.append(obs_time.mjd)
    time_list.append(obs_time)
    PA_list.append(header['posangle'])
    rotator_list.append(header['rotator'])
    parallactic_list.append(parallactic.value)
    exptimes.append(header['exptime'])
    bmjd_list.append(get_bmjd(header))

def get_bmjd(header):
    target_coord = coord.SkyCoord(header['ra'], header['dec'], frame = 'icrs', unit= (u.hourangle, u.deg) )
    obs_time=Time(header['DATE-OBS'],location=cerro_pachon_location,scale='tai')
    obs_time=obs_time+(header['exptime']*u.s/2.)
    bary_corr =obs_time.tdb.light_travel_time(target_coord)
    bmjd = (obs_time.tdb+ bary_corr.tdb).isot
    #bmjd = (obs_time.tdb-bary_corr.tdb).isot #wrong way to do this. You're supposed to add it but I'm trying to match timestamps at this point.
    
    return bmjd
    
starts=[]
ends=[]
i=0
for filename in speclist:
    get_parallactic(fits.getheader(filename))
    name_core=filename.split('.')[0][4:]
    print('\n======')
    print('name_core',name_core)
    print('name_string',name_string)
    print('color_index',color_index)
    
    if name_core==name_string:
        pass
    else:
        name_string=name_core
        color_index+=1
        print('New name_string',name_string)
        print('New color_index', color_index)
        starts.append(i)
        ends.append(i)
    print('color',clist[color_index])
    point_colors.append(clist[color_index])
    i+=1

ends=ends[1:]
ends.append(i)
print(starts)
print(ends)
    
point_colors=np.array(point_colors)
time_array=Time(time_list,scale='utc')
exp_array=np.array(exptimes)*u.s
#time_array=time_array+exp_array/2.

bmjd_array=Time(bmjd_list)
#print('type(bmjd_array)',type(bmjd_array))
pa_array=np.array(PA_list)
parallactic_array=np.array(parallactic_list)
neg_inds=np.where(parallactic_array<0)
parallactic_array[neg_inds]=360+parallactic_array[neg_inds]

rotator_array=np.array(rotator_list)
#print(time_array)
#print(pa_array)
#print(parallactic_array)
starttime=Time('2024-01-21T00:00:00')
endtime=Time('2024-01-21T23:59:00')



teff_times=Time(teff_table['timestamp_unix'],format='unix_tai')
logg_times=Time(logg_table['timestamp_unix'],format='unix')

inbounds=np.where((teff_times>starttime) & (teff_times<endtime))
#print('\n\nTimes for the correct night')
#print(teff_times[inbounds].isot)
#print('len(teff_times)',len(teff_times[inbounds]))
#print('time_array.shape',time_array.shape)
#print(teff_times.isot)
teff_times_diff=teff_times-np.roll(teff_times,1)
time_header_diffs=time_array-np.roll(time_array,1)
#for index in range(0,len(teff_times)):
    #print(teff_times[index].isot,teff_times_diff[index].to(u.s))
#for index in range(0,time_array.shape[0]):
    #print(time_array[index],time_header_diffs[index].to(u.s))
teff_list=[]
logg_list=[]
index=0
#for row in time_array:
#for row in bmjd_array:

    #time_diffs=row-teff_times
    #mintime=np.nanmin(np.abs(time_diffs))
    #timematch=np.argmin(np.abs(time_diffs))
    #print('minimum time diff', mintime.to(u.s),'exptime',exp_array[index])
    #print('min time index',timematch)
    #print(row,teff_times[timematch].isot,(row-teff_times[timematch]).to(u.s))
    #print('teff',teff_table['c1_teff_difference'][timematch])
    #print('logg',logg_table['c1_logg_difference'][timematch])
    #print('filename',speclist[index],'\n')
    #teff_list.append(teff_table['c1_teff_difference'][timematch])
    #logg_list.append(logg_table['c1_logg_difference'][timematch])
    #index+=1

teff_subtable=teff_table[inbounds]
logg_subtable=logg_table[inbounds]


test_header=fits.getheader(speclist[0])
test_target=coord.SkyCoord(ra=test_header['ra'], dec=test_header['dec'],unit=(u.hourangle,u.deg))
test_time=Time(test_header['date-obs'])
#ap.plots.plot_parallactic(test_target,soar,time_array)
plt.show()

#plt.scatter(time_array,pa_array, label='Position Angle')
#plt.scatter(time_array,parallactic_array, label='Parallactic Angle')

#plt.legend()
#plt.xlabel('Time')
#plt.ylabel('Angle (degrees)')
#plt.xticks(rotation=90)
#plt.show()


#fig, ax=plt.subplots(1,1)
#ax.scatter(time_array,pa_array-parallactic_array,c=point_colors)
#ax.set_xlabel('Time')
#ax.set_ylabel('Position Angle - Parallactic Angle (degrees)')
##ax.set_xticks(rotation=90)
#ax.xaxis.set(major_locator=mdates.YearLocator(),
             #major_formatter=mdates.DateFormatter("%Y"))
             
#teff_array=np.array(teff_list)
#logg_array=np.array(logg_list)
teff_subtimes=teff_times[inbounds]
logg_subtimes=logg_times[inbounds]
subteffs=teff_subtable['c1_teff_difference']
subloggs=logg_subtable['c1_logg_difference']
angle_off=pa_array-parallactic_array
#hold_teffs=np.copy(subteffs)

for first,last in zip(starts,ends):
    subteffs[first:last]=subteffs[first:last]-subteffs[first]
    subloggs[first:last]=subloggs[first:last]-subloggs[first]

inlier_inds=np.where(np.abs(subteffs)<3000.)

def trim(array):
    newarray=array[inlier_inds]
    return(newarray[:-4])

teff_subtimes=trim(teff_subtimes)
logg_subtimes=trim(logg_subtimes)
subteffs=trim(subteffs)
subloggs=trim(subloggs)
angle_off=trim(angle_off)
point_colors=trim(point_colors)




#angle_off=np.sin(angle_off/180.*np.pi)


sorted_order=np.argsort(angle_off)
angle_off=angle_off[sorted_order]
subteffs=subteffs[sorted_order]
subloggs=subloggs[sorted_order]
teff_subtimes=teff_subtimes[sorted_order]
logg_subtimes=logg_subtimes[sorted_order]
point_colors=point_colors[sorted_order]
#print('hold_teffs-sub_teffs',hold_teffs-subteffs)
nonnan_indices=~np.isnan(subteffs)
#nonnan_indices=nan_indices*-1
#print(nan_indices)
print(nonnan_indices)
print('subteffs')
print(subteffs)
print('subteffs nonnan_indices')
print(subteffs[nonnan_indices])


angle_teff_r=scistats.pearsonr(angle_off[nonnan_indices],subteffs[nonnan_indices])
angle_logg_r=scistats.pearsonr(angle_off[nonnan_indices],subloggs[nonnan_indices])
print('angle off and Teff R:',angle_teff_r)
print('angle off and logg R:',angle_logg_r)

plt.scatter(teff_subtimes,subteffs,c=point_colors)
#plt.scatter(bmjd_array,teff_array,c=point_colors,marker='*')



plt.xlabel('Time')
plt.ylabel('Drift in Teff Differences for a given Target (K)')
plt.xticks(rotation=90)

plt.show()

plt.scatter(logg_subtimes,subloggs,c=point_colors)



plt.xlabel('Time')
plt.ylabel('Drift in log(g) differences for a given Target')
plt.xticks(rotation=90)

plt.show()
             

plt.scatter(teff_subtimes,angle_off,c=point_colors)



plt.xlabel('Time')
plt.ylabel('Position Angle - Parallactic Angle (degrees)')
plt.xticks(rotation=90)

plt.show()

plt.scatter(angle_off,subloggs,c=point_colors)
plt.xlabel('Position Angle - Parallactic Angle (degrees)')
plt.ylabel('Drift in log(g) differences for a given Target')
plt.text(5,-0.15,'R:'+str(angle_logg_r[0])+'\nP='+str(angle_logg_r[1]))
plt.show()

plt.scatter(angle_off,subteffs,c=point_colors)
plt.xlabel('Position Angle - Parallactic Angle (degrees)')
plt.ylabel('Drift in Teff Differences for a given Target (K)')
plt.text(5,-1500,'R:'+str(angle_teff_r[0])+'\nP='+str(angle_teff_r[1]))
plt.show()




#plt.scatter(time_array,exptimes,c=point_colors)
#plt.xlabel('Time')
#plt.ylabel('Exposure times (seconds)')
#plt.xticks(rotation=90)

#plt.show()

#plt.scatter(time_array, rotator_list)
#plt.ylabel('rotator angle (degrees)')
#plt.xlabel('Time')
#plt.show()
    
    
    
    
    
    
