"""
written by Ben Kaiser (UNC - Chapel Hill)

:INPUTS:
    RA : decimal degree value provided as commandline input
    DEC: decimal degree value provided as commandline input

:OUTPUTS:
    print plots showing the times the target should be observed to get the radial velocity values in quadrature
    


*** There are orbital parameters that are contained within the variables in this file that have to be edited in file in order for this script to work with a system other than the one it was originally written for. The RA and DEC are required as inputs so the target isn't obvious to  anyone but the person running it.****
"""

import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
from astropy.units import cds
cds.enable()
#plt.rc('font', size =18)
#plt.rc('lines', markersize=12)
plt.rc('font', size = 12)
plt.rc('lines', markersize = 12)

ra = float(sys.argv[1]) #values in decimal degrees
dec = float(sys.argv[2])

parkes_location = coord.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
cerro_pachon_location = coord.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)

#times = ['2018-03-19T02:19:00', '2018-03-19T10:11:00']
#times =['2018-04-23T02:00:00','2018-04-23T10:20:00']
#times =['2018-04-24T02:00:00','2018-04-24T10:20:00']
#times = ['2018-03-27T01:45:00', '2018-03-27T10:04:00']
#times =['2018-07-3T22:30:00','2018-07-04T06:10:00']
#times =['2018-07-08T23:00:00','2018-07-09T06:00:00']
times =['2018-08-07T23:00:00','2018-08-08T04:15:00']
#times =['2018-08-16T23:00:00','2018-08-17T06:00:00']


#time_of_interest=[['2018-07-08T23:30:00','2018-07-09T00:45:00'],['2018-07-09T04:15:00','2018-07-09T05:30:00']]
time_of_interest=[['2018-08-08T00:15:00','2018-08-08T01:15:00'],['2018-08-08T04:15:00','2018-08-08T04:16:00']]
#time_of_interest=[['2018-08-16T23:30:00','2018-08-17T00:45:00'],['2018-08-17T04:15:00','2018-08-17T05:30:00']]


obs_times= Time(times, format = 'isot', scale='utc', location = cerro_pachon_location)
target_coord = coord.SkyCoord(ra, dec, unit= (u.deg, u.deg), frame= 'icrs')

def to_barycenter(input_times):
    bary_corr =input_times.tdb.light_travel_time(target_coord)
    return (input_times.tdb+ bary_corr.tdb).mjd

q_value = 0.88/1.35
lam_rest = 6562.81*u.angstrom #angstroms
#m1 = (0.14 *u.Msun).si
#m2 = (1.4 *u.Msun).si
#m1 = (1.4*u.Msun).si
m1= (1.4*u.Msun).si
#m2 = (0.14*u.Msun).si #switched the masses to change the phase by 180 degrees.
#m2=( 1.58 *u.Msun).si
m2 = (m1*q_value).si
#m2 = (0.12*u.Msun).si #switched the masses to change the phase by 180 degrees.

#e=1e-5
e = 0.000010 #PSR J1435-6100
#e = 0.000023# median
#e = 0.000031# max
#e=0.5
#e=0.0011494 #1227
#omega = 97 * u.degree
#omega = 10 *u.degree #PSR J1435-6100
omega = 4 *u.degree #PSR J1435-6100

#period = (0.4497391377 * u.day).to(u.second) #median
#period = (0.449739137 * u.day).to(u.second) #low
#period = (0.4497391384 * u.day).to(u.second) #high
#period = (6.721013337  *u.day).to(u.second)# for a different target
period = (1.354885217 *u.day).to(u.second)  #PSR J1435-6100

#t0 = Time(55756.23, format = 'mjd', scale= 'utc', location = parkes_location)
#tasc = Time(55756.1047771, format = 'mjd', scale= 'utc', location = parkes_location) #median
##tasc = Time(55756.1047767, format = 'mjd', scale= 'utc', location = parkes_location) #min
#tconj= Time(55756.21712, format = 'mjd', scale ='utc', location = parkes_location)

#PSR J1435-6100
#t0 = Time(55756.23, format = 'mjd', scale= 'utc', location = parkes_location)
tasc = Time(51270.6084449, format = 'mjd', scale= 'utc', location = parkes_location) #median
#tasc = Time(55756.1047767, format = 'mjd', scale= 'utc', location = parkes_location) #min
#tconj= Time(55756.21712, format = 'mjd', scale ='utc', location = parkes_location)

def get_tconj(tasc, pb, e, omega):
    omega = omega.to(u.radian)
    return tasc + pb/4. + 2*(e*np.cos(omega))/(2*np.pi/pb)

tconj= get_tconj(tasc, period, e, omega)
print "tconj: ", tconj
print "tasc: ", tasc
#epoch_difference = obs_times[0].mjd-tasc.mjd
epoch_difference = obs_times[0].mjd-tconj.mjd
print "epoch difference: ", epoch_difference


bmjd_obs = to_barycenter(obs_times) #corrected to barycenter to use against the rv curve
#bmjd_t0 = to_barycenter(t0) #corrected initial epoch
bmjd_tasc= to_barycenter(tasc)
bmjd_tconj = to_barycenter(tconj)

#calculate nearest bmjd to start for calculating radial velocities around the observations
#time_dif= bmjd_obs[0]-bmjd_tasc #days between epoch and beginning of observability
time_dif = bmjd_obs[0]-bmjd_tconj
#nearest_time = int(time_dif/period.to(u.day).value)*period.to(u.day).value+bmjd_tasc #beginning of the orbital cycle
nearest_time = int(time_dif/period.to(u.day).value)*period.to(u.day).value+bmjd_tconj #beginning of the orbital cycle

print ("obs start", bmjd_obs[0])
print ("nearest start time", nearest_time)
G = u.cds.G.si
print (G)
def get_a(period, m1, m2):
    return ((period**2/(4*np.pi**2)*G*(m1+m2))**(1./3)).to(u.meter)

print (get_a(period, m1, m2))

def L(m1,m2,a,e):
    return np.sqrt(G*(m1+m2)*a*(1-e**2))

def r(th,a,e):
    return a*(1-e**2)/(1+e*np.cos(th))

def p(a,m1,m2):
    return 2*np.pi*np.sqrt(a**3/(G*(m1+m2)))

def x1(th,a,e,m1,m2):
    return (-m2/(m1+m2)*r(th,a,e)*np.cos(th)).to(u.meter)

#set the range to plot
th_max =(bmjd_obs[1]-nearest_time)/period.to(u.day).value*2*np.pi
print ("max phase:", th_max)
N= 10000
#th= np.linspace(0.,th_max,N)
th= np.linspace(0, 4*np.pi, N)
a_thing=get_a(period,m1,m2)
dt= (r(th,a_thing,e)**2/L(m1,m2,a_thing,e)*4*np.pi/N).to(u.second)
print ("dt.unit", dt.unit)
# print (dt)
# def v(a,b,dt):
#     return(a-b)/dt
#actually I don't want to use this velocity definition. I'm gonna do it with arrays. No for loops.

days=np.cumsum(dt).to(u.day)
def x2(th,a,e,m1,m2):
    return (m1/(m1+m2)*r(th,a,e)*np.cos(th)).to(u.meter)
def v(th,a,e,m1,m2,dt, body_num):
    """velocity function that works for either body, but it requires an input of body number"""
    if body_num== 1:
        x_i= x1(th,a,e,m1,m2)
    elif body_num == 2:
        x_i= x2(th,a,e,m1,m2)
    temp_cop= np.copy(x_i)
    x_i1= np.roll(x_i,1) #shifts all entries down one index
    vel= (x_i1-x_i)/dt
    return vel
    
def get_quad_points(velocities):
    zeros = np.argsort(velocities**2)[0:6]
    first_max = np.argsort(velocities)[-3:]
    first_min= np.argsort(velocities)[0:3]
    return np.hstack([zeros,first_max, first_min])

v1= v(th, a_thing, e, m1, m2, dt, 1)
v1 = v1.to(u.km/u.s)
v2= v(th, a_thing , e, m1, m2, dt, 2)
v2 = v2.to(u.km/u.s)

print ("Maximum radial velocity companion:", np.nanmax(v2))
bmjd_times =Time(nearest_time+ days.value,  format = 'mjd', scale='tdb', location = cerro_pachon_location)
mjd_times= bmjd_times+ bmjd_times.light_travel_time(target_coord)
utc_times = mjd_times.utc.mjd
utc_difs =((utc_times-obs_times[0].mjd)*u.day).to(u.hour)


def get_points_of_interest(time_of_interest=time_of_interest):
    utc_ranges= []
    v2_ranges = []
    for times1 in time_of_interest:
        times1= Time(times1, format = 'isot', scale='utc', location = cerro_pachon_location)
        allowed= np.where(utc_times < times1[1].utc.mjd)
        utc_times_interest= np.copy(utc_times[allowed])
        v2_interest = np.copy(v2[allowed])
        allowed2= np.where(utc_times_interest > times1[0].utc.mjd)
        utc_times_interest= utc_times_interest[allowed2]
        utc_times_interest =((utc_times_interest-obs_times[0].mjd)*u.day).to(u.hour)
        v2_interest= v2_interest[allowed2]
        utc_ranges.append([utc_times_interest])
        v2_ranges.append([v2_interest])
    return utc_ranges, v2_ranges


quad_points = get_quad_points(v2)
quad_days = utc_times[quad_points]
quad_hours = utc_difs[quad_points]
quad_times = Time(((quad_hours.to(u.day) +obs_times[0].mjd*u.day).to(u.day)).value, format='mjd', scale= 'utc').iso

for plusmin,actual_time, day_value  in zip(quad_hours, quad_times, quad_days):
    print plusmin, actual_time, day_value
plt.xlabel(r't (hours)')
plt.ylabel(r'v ('+str(v1.unit)+')')
plt.axvline(x =( (obs_times[0].mjd-obs_times[0].mjd)*u.day).to(u.hour).value, color = 'k', label = times[0], linestyle = '--')
plt.axvline(x=( (obs_times[1].mjd-obs_times[0].mjd)*u.day).to(u.hour).value, color = 'r', label = times[1], linestyle = '--')
plt.plot(utc_difs,v1,label='NS');
plt.plot(utc_difs,v2,label= 'Comp');
plt.plot(utc_difs[quad_points], v2[quad_points], marker = '*', color ='k', linestyle = 'None')


utc_times_of_interest, v2_of_interest= get_points_of_interest()
#print("utc_times_of_interest[0]", utc_times_of_interest[0])
#print(utc_times_of_interest)
for this_time, this_v2 in zip(utc_times_of_interest[0], v2_of_interest[0]):
    plt.plot(this_time, this_v2, label = 'Time of Interest', color = 'r')
    print("Max v2: ", np.nanmax(this_v2), " in ", this_time[0], this_time[-1])
    print("Min v2: ", np.nanmin(this_v2), " in ", this_time[0], this_time[-1])
    
for this_time, this_v2 in zip(utc_times_of_interest[1], v2_of_interest[1]):
    plt.plot(this_time, this_v2, label = 'Time of Interest', color = 'r')
    print("Max v2: ", np.nanmax(this_v2), " in ", this_time[0], this_time[-1])
    print("Min v2: ", np.nanmin(this_v2), " in ", this_time[0], this_time[-1])


plt.legend();
plt.title('e= '+ str(e)+ ',  a= '+str(a_thing.to(u.au))+ ',  $m_1$= '+str(m1.to(u.Msun))+',  $m_2$= '+str(m2.to(u.Msun)))
plt.show()
print (a_thing.to(u.au))
print ("a in ls", (a_thing/(u.cds.c.si)).to(u.second))


def plot_datapoints():
    try:
        plotfile = glob(sys.argv[3])
        
        all_array = np.genfromtxt(plotfile[0]).T
        mjd_array = all_array[0]
        H_delta = all_array[1]
        H_gamma = all_array[2]
        H_beta = all_array[3]
        H_delta_s = all_array[4]
        H_gamma_s = all_array[5]
        H_beta_s = all_array[6]
        #print Time(mjd_array, format = 'mjd').utc.isot
        #print all_array
        mean_rv = np.copy(np.mean(all_array[1:, :], axis = 0))
        std_dev =np.copy( np.std(all_array[1:, :], axis =0))
        #print mean_rv
        def zero_rvs(rv_array):
            print "systemic velocity (includes Earth's motion):", np.mean([rv_array.max(),rv_array.min()])
            return rv_array-np.mean([rv_array.max(),rv_array.min()])
        H_delta = zero_rvs(H_delta)
        H_gamma = zero_rvs(H_gamma)
        H_beta = zero_rvs(H_beta)
        remean_rv = np.mean([H_delta, H_gamma, H_beta], axis = 0)
        remean_std = np.std([H_delta,H_gamma, H_beta], axis=0)
        mean_rv = zero_rvs(mean_rv)
        #plt.plot(mjd_array, H_delta, label = r"H-$\delta$", linestyle = 'none', marker = '*')
        #plt.plot(mjd_array, H_gamma, label = r"H-$\gamma$", linestyle = 'none', marker = '*')
        #plt.plot(mjd_array, H_beta, label = r"H-$\beta$", linestyle = 'none', marker = '*')
        plt.errorbar(mjd_array, H_delta, H_delta_s, label = r"H-$\delta$", linestyle = 'none', marker = '*')
        plt.errorbar(mjd_array, H_gamma, H_gamma_s, label = r"H-$\gamma$", linestyle = 'none', marker = '*')
        plt.errorbar(mjd_array, H_beta, H_beta_s, label = r"H-$\beta$", linestyle = 'none', marker = '*')
        #plt.errorbar(mjd_array, mean_rv, std_dev, label = 'Mean RV', linestyle = 'none', marker = 'o')
        plt.errorbar(mjd_array, remean_rv, remean_std, label = r"Mean of zeroed RV's", linestyle = 'none', marker = 'o')
    except IndexError as error:
        print "No data file provided, so just outputting the predictions."
    return


plt.xlabel(r'MJD')
plt.ylabel(r'v ('+str(v1.unit)+')')
plt.axvline(x =obs_times[0].mjd, color = 'k', label = times[0])
plt.axvline(x= obs_times[1].mjd, color = 'r', label = times[1])
plt.plot(utc_times,v1,label='NS');
plt.plot(utc_times,v2,label= 'Comp');
#plt.plot(utc_difs[quad_points], v1[quad_points], marker = '*', color ='k', linestyle = 'None')
plot_datapoints()
plt.legend();
plt.title('e= '+ str(e)+ ',  a= '+str(a_thing.to(u.au))+ ',  $m_1$= '+str(m1.to(u.Msun))+',  $m_2$= '+str(m2.to(u.Msun)))
plt.show()

def calc_lam_obs(velocity, lam_rest):
    return (velocity/(u.cds.c.si)*lam_rest+lam_rest).to(u.angstrom)

lam_obs= calc_lam_obs(v2,lam_rest)
print ("rest wavelength:", lam_rest)
print ("Maximum wavelength observed (neglecting Earth-induced shifts and systemic velocity): ", np.nanmax(lam_obs))
plt.xlabel(r't (hours)')
plt.ylabel(r'$\lambda_{obs}$ ('+str(lam_obs.unit)+')')
plt.axvline(x =( (obs_times[0].mjd-obs_times[0].mjd)*u.day).to(u.hour).value, color = 'k', label = times[0])
plt.axvline(x=( (obs_times[1].mjd-obs_times[0].mjd)*u.day).to(u.hour).value, color = 'r', label = times[1])
plt.plot(utc_difs,lam_obs,label='Comp');
plt.legend();
plt.title('e= '+ str(e)+ ',  a= '+str(a_thing.to(u.au))+ ',  $m_1$= '+str(m1.to(u.Msun))+',  $m_2$= '+str(m2.to(u.Msun))+ r' $\lambda_{rest}$ = ' + str(lam_rest))
plt.show()


lam_range = np.linspace(3000., 7000, 8000)*u.angstrom
all_lams = calc_lam_obs(np.max(v2), lam_range)
plt.xlabel(r'$\lambda_{rest}$ ('+ str(lam_range.unit)+')')
plt.ylabel(r'$\lambda_{obs}$ ('+str(all_lams.unit)+') at max velocity of'+ str(np.max(v2)))
plt.plot(lam_range, all_lams, label = 'wavelengths');
plt.legend();
plt.title('e= '+ str(e)+ ',  a= '+str(a_thing.to(u.au))+ ',  $m_1$= '+str(m1.to(u.Msun))+',  $m_2$= '+str(m2.to(u.Msun))+ r'$v_{companion}$' + str(np.max(v2) ))
plt.show()
