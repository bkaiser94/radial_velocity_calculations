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

ra = float(sys.argv[1]) #values in decimal degrees
dec = float(sys.argv[2])

parkes_location = coord.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
cerro_pachon_location = coord.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)

times = ['2018-03-19T02:19:00', '2018-03-19T10:11:00']
obs_times= Time(times, format = 'isot', scale='utc', location = cerro_pachon_location)
target_coord = coord.SkyCoord(ra, dec, unit= (u.deg, u.deg), frame= 'icrs')

def to_barycenter(input_times):
    bary_corr =input_times.light_travel_time(target_coord)
    return (input_times+ bary_corr).tdb.mjd


lam_rest = 6562.81*u.angstrom #angstroms
m1 = (0.12 *u.Msun).si
m2 = (1.4 *u.Msun).si
e = 0.000023
omega = 97 * u.degree
period = (0.4497391377 * u.day).to(u.second)
t0 = Time(55756.23, format = 'mjd', scale= 'utc', location = parkes_location)
tasc = Time(55756.1047771, format = 'mjd', scale= 'utc', location = parkes_location)


bmjd_obs = to_barycenter(obs_times) #corrected to barycenter to use against the rv curve
bmjd_t0 = to_barycenter(t0) #corrected initial epoch
bmjd_tasc= to_barycenter(tasc)

#calculate nearest bmjd to start for calculating radial velocities around the observations
time_dif= bmjd_obs[0]-bmjd_tasc #days between epoch and beginning of observability
nearest_time = int(time_dif/period.to(u.day).value)*period.to(u.day).value+bmjd_tasc #beginning of the orbital cycle
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
N= 3000
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
    
v1= v(th, a_thing, e, m1, m2, dt, 1)

v2= v(th, a_thing , e, m1, m2, dt, 2)

print ("Maximum radial velocity companion:", np.nanmax(v1))
bmjd_times =Time(nearest_time+ days.value,  format = 'mjd', scale='tdb', location = cerro_pachon_location)
mjd_times= bmjd_times+ bmjd_times.light_travel_time(target_coord)
utc_times = mjd_times.utc.mjd
utc_difs =((utc_times-obs_times[0].mjd)*u.day).to(u.hour)
plt.xlabel(r't (hours)')
plt.ylabel(r'v ('+str(v1.unit)+')')
plt.axvline(x =( (obs_times[0].mjd-obs_times[0].mjd)*u.day).to(u.hour).value, color = 'k', label = times[0])
plt.axvline(x=( (obs_times[1].mjd-obs_times[0].mjd)*u.day).to(u.hour).value, color = 'r', label = times[1])
plt.plot(utc_difs,v1,label='Comp');
plt.plot(utc_difs,v2,label= 'NS');
plt.legend();
plt.title('e= '+ str(e)+ ',  a= '+str(a_thing.to(u.au))+ ',  $m_1$= '+str(m1.to(u.Msun))+',  $m_2$= '+str(m2.to(u.Msun)))
plt.show()
print(N)
print (p(a_thing,m1,m2).to(u.day))
print (v1.unit)
print (a_thing.to(u.au))
print ("a in ls", (a_thing/(u.cds.c.si)).to(u.second))


def calc_lam_obs(velocity, lam_rest):
    return (velocity/(u.cds.c.si)*lam_rest+lam_rest).to(u.angstrom)

lam_obs= calc_lam_obs(v1,lam_rest)
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
all_lams = calc_lam_obs(np.max(v1), lam_range)
plt.xlabel(r'$\lambda_{rest}$ ('+ str(lam_range.unit)+')')
plt.ylabel(r'$\lambda_{obs}$ ('+str(all_lams.unit)+') at max velocity of'+ str(np.max(v1)))
plt.plot(lam_range, all_lams);
plt.legend();
plt.title('e= '+ str(e)+ ',  a= '+str(a_thing.to(u.au))+ ',  $m_1$= '+str(m1.to(u.Msun))+',  $m_2$= '+str(m2.to(u.Msun))+ r'$v_{companion}$' + str(np.max(v1) ))
plt.show()
