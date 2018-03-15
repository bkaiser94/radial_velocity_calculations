import numpy as np
import os
from glob import glob
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
from astropy.units import cds
cds.enable()

ra = sys.argv[1]
dec = sys.argv[2]
times = ['2018-03-19T02:19:00', '2018-03-19T10:01:00']
#obs_times= Time(times, format = 'isot', scale='utc', location = coord.EarthLocation.of_site('Cerro Tololo'))
parkes_location = coord.EarthLocation.from_geocentric(x = -4554231.533*u.m,y= 2816759.109*u.m, z =  -3454036.323*u.m) # from http://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html 
cerro_pachon_location = coord.EarthLocation.from_geodetic(lat = '-30 14 16.41', lon = '-70 44 01.11', height = 2748* u.m)
m1 = (0.14 *u.Msun).si
m2 = (1.4 *u.Msun).si
e = 0.000023
omega = 97 * u.degree
period = (0.4497391377 * u.day).to(u.second)
t0 = Time(55756.23, format = 'mjd', scale= 'utc', location = parkes_location)

for thing in coord.EarthLocation.get_site_names():
    print (thing)
print (u.Unit('m'))
print (parkes_location)
print (period)
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

N =1000
th= np.linspace(0.,4*np.pi,N)
a_thing=get_a(period,m1,m2)
dt= (r(th,a_thing,e)**2/L(m1,m2,a_thing,e)*4.*np.pi/N).to(u.second)
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

plt.xlabel(r't (days)')
plt.ylabel(r'v ('+str(v1.unit)+')')
plt.plot(days,v1,label='Comp');
plt.plot(days,v2,label= 'NS');
plt.legend();
plt.title('e= '+ str(e)+ ',  a= '+str(a_thing.to(u.au))+ ',  $m_1$= '+str(m1.to(u.Msun))+',  $m_2$= '+str(m2.to(u.Msun)))
plt.show()
print(N)
print (p(a_thing,m1,m2).to(u.day))
print (v1.unit)
