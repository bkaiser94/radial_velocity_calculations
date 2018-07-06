
import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt
#from astropy.io import fits
import sys
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
from astropy.units import cds
import astropy.constants as const
cds.enable()
#plt.rc('font', size =18)
#plt.rc('lines', markersize=12)
plt.rc('font', size = 12)
plt.rc('lines', markersize = 12)

precision  = 3

#P_orbital = 1.354885217 *u.day  #PSR J1435-6100
P_orbital = 0.4497391377 *u.day  #PSR J1431-4715

#q= 0.90/1.35 #inclination 90 
#q= 1.35/0.90 #PSR J1435-6100
q= 10.582 #PSR J1431-4715 from RV fits

def mass_function(m1, m2, inclination):
    return (m2*np.sin(inclination))**3 / (m1+m2)**2
#def get_m2(f, 
print const.G.to(u.Rsun**3/(u.Msun * u.day**2))
def get_roche_lobe(m_ns, P, q):
    #first_term = (m_ns.si * P.si**2 * const.G.si*(1+q)/(4*np.pi**2))**(1./3)
    first_term = (m_ns.si * P.si**2 * const.G.si*(1+1/q)/(4*np.pi**2))**(1./3)
    second_term = (0.49*q**(-2./3))/(0.6*q**(-2./3)+np.log(1+q**(-1./3)))
    return (first_term*second_term).to(u.Rsun)

def get_mean_density(period):
    """
    return the mean density in g cm^-3
    """
    return 107*(period.to(u.hour).value)**(-2)

m_ns_range = np.linspace(1,3.4, 1000)*u.Msun

roche_lobes = get_roche_lobe(m_ns_range, P_orbital, q)

period_range = (10**np.linspace(np.log10(0.01), np.log10(90), 1000))*u.day # exponentially spaced periods to use
print "Mean density: ", get_mean_density(P_orbital) , " g/cm^3"
#print 0.88/1.35
print 1/q

plt.plot(period_range.value, get_mean_density(period_range))
plt.xlabel('Period (Days)')
plt.ylabel('Density (g/cm^3)')
plt.xscale('log')
plt.yscale('log')
plt.show()

plt.plot(m_ns_range, roche_lobes)
plt.xlabel(r"$M_{NS} (M_{\odot})$")
plt.ylabel(r"$R_L(R_{\odot})$" )
plt.title("Roche Lobe vs. Neutron Star Mass for Period=" + str(np.round(P_orbital, precision)) + " , q=" + str(np.round(q,precision)) + r" , and $i= 90\degree$")

plt.show()

plt.plot(m_ns_range/q, roche_lobes)
plt.xlabel(r"$M_{comp} (M_{\odot})$")
plt.title("Roche Lobe vs. Companion Star Mass for Period=" + str(np.round(P_orbital, precision)) + " , q=" + str(np.round(q,precision)) + r" , and $i= 90\degree$")

plt.ylabel(r"$R_L(R_{\odot})$" )

plt.show()





