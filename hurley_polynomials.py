"""
Created by Ben Kaiser (UNC-Chapel Hill) 2020-01-06

This is essentially a Python implementation of some of the routines/equations presented in Hurley, Pols, and 
Tout 2000 with an emphasis (and probably entirety of implementation) focused on being used for white dwarf 
progenitor MS lifetimes. 

Actually not so much just the main-sequence lifetimes; I'm actually focusing on the total time spent by the 
progenitor before being a white dwarf, so that would include the time to giant branch, the time as a giant, and 
then also the core-helium burning time. I'm most likely still going to skip the AGB contribution since it's 
comparatively brief according to the Hurley et al. 2000 paper (assuming I read it correctly). 


"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy import constants as const
from astropy.time import Time
from astropy.table import Table
import scipy.interpolate as scinterp



import cal_params as cp
from a_coeffs import a_coeffs
from b_coeffs import b_coeffs

#nothing

#default_z=0.0001 #allegedly thick disk value
default_z=0.02 #approximately solar


def make_match(value, array):
    try:
        output_val=np.ones(array.shape)*value
    except AttributeError:
        output_val=value
    return output_val

def zeta(z):
    return np.log10(z/0.02)

def a(z,index):
    """
    the a_n equation from the beginning of Appendix A of Hurley et al. 2000
    
    I just left out the mu term because none of these variables have a mu value, and I didn't bother putting it in the
    a_coeffs list either as a result.
    """
    zeta_val=zeta(z)
    a_string="a"+str(index)
    def a_func(a_index):
        return a_index[0]+a_index[1]*zeta_val+a_index[2]*zeta_val**2+a_index[3]*zeta_val**3
    try:
        a_index=a_coeffs[a_string]
        return a_func(a_index)
    except KeyError:
        if index== 29:
            a_string="a'29"
            a_index=a_coeffs[a_string]
            return a_func(a_index)**a_func(a_coeffs["a32"])
        else:
            print("\n\n************")
            print("\n\nWARNING\n\n")
            print("a with index"+str(index) +" does not appear to exist in the a_coeffs dict and it isn't currently coded as an exception in a() in hurley_polynomials.py. You need to go look at a_coeffs.py and hurley_polynomials.py")
            print("\n\n************")
            return
    #a_out=a_index[0]+a_index[1]*zeta_val+a_index[2]*zeta_val**2+a_index[3]*zeta_val**3
    #print('a_'+str(index),a_out)
    return a_out

def b(z,index):
    b_string="b"+str(index)
    zeta_val=zeta(z)
    def b_func(b_row):
        return b_row[0]+b_row[1]*zeta_val+b_row[2]*zeta_val**2+b_row[3]*zeta_val**3
    try:
        b_row=b_coeffs[b_string]
        b_val=b_func(b_row)
        return b_val
    except KeyError:
        if index==40:
            return max(b_func(b_coeffs["b40_static"]),1.0)
        elif index==41:
            bprime=b_func(b_coeffs["b'41"])
            return bprime**(b_func(b_coeffs["b42"]))
        elif index==44:
            bprime= b_func(b_coeffs["b'44"])
            return bprime**5
        elif index==11:
            bprime=b_func(b_coeffs["b'11"])
            return bprime**2
        elif index==13:
            bprime=b_func(b_coeffs["b'13"])
            return bprime**2
        else:
            print("\n\n************")
            print("\n\nWARNING\n\n")
            print("b with index"+str(index) +" does not appear to exist in the b_coeffs dict and it isn't currently coded as an exception in b() in hurley_polynomials.py. You need to go look at b_coeffs.py and hurley_polynomials.py")
            print("\n\n************")
            return

def get_mass_HeF(z):
    """
    eq(2) from Hurley et al. 2000
    
    Input "z", which is 'metallicity' but actually the metal mass fraction of the star.
    
    """
    
    return 1.995+0.25*zeta(z)+0.087*zeta(z)**2

def get_t_bgb(z,mass):
    #print('a_1',a(z,1))
    #print('a(z,2)*mass**4',a(z,2)*mass**4)
    #print('a(z,3)*mass**5.5 ',a(z,3)*mass**5.5 )
    #print('mass**7',mass**7)
    #print('(a(z,1)+a(z,2)*mass**4+a(z,3)*mass**5.5 + mass**7)',(a(z,1)+a(z,2)*mass**4+a(z,3)*mass**5.5 + mass**7))
    #print('a(z,4)*mass**2',a(z,4)*mass**2)
    #print('a(z,5)*mass**7',a(z,5)*mass**7)
    #print('(a(z,4)*mass**2 + a(z,5)*mass**7)',(a(z,4)*mass**2 + a(z,5)*mass**7))
    return(a(z,1)+a(z,2)*mass**4+a(z,3)*mass**5.5 + mass**7)/(a(z,4)*mass**2 + a(z,5)*mass**7)

def get_x(z):
    first_val=make_match(0.95,z)
    last_val=make_match(0.99, z)
    return np.max([first_val,
                   np.min([
                       0.95-0.03*(zeta(z)+0.30103),
                       last_val
                       ],axis=0)],axis=0)

def get_mu(z,mass):
    #print(np.max([a(z,6)/(mass**a(z,7)),a(z,8)+(a(z,9)/mass**(a(z,10)))],axis=0))
    first_val=make_match(0.5, mass)
    return np.max([
        first_val,
        1.0-0.01 * np.max([
            a(z,6)/(mass**a(z,7)),
            a(z,8)+(a(z,9)/(mass**a(z,10)))
            ],axis=0)],axis=0)


def get_t_ms(mass, z=default_z):
    """
    Default is what I'm pretty sure is solar metallicity
    """
    t_bgb=get_t_bgb(z,mass)
    return np.max([
        get_mu(z,mass)*t_bgb,
        get_x(z)*t_bgb
        ],axis=0)*1e-3


def get_L_bgb(mass, z=default_z):
    """
    eq (10)
    
    """
    c2=9.301992
    c3=4.637345
    numerator=a(z,27)*mass**a(z,31)+a(z,28)*mass**c2
    denominator=a(z,29)+a(z,30)*mass**c3+mass**a(z,32)
    return numerator/denominator





def get_t_hems(mass):
    """
    Eq. (79) from Hurley et al. 2000
    """
    return (0.4129 + 18.81 * mass**4+1.853*mass**6)/mass**6.5

def get_t_he(mass,z):
    """
    equation (57) from Hurley et al. 2000
    
    
    Also obtains t_BGB as an intermediate step, but I'll have it return that as well I think since it's already computed.
    
    Returns output_t_he, t_bgb
    
    so it's a tuple, so you should receive it with 2 variables that will then each be arrays.
    """
    mass_HeF=get_mass_HeF(z)
    #print("Z=",z)
    #print('mass_HeF', mass_HeF)
    zeta_val=zeta(z)
    t_bgb=get_t_bgb(z,mass)
    
    def get_B(mass):
        """
        un-numbered equation below equation(38)
        """
        #print('mass.shape', mass.shape)
        first_entry=np.ones(mass.shape)*3e4
        the_array=np.array([first_entry, 500+1.75e4 * mass**0.6])
        #print("the_array.shape",the_array.shape)
        #return np.max(np.array([3e4, 500+1.75e4 * mass**0.6]), axis=0)
        return np.max(the_array, axis=0)


    def get_D(mass):
        """
        
        un-numbered equation above equation (39) 
        
        I'm pretty sure I need to add in a linear interpolation element to cover the gap (if there is one) between mass_HeF and 2.5
        """
        D_lo=5.37+0.135*zeta_val
        def get_D_hi(mass):
            try:
                the_array=np.array([-1.0*np.ones(mass.shape),0.975*D_lo-0.18*mass, 0.5*D_lo-0.06*mass])
            except AttributeError:
                the_array=np.array([-1.0,0.975*D_lo-0.18*mass, 0.5*D_lo-0.06*mass])
            return np.max(the_array,axis=0)
        D_hi= get_D_hi(mass)
        #Need D intermediate to cover values between mass_HeF and 2.5, which are not the same value, so I do need the linear interpolation over this range.
        mid_inds=np.where((mass > mass_HeF) & (mass < 2.5))
        slope=(get_D_hi(2.5)-D_lo)/(2.5-mass_HeF)
        D_mid= slope*(mass[mid_inds]-mass_HeF)+D_lo
        lo_inds=np.where(mass <= mass_HeF)
        hi_inds= np.where(mass >= 2.5)
        D_hi[lo_inds]=D_lo
        D_hi[mid_inds]=D_mid
        return 10**D_hi
    
    
    def get_m_x(mass):
        """
        equation (38). Exponent can be set to 1/3 because the difference of p and q is a constant 3 actually.
        """
        print(get_B(mass).shape, 'B shape')
        return (get_B(mass)/get_D(mass))**(1./3)
    def get_p(mass):
        """
        un-numbered equation below eq 38
        
        only the mass <= mass_HeF should matter as of 2020-04-12 because the t_he timescales only care about p for that mass range.
        """
        p_array=np.ones(mass.shape)
        hi_inds=np.where(mass >= 2.5)
        lo_inds=np.where(mass <= mass_HeF)
        p_array[hi_inds]=5.
        p_array[lo_inds]=6.
        return p_array
    
    def get_q(mass):
        """
        un-numbered equation below eq 38
        
        
        
        """
        q_array=np.ones(mass.shape)
        hi_inds= np.where(mass >= 2.5)
        lo_inds= np.where(mass <= mass_HeF)
        q_array[hi_inds]=2.
        q_array[lo_inds]=3.
        return q_array
    
    def get_A_H_prime(mass):
        """
        un-numbered equation above (44) 
        
        """
        first_entry=np.ones(mass.shape)*-4.8
        the_array=np.array([first_entry, np.min([-5.7+0.8*mass, -4.1+0.14*mass],axis=0)])
        return 10.**np.max(the_array, axis=0)
    
    def get_L_HeI(mass):
        def get_hi_L_HeI(mass):
            numerator= b(z,11)+b(z,12)*mass**3.8
            denominator=b(z,13)+mass**2
            return numerator/denominator
        alpha_1= (b(z,9)*mass_HeF**b(z,10)-get_hi_L_HeI(mass_HeF))/(get_hi_L_HeI(mass_HeF))
        def get_lo_L_HeI(mass):
            numerator=b(z,9)*mass**b(z,10)
            denominator=1+alpha_1*np.exp(15.*(mass-mass_HeF))
            return numerator/denominator
        hi_inds= np.where(mass >= mass_HeF)
        lo_inds=np.where(mass < mass_HeF)
        L_array= np.ones(mass.shape)
        L_array[hi_inds]=get_hi_L_HeI(mass[hi_inds])
        L_array[lo_inds]=get_lo_L_HeI(mass[lo_inds])
        return L_array
    
    p= get_p(mass) #I want this to be defined as a local variable that is "global" within the function so I don't have to call it in all of the functions that are coming up.
    q= get_q(mass)
    D = get_D(mass) # I decided I also want this essentially indefinitely defined.
    B= get_B(mass)
    A_H_prime=get_A_H_prime(mass)
    #print('A_H_prime:', A_H_prime)
    #L_bgb=1.
    L_bgb=get_L_bgb(mass)
    #L_x= 1.
    L_HeI=get_L_HeI(mass)
    m_x=get_m_x(mass)
    
    
    def get_core_luminosity_rel(mass_core):
        """
        equation (37)
        """
        return np.min([B*mass_core**q, D*mass_core**p], axis=0)
    
    L_x= get_core_luminosity_rel(m_x)
    
    def get_t_inf1(mass=mass):
        """
        Equation (40) of Hurley et al. 2000
        
        """
        first_term=1/((p-1)*A_H_prime*D)
        second_term= (D/L_bgb)**((p-1.)/p)
        return t_bgb+first_term*second_term
    
    def get_t_x(mass=mass):
        
        """
        equation (41)
        
        """
        
        t_inf1= get_t_inf1(mass=mass)
        first_term=t_inf1-t_bgb
        second_term=(L_bgb/L_x)**((p-1.)/p)
        
        return t_inf1-(first_term)*second_term
    
    def get_t_inf2(mass=mass):
        """
        equation (42)
        
        """
        t_x=get_t_x(mass=mass)
        first_term=1./((q-1)*A_H_prime*B)
        second_term=(B/L_x)**((q-1.)/q)
        return t_x + first_term*second_term
    
    def get_t_HeI(mass):
        """
        equation (43)
        """
        #L_HeI= np.ones(mass.shape) # *****to be replaced with the actual formula in the future
        t_array=np.ones(mass.shape)
        def get_t_lo():
            first_term= 1./((p-1)*A_H_prime*D)
            second_term=(D/L_HeI)**((p-1)/p)
            return get_t_inf1()-first_term*second_term
        def get_t_hi():
            first_term=1/((q-1)*A_H_prime*B)
            second_term=(B/L_HeI)**((q-1)/q)
            
            return get_t_inf2()-first_term*second_term
        hi_t_vals= get_t_hi()
        lo_t_vals=get_t_lo()
        hi_inds= np.where(L_HeI > L_x)
        lo_inds= np.where(L_HeI <= L_x)
        t_array[hi_inds]= hi_t_vals[hi_inds]
        t_array[lo_inds]=lo_t_vals[lo_inds]
        return t_array
    
    def get_mass_core_gb(time):
        """
        Equation (39) from Hurley et al. 2000
        
        """
        def get_lo_core_mass(time):
            
            return ((p-1)*A_H_prime*D*(get_t_inf1()-time))**(1./(1-p))
        
        def get_hi_core_mass(time):
            
            return ((q-1)*A_H_prime*B*(get_t_inf2() -time))**(1./(1-q))
        
        t_x= get_t_x()
        hi_inds= np.where(time > t_x)
        lo_inds= np.where(time <= t_x)
        t_array=np.ones(time.shape)
        t_array[lo_inds]=get_lo_core_mass(time)[lo_inds]
        t_array[hi_inds]=get_hi_core_mass(time)[hi_inds]
        
        return t_array
    
    def get_mass_core_alt(mass, t):
        """
        eq (34)
        
        """
        
        t_inf=t_bgb+1/(A_H_prime*D*(p-1))*(D/L_bgb)**((p-1)/p)
        
        return ((p-1)*A_H_prime*D*(t_inf-t))**(1/(1-p))
    
    ########
    
    t_HeI= get_t_HeI(mass)
    mass_core= get_mass_core_gb(t_HeI)
    alt_HeI_mass_core= get_mass_core_alt(mass, t_HeI)
    alt_bgb_mass_core = get_mass_core_alt(mass, t_bgb)
    
    #print('\n\n')
    #print('np.log10(t_bgb):', np.log10(t_bgb))
    #print('(t_HeI-t_bgb)/t_bgb', (t_HeI-t_bgb)/t_bgb)
    #print("t_HeI", t_HeI)
    #print('t_inf1', get_t_inf1())
    #print('t_x', get_t_x())
    #print('t_inf2', get_t_inf2())
    #print('L_bgb', L_bgb)
    #print('L_HeI', L_HeI, 'np.log10(L_HeI)', np.log10(L_HeI))
    #print('L_x', L_x)
    #print("D", D)
    #print('b(z,39)', b(z, 39))
    #print('b(z,40)', b(z,40))
    #print('\n')
    #print('mass_core (from eq 39)', mass_core)
    #print('alt_HeI_mass_core (from eq34 with t_HeI)', alt_HeI_mass_core)
    #print("alt_bgb_mass_core (from eq34 with t_bgb", alt_bgb_mass_core)
    
    #print('\n\n')
    
    #########3
    
    def hi_t_he(mass, t_bgb):
        """
        high-mass part of the piecewise function in equation (57)
        """
        numerator=t_bgb*b(z,41)*mass**b(z,42)+b(z,43)*mass**5
        denominator=b(z,44)+mass**5
        return numerator/denominator
    
    def lo_t_he(mass,mass_core, t_bgb):
        """
        low-mass part of the piecewise function in equation (57)
        """
        alpha4=(hi_t_he(mass_HeF, t_bgb)-b(z,39))/b(z,39)
        mu=(mass-mass_core)/(mass_HeF-mass_core)
        first_term=b(z,39)+(get_t_hems(mass_core)-b(z,39))*(1-mu)**b(z,40)
        second_term=1+alpha4*np.exp(15*(mass-mass_HeF))
        return first_term*second_term
    
    #m_x=1
    
    #print('hi_t_he(mass_HeF)', hi_t_he(mass_HeF, get_t_bgb(z, mass_HeF)))
    
    output_t_he=np.ones(mass.shape)
    lo_inds=np.where(mass< mass_HeF)
    hi_inds=np.where(mass >= mass_HeF)
    #hi_times= hi_t_he(mass[hi_inds], t_bgb[hi_inds])
    #lo_times= lo_t_he(mass[lo_inds],mass_core[lo_inds],t_bgb[lo_inds])
    hi_times= hi_t_he(mass, t_bgb)[hi_inds]
    lo_times= lo_t_he(mass,mass_core,t_bgb)[lo_inds]
    output_t_he[lo_inds]=lo_times
    output_t_he[hi_inds]=hi_times
    
    #plt.scatter(mass, A_H_prime)
    #plt.xlabel('Mass')
    #plt.ylabel('A_H_prime')
    #plt.show()
    
    #plt.scatter(mass, D)
    #plt.xlabel('mass')
    #plt.ylabel('D')
    #plt.axvline(mass_HeF, linestyle='--', color='k', label='M_HeF')
    #plt.axvline(2.5, label='M=2.5', color='r', linestyle='--')
    #plt.legend()
    #plt.show()
    
    #plt.plot(mass, mass_core, label='Core Mass')
    #plt.plot(mass, m_x, label='M_x')
    #plt.legend()
    #plt.xlabel('Mass')
    #plt.show()
    
    #plt.scatter(mass, np.log10(L_HeI))
    #plt.xlabel('Mass')
    #plt.ylabel('L_HeI')
    #plt.show()
    
    #plt.plot(mass, mass_core, label='eq 39')
    #plt.plot(mass, alt_HeI_mass_core, label='eq34 t_HeI')
    #plt.plot(mass, alt_bgb_mass_core, label='eq 34 t_bgb')
    #plt.xlabel('Mass')
    #plt.ylabel('Core Mass')
    #plt.legend()
    #plt.show()
    
    #plt.plot(mass, mass_core/mass, label='eq 39')
    #plt.plot(mass, alt_HeI_mass_core/mass, label='eq34 t_HeI')
    #plt.plot(mass, alt_bgb_mass_core/mass, label='eq 34 t_bgb')
    #plt.xlabel('Mass')
    #plt.ylabel('Core Mass/Mass')
    #plt.legend()
    #plt.show()
    return (output_t_he+t_bgb)*1e-3
    #return output_t_he, t_bgb
#test_mass=1.
#print('zeta', zeta(default_z))
#print('t_bgb', get_t_bgb(default_z, test_mass))
#print('mu', get_mu(default_z, test_mass))
#print('x', get_x(default_z))
#print('t_ms',get_t_ms(test_mass))


#mass_vals=np.linspace(0.1,13.,100)
#log_mass_vals=np.linspace(-0.4, 1.9, 100)
#mass_vals= 10.**log_mass_vals

#t_ms_vals= get_t_ms(mass_vals)
#t_bgb_vals=get_t_bgb(default_z, mass_vals)


#plt.plot(np.log10(mass_vals), np.log10(t_ms_vals))
##plt.plot(np.log10(mass_vals), np.log10(t_bgb_vals))
#plt.plot(np.log10(mass_vals), np.log10(get_t_ms(mass_vals,z=0.0001)))
#plt.xlim(-0.4, 1.9)
#plt.show()



if __name__ == '__main__':
    #test_masses=np.random.normal(2,0.2, 1000)
    test_masses=np.linspace(0.8,10,1000)
    #test_masses=np.array([1.0,2.0,5.0])
    test_t_he, test_t_bgb= get_t_he(test_masses, default_z)
    print('average t_he:', np.mean(test_t_he), 'average t_bgb:', np.mean(test_t_bgb))
    #print('average m_x:', np.mean(test_m_x))
    
    #plt.scatter(test_masses, test_m_x)
    #plt.xlabel('Mass')
    #plt.ylabel(r'$M_x$')
    #plt.axvline(get_mass_HeF(default_z), linestyle='--', color='k', label='M_HeF')
    #plt.axvline(2.5, label='M=2.5', color='r', linestyle='--')
    #plt.legend()
    #plt.show()
    
    plt.scatter(test_masses, test_t_he)
    plt.xlabel('Mass')
    plt.ylabel('t_he (Myr)')
    plt.show()
    
    
    plt.scatter(test_masses, test_t_bgb)
    plt.xlabel('Mass')
    plt.ylabel('t_bgb (Myr)')
    plt.show()
    
    plt.scatter(test_t_bgb, test_t_he)
    plt.ylabel('t_he(Myr)')
    plt.xlabel('t_bgb (Myr)')
    plt.show()
    
    plt.scatter(test_masses, test_t_he/test_t_bgb)
    plt.xlabel('Mass')
    plt.ylabel('t_he/t_bgb')
    plt.axvline(get_mass_HeF(default_z), linestyle='--', color='k')
    plt.show()
    
    plt.hist(test_masses)
    plt.title('Masses')
    plt.show()
    
    plt.title('t_he')
    plt.hist(test_t_he)
    plt.show()
    
    plt.title('t_bgb')
    plt.hist(test_t_bgb)
    plt.show()
    
