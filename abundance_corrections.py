"""
Created by Ben Kaiser (UNC Chapel Hill) 2020-01-19

This should contain all of the functions that perform the diffusion corrections on the abundances

All of the logg values that are allowed as kwargs are required to be model grid points, meaning 7.5, 8.0, 8.5, or 9.0. I was worried about trying to interpolate because sodium starts doing weird stuff at the high gravities at the low temperatures of the model grid.
"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
import periodictable as pt
from astropy import units as u
from astropy import constants as const

import interp_tau as itau
import cal_params as cp

def declining_phase(target_teff, log_el1_over_el2, time,el1, el2, logg=8.0, steady_state_start=False, cross_extrap=True):
    """
    provide log_el1_over_el2 as an absolute number ratio not the one normalized to solar abundances.
    
    the output value will be in log10 of absolute number abundances... hopefully
    """
    if cross_extrap:
        el1_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el1)
        el2_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el2)
    else:
        el1_logtau= itau.extrapolate_single_el_tau(target_teff, el1, input_logg=logg)
        el2_logtau= itau.extrapolate_single_el_tau(target_teff, el2, input_logg=logg)
    
    el1_tau=10.**el1_logtau
    el2_tau=10.**el2_logtau
    exp_term= np.exp(time*((el2_tau-el1_tau)/(el1_tau*el2_tau))) #from equation 3 of Harrison et al. 2018
    if steady_state_start:
        coeff= 10.**(log_el1_over_el2) * el2_tau/el1_tau
    else:
        coeff= 10.**(log_el1_over_el2)
    dp_el1el2= coeff*exp_term
    return np.log10(dp_el1el2)


def get_time_since_accretion(target_teff, log_atm_ratio, log_desired_ratio, el1, el2,logg=8.0, steady_state_start=False,cross_extrap=True):
    """
    track back some abundance ratio from the present atmospheric abundance ratio to some expected abundance 
    ratio based on a solar system object (most likely). This gives the amount of time that must have passed since 
    accreting whatever body of that abundance in order to get the present day abundance.
    
    """
    if cross_extrap:
        print('cross_extrap')
        el1_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el1)
        el2_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el2)
    else:
        el1_logtau= itau.extrapolate_single_el_tau(target_teff, el1, input_logg=logg)
        el2_logtau= itau.extrapolate_single_el_tau(target_teff, el2, input_logg=logg)
    el1_tau=10.**el1_logtau
    el2_tau=10.**el2_logtau
    t_coeffs= (el1_tau*el2_tau)/(el2_tau-el1_tau)
    pollutant_term=log_desired_ratio- log_atm_ratio
    if steady_state_start:
        print('\nusing steady state\n')
        pollutant_term=pollutant_term+el1_logtau-el2_logtau
    else:
        pass
    pollutant_term=pollutant_term*np.log(10.)
    #print("pollutant_term", pollutant_term)
    #print("t_coeffs", t_coeffs)
    return t_coeffs*pollutant_term


def LiCa_DP_NaCa(target_teff, log_LiCa, log_NaCa, desired_log_NaCa, logg=8.0, steady_state_start=False, cross_extrap=True):
    """
    Take the log(Na/Ca ) in the atmosphere and an expected log(Na/Ca) for some sort of solar system object (most 
    likely) (desired_log_NaCa), and then using the time for since accretion from the function 
    get_time_since_accretion(), which is called internally, then un-decline the Li/Ca abundance to what the 
    log(Li/Ca ) would have been for the body that was accreted.
    
    """
    t_NaCa= get_time_since_accretion(target_teff, log_NaCa, desired_log_NaCa, "Na", "Ca",  logg=logg, steady_state_start=steady_state_start, cross_extrap=cross_extrap)
    dp_LiCa= declining_phase(target_teff, log_LiCa,t_NaCa,  'Li', "Ca", logg=logg, steady_state_start=steady_state_start, cross_extrap=cross_extrap)
    
    return t_NaCa, dp_LiCa


def recover_lost_element_number(t_passed, log_elHe_atm, log_m_cvz, el, el_tau):
    """
    Get back how much of whatever element must have fallen out of the convection zone for the given amount of 
    time that has passed since accretion
    
    """
    return

def get_accreted_mass( el, log_elHe,t_passed, teff=5000., logg=8.0,log_q=-5.0, m_wd=0.56, cross_extrap=True):
    """
    Assumes all of the convective mass can be treated as being helium and also that the helium isotope abundances are the same as that found on Earth(?, whatever the default periodic table mean molecular weight is). 
    
    returns the mass of whatever element that was present in the photosphere at the beginning of decline. 
    
    This should then be taken and divided by the mass of that element as a fraction of something like a chondrite 
    or eucrite or whatever to get the total mass of the accreted body.
    """
    el_num=cp.el_nums[el]
    if cross_extrap:
        log_el_tau=itau.extrapolate_tau_x_logg(teff, logg, el)
    else:
        log_el_tau=itau.extrapolate_single_el_tau(teff, el, input_logg=logg)
    #print('log_el_tau', log_el_tau)
    #print("t_passed", t_passed)
    #print('(t_passed/(10.**log_el_tau))',(t_passed/(10.**log_el_tau)))
    m_wd= (m_wd*const.M_sun).to(u.kg).value #converting the white dwarf mass to kg but making it a float again
    #print("np.log10(pt.elements[el_num].mass/pt.elements[2].mass)",np.log10(pt.elements[el_num].mass/pt.elements[2].mass))
    #print("pt.elements[el_num].mass", pt.elements[el_num].mass)
    #print("pt.elements[2].mass", pt.elements[2].mass)
    #print("log_elHe", log_elHe)
    #print("log_q", log_q)
    #print("np.log10(m_wd)", np.log10(m_wd))
    #print('====')
    log_m_acc= np.log10(pt.elements[el_num].mass/pt.elements[2].mass)+log_elHe+log_q+np.log10(m_wd)+(t_passed/(10.**log_el_tau))*np.log10(np.e)
    return log_m_acc



if __name__ == '__main__':
    print(cp.el_nums)
    time_range=np.linspace(0, 40, 20)
    time_range=time_range*1e6

    dp_lica= declining_phase(3830., 1.7,time_range, 'Li', 'Ca')
    print('DP log(Li/Ca)', dp_lica)

    #plt.plot(np.log10(time_range), dp_lica)
    #plt.plot(time_range, dp_lica)
    #plt.plot(time_range, declining_phase(3830., 1.7,time_range, 'Li', 'Ca', steady_state_start=True), label='Steady state start')
    #plt.legend()
    #plt.show()

    t_NaCa= get_time_since_accretion(3830., 0.0, -1.1, 'Na', 'Ca', steady_state_start=False, cross_extrap=True, logg=7.77)
    print("t_NaCa", t_NaCa, np.log10(t_NaCa))
    print(declining_phase(3830., 1.7,t_NaCa, 'Li', 'Ca', steady_state_start=False))

    #t_NaCa= get_time_since_accretion(3830., 0.0, -1.0, 'Na', 'Ca', logg=7.5)
    print("t_NaCa", t_NaCa, np.log10(t_NaCa))
    print(declining_phase(3830., 1.7,t_NaCa, 'Li', 'Ca', steady_state_start=False, logg=7.5))
    log_mCa= get_accreted_mass('Ca', -9.5, t_NaCa,  teff=3830., logg=8.0, log_q=-4.88, m_wd=0.56)
    print(log_mCa)
    mCa=10.**log_mCa
    print(mCa)
    m_total=mCa/0.073
    print('total accreted mass', m_total, 'kg')
    log_mNa= get_accreted_mass('Na', -9.5, t_NaCa,  teff=3830., logg=8.0, log_q=-4.88, m_wd=0.56)
    print(log_mNa)
    mNa=10.**log_mNa
    m_total_Na=mNa/0.0033
    print('total accreted mass', m_total_Na, 'kg')
    print("mNa/mCa", mNa/mCa)
    print(np.log10(pt.elements[cp.el_nums['Na']].mass/pt.elements[cp.el_nums['Ca']].mass*mNa/mCa))
