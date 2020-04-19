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

n_points=1e5

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

def get_t_relHe_fwd(el, target_teff, log_elHe_atm, log_elHe_des, logg=8.0, cross_extrap=True):
    """
    Take the present-day value for log10(el/He) and figure out how long would have to pass for diffusion to lower
    the overall photospheric abundance to log_elHe_des
    
    """
    if cross_extrap:
        el_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el)
    else:
        el_logtau= itau.extrapolate_single_el_tau(target_teff, el, input_logg=logg)
    return (log_elHe_atm- log_elHe_des)*np.log(10.)*10.**(el_logtau)

def get_relHe_fwd(el, time,  target_teff, log_elHe_atm, logg=8.0, cross_extrap=True):
    """
    time input in Myr, it will be converted to years inside this function
    """
    if cross_extrap:
        el_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el)
    else:
        el_logtau= itau.extrapolate_single_el_tau(target_teff, el, input_logg=logg)
    time=time*1e6
    return log_elHe_atm+np.log10(np.e)*(-time/(10.**el_logtau))

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

def el1el2_DP_el3el2(target_teff, log_el1el2, log_el3el2, desired_log_el3el2, el1, el2, el3,  logg=8.0, steady_state_start=False, cross_extrap=True):
    """
    Take the log(el3/el2) in the atmosphere and an expected log(el3/el2) for some sort of solar system object (most likely) (desired_log_el3el2), and then using the time for since accretion from the function 
    get_time_since_accretion(), which is called internally, then un-decline the el1/el2 abundance to what the 
    log(el1 / el2) would have been for the body that was accreted.
    
    This is the more generalized form of LiCa_DP_NaCa(), hopefully.
    
    """
    t_el3el2= get_time_since_accretion(target_teff, log_el3el2, desired_log_el3el2, el3, el2,  logg=logg, steady_state_start=steady_state_start, cross_extrap=cross_extrap)
    dp_el1el2= declining_phase(target_teff, log_el1el2, t_el3el2,  el1,  el2, logg=logg, steady_state_start=steady_state_start, cross_extrap=cross_extrap)
    
    return t_el3el2, dp_el1el2

def get_el1el2_wrt_time(log_el1el2, time, el1_tau, el2_tau):
    """
    time in Myr
    
    Not for steady state start. This assumes declining from early phase abundance
    """
    return log_el1el2 + time*1e6*np.log10(np.e)*((10.**el2_tau-10.**el1_tau)/(10.**el1_tau * 10.**el2_tau))

def el1el2_DP_el3el2_ftimes(target_teff, log_el1el2, log_el3el2, time,  el1, el2, el3,  logg=8.0, steady_state_start=False, cross_extrap=True):
    if cross_extrap:
        el1_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el1)
        el2_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el2)
        el3_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el3)
    else:
        el1_logtau= itau.extrapolate_single_el_tau(target_teff, el1, input_logg=logg)
        el2_logtau= itau.extrapolate_single_el_tau(target_teff, el2, input_logg=logg)
        el3_logtau= itau.extrapolate_single_el_tau(target_teff, el3, input_logg=logg)
    print(el1, 'log tau', el1_logtau)
    print(el2, 'log tau', el2_logtau)
    print(el3, 'log tau', el3_logtau)
    
    log_el1el2_dp= get_el1el2_wrt_time(log_el1el2, time, el1_logtau, el2_logtau)
    log_el3el2_dp= get_el1el2_wrt_time(log_el3el2, time, el3_logtau, el2_logtau)
    return log_el1el2_dp, log_el3el2_dp

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


def easy_dist_decline(wd_row, el1, el2, el3, desired_log_el3el2, n_points=n_points, plot_all=False, start_he=False):
    teff_dist=np.random.normal(loc=wd_row['teff'], scale=wd_row['teff_err'], size=n_points)
    logg_dist= np.random.normal(loc=wd_row['logg'], scale=wd_row['logg_err'], size=n_points)
    def get_el(el, el2=el2):
        if start_he:
            el_stem= el.lower()+'/he'
        else:
            el_stem=el.lower()+'/'+el2.lower()
        el_abund= wd_row[el_stem]
        el_err= wd_row[el_stem+'_err']
        return el_abund, el_err
    def make_el_dist(el, el2=el2):
        el_abund, el_err= get_el(el, el2)
        el_dist= np.random.normal(loc=el_abund, scale= el_err, size=n_points)
        return el_dist
    if start_he:
        el1_dist= make_el_dist(el1)
        el2_dist= make_el_dist(el2)
        el3_dist= make_el_dist(el3)
        log_el1el2= el1_dist-el2_dist
        log_el3el2= el3_dist-el2_dist
    else:
        print('doing simultaneous decline')
        log_el1el2= make_el_dist(el1, el2=el2)
        log_el3el2=make_el_dist(el3,el2=el2)
    t_decline, dp_el1el2= el1el2_DP_el3el2(teff_dist, log_el1el2, log_el3el2, desired_log_el3el2, el1, el2, el3, logg=logg_dist, cross_extrap=True)
    
    dp_el1el2=dp_el1el2[~np.isinf(dp_el1el2)]
    dp_el1el2=dp_el1el2[~np.isnan(dp_el1el2)]
    
    t_decline=t_decline[~np.isinf(t_decline)]
    t_decline=t_decline[~np.isnan(t_decline)]
    
    print(wd_row['name'], 'mean log('+el1+'/'+el2+')' , np.mean(dp_el1el2), '+/-',np.std(dp_el1el2))
    print(wd_row['name'], 'median log('+el1+'/'+el2+')' , np.median(dp_el1el2), 'up/down',np.percentile(dp_el1el2, 84),np.percentile(dp_el1el2, 16))
    print('mean t_decline', np.mean(t_decline*1e-6),'Myr')
    if plot_all:
        plt.hist(dp_el1el2, bins=101, normed=True)
        plt.title(wd_row['name'])
        plt.show()
        plt.hist(t_decline*1e-6, bins=101, normed=True)
        plt.xlabel('t_decline (Myr)')
        plt.title(wd_row['name'])
        plt.show()
    else:
        pass
    return


if __name__ == '__main__':
    #wd_name='WDJ2356-209'
    target_logg=7.98
    target_logg_err=0.07
    target_teff= 4040. #K
    target_teff_err=110.
    
    logg_dist=np.random.normal(loc=target_logg, scale=target_logg_err, size=n_points)
    teff_dist= np.random.normal(loc=target_teff, scale=target_teff_err, size=n_points)
    #fe_dist= np.random.normal(loc=-8.6, scale=0.2, size=n_points)
    fe_dist= np.random.normal(loc=-11.7, scale=0.0001, size=n_points)
    ca_dist=np.random.normal(loc=-9.4, scale=0.2, size=n_points)
    na_dist= np.random.normal(loc=-8.3, scale=0.2, size=n_points)
    naca_line= np.linspace(0.25,1.25, 100)
    feca_dist= fe_dist-ca_dist
    naca_dist= na_dist-ca_dist
    #tFeCa, dp_FeCa= el1el2_DP_el3el2(teff_dist, 0.8, 1.1, -0.01, 'Fe', 'Ca', 'Na', logg=logg_dist, cross_extrap=True)
    
    #tFeCa, dp_FeCa= el1el2_DP_el3el2(target_teff,0.8, naca_line, -0.01, 'Fe', 'Ca', 'Na', logg=target_logg, cross_extrap=True)
    #plt.plot(naca_line, dp_FeCa)
    #plt.xlim(-2.0, 1.25)
    #plt.ylim(0.0, 2.0)
    #plt.show()
    tFeCa, dp_FeCa= el1el2_DP_el3el2(teff_dist,feca_dist, naca_dist, -0.01, 'Li', 'Ca', 'Na', logg=logg_dist, cross_extrap=True)
    #plt.scatter(teff_dist, dp_FeCa)
    #plt.show()
    print("mean", np.mean(dp_FeCa), np.std(dp_FeCa))
    plt.hist(dp_FeCa, bins=101, normed=True, label='log(Na/Ca)=-0.01', alpha=0.2)
    tFeCa, dp_FeCa= el1el2_DP_el3el2(teff_dist,feca_dist, naca_dist, -1.1, 'Li', 'Ca', 'Na', logg=logg_dist, cross_extrap=True)
    plt.hist(dp_FeCa, bins=101, normed=True, label='log(Na/Ca)=-1.1', alpha=0.2)
    print("mean", np.mean(dp_FeCa), np.std(dp_FeCa))
    plt.legend()
    plt.show()
    
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
