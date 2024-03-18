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

n_points=int(1e4)
limit_indicator=99. #value above which if the absolute value of the error on a measurement is above it indicates it should be a limit

#default values to populate definitions of kwargs

#default_modeler='Fontaine2015'
#default_atm_type='He'
#default_overshoot=0.0

#default_modeler='Koester2020'
#default_atm_type='H'
#default_overshoot=1.0

default_modeler='Koester2020'
default_atm_type='Nonsense' #This should make it crash if I am not correctly passing around kwargs
default_overshoot=1.0


default_cross_extrap=True

def get_el1el2_full_err(abund_table,el1, el2,n_sigma=1.):
    """
    Basically equation A2 of Klein et al. 2021 (the Beryllium paper), which combines the errors on each element abundance without double-counting the T_eff uncertainties.
    
    I'm going to have it receive two rows from an astropy table as the input so I don't have to specify a ton of specific variables
    
    """
    first_term=(abund_table[el1+'/he_spread_err']/abund_table[el1+'/he'])**2.
    second_term=(abund_table[el2+'/he_spread_err']/abund_table[el2+'/he'])**2.
    #third term=((abund_table[el1+'/he_teff_err']/abund_table[el1+'/he'])-(abund_table[el2+'/he_teff_err']/abund_table[el2+'/he_teff_err']))**2.
    third_term_first_half=abund_table[el1+'/he_teff_err']/abund_table[el1+'/he']
    third_term_second_half=abund_table[el2+'/he_teff_err']/abund_table[el2+'/he']
    third_term=(third_term_first_half-third_term_second_half)**2.
    el1el2=abund_table[el1+'/he']/abund_table[el2+'/he']
    full_error=np.sqrt(first_term+second_term+third_term)*el1el2
    print(el1+'/'+el2, el1el2, '+/-',full_error)
    #upper_bound=el1el2+full_error
    #lower_bound=el1el2-full_error
    upper_bound=el1el2+(full_error*n_sigma)
    lower_bound=el1el2-(full_error*n_sigma)
    log_hi_bound=np.log10(upper_bound)
    log_lo_bound=np.log10(lower_bound)
    print('log10('+ el1+'/'+el2+'):', np.log10(el1el2),',upper bound:',log_hi_bound, ',lower bound:', log_lo_bound)
    return np.log10(el1el2), log_lo_bound, log_hi_bound

def declining_phase(target_teff, log_el1_over_el2, time,el1, el2, logg=8.0, steady_state_start=False, cross_extrap=True, modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot):
    """
    provide log_el1_over_el2 as an absolute number ratio not the one normalized to solar abundances.
    
    the output value will be in log10 of absolute number abundances... hopefully
    """
    if cross_extrap:
        el1_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el1, modeler=modeler, atm_type=atm_type, overshoot=overshoot)
        el2_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el2, modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    else:
        el1_logtau= itau.extrapolate_single_el_tau(target_teff, el1, input_logg=logg, modeler=modeler, atm_type=atm_type, overshoot=overshoot)
        el2_logtau= itau.extrapolate_single_el_tau(target_teff, el2, input_logg=logg, modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    
    el1_tau=10.**el1_logtau
    el2_tau=10.**el2_logtau
    exp_term= np.exp(time*((el2_tau-el1_tau)/(el1_tau*el2_tau))) #from equation 3 of Harrison et al. 2018
    if steady_state_start:
        coeff= 10.**(log_el1_over_el2) * el2_tau/el1_tau
    else:
        coeff= 10.**(log_el1_over_el2)
    dp_el1el2= coeff*exp_term
    return np.log10(dp_el1el2)


def get_time_since_accretion(target_teff, log_atm_ratio, log_desired_ratio, el1, el2,logg=8.0, steady_state_start=False,cross_extrap=True, modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot):
    """
    track back some abundance ratio from the present atmospheric abundance ratio to some expected abundance 
    ratio based on a solar system object (most likely). This gives the amount of time that must have passed since 
    accreting whatever body of that abundance in order to get the present day abundance.
    
    """
    if cross_extrap:
        print('cross_extrap')
        el1_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el1, modeler=modeler, atm_type=atm_type, overshoot=overshoot)
        el2_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el2)
    else:
        el1_logtau= itau.extrapolate_single_el_tau(target_teff, el1, input_logg=logg, modeler=modeler, atm_type=atm_type, overshoot=overshoot)
        el2_logtau= itau.extrapolate_single_el_tau(target_teff, el2, input_logg=logg, modeler=modeler, atm_type=atm_type, overshoot=overshoot)
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

def get_t_relHe_fwd(el, target_teff, log_elHe_atm, log_elHe_des, logg=8.0, cross_extrap=True,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot):
    """
    Take the present-day value for log10(el/He) and figure out how long would have to pass for diffusion to lower
    the overall photospheric abundance to log_elHe_des
    
    """
    if cross_extrap:
        el_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    else:
        el_logtau= itau.extrapolate_single_el_tau(target_teff, el, input_logg=logg,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    return (log_elHe_atm- log_elHe_des)*np.log(10.)*10.**(el_logtau)

def get_relHe_fwd(el, time,  target_teff, log_elHe_atm, logg=8.0, cross_extrap=True,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot):
    """
    time input in Myr, it will be converted to years inside this function
    """
    if cross_extrap:
        el_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    else:
        el_logtau= itau.extrapolate_single_el_tau(target_teff, el, input_logg=logg,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    time=time*1e6
    return log_elHe_atm+np.log10(np.e)*(-time/(10.**el_logtau))

def LiCa_DP_NaCa(target_teff, log_LiCa, log_NaCa, desired_log_NaCa, logg=8.0, steady_state_start=False, cross_extrap=True,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot):
    """
    Take the log(Na/Ca ) in the atmosphere and an expected log(Na/Ca) for some sort of solar system object (most 
    likely) (desired_log_NaCa), and then using the time for since accretion from the function 
    get_time_since_accretion(), which is called internally, then un-decline the Li/Ca abundance to what the 
    log(Li/Ca ) would have been for the body that was accreted.
    
    """
    t_NaCa= get_time_since_accretion(target_teff, log_NaCa, desired_log_NaCa, "Na", "Ca",  logg=logg, steady_state_start=steady_state_start, cross_extrap=cross_extrap,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    dp_LiCa= declining_phase(target_teff, log_LiCa,t_NaCa,  'Li', "Ca", logg=logg, steady_state_start=steady_state_start, cross_extrap=cross_extrap,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    
    return t_NaCa, dp_LiCa

def el1el2_DP_el3el2(target_teff, log_el1el2, log_el3el2, desired_log_el3el2, el1, el2, el3,  logg=8.0, steady_state_start=False, cross_extrap=True,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot):
    """
    Take the log(el3/el2) in the atmosphere and an expected log(el3/el2) for some sort of solar system object (most likely) (desired_log_el3el2), and then using the time for since accretion from the function 
    get_time_since_accretion(), which is called internally, then un-decline the el1/el2 abundance to what the 
    log(el1 / el2) would have been for the body that was accreted.
    
    This is the more generalized form of LiCa_DP_NaCa(), hopefully.
    
    """
    t_el3el2= get_time_since_accretion(target_teff, log_el3el2, desired_log_el3el2, el3, el2,  logg=logg, steady_state_start=steady_state_start, cross_extrap=cross_extrap,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    dp_el1el2= declining_phase(target_teff, log_el1el2, t_el3el2,  el1,  el2, logg=logg, steady_state_start=steady_state_start, cross_extrap=cross_extrap,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    
    return t_el3el2, dp_el1el2

def get_el1el2_wrt_time(log_el1el2, time, el1_tau, el2_tau):
    """
    time in Myr
    
    Not for steady state start. This assumes declining from early phase abundance
    """
    return log_el1el2 + time*1e6*np.log10(np.e)*((10.**el2_tau-10.**el1_tau)/(10.**el1_tau * 10.**el2_tau))

def el1el2_DP_el3el2_ftimes(target_teff, log_el1el2, log_el3el2, time,  el1, el2, el3,  logg=8.0, steady_state_start=False, cross_extrap=True,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot):
    if cross_extrap:
        el1_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el1,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
        el2_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el2,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
        el3_logtau=itau.extrapolate_tau_x_logg(target_teff, logg, el3,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    else:
        el1_logtau= itau.extrapolate_single_el_tau(target_teff, el1, input_logg=logg,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
        el2_logtau= itau.extrapolate_single_el_tau(target_teff, el2, input_logg=logg,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
        el3_logtau= itau.extrapolate_single_el_tau(target_teff, el3, input_logg=logg,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    print(el1, 'log tau', el1_logtau,"{:e}".format(10.**el1_logtau))
    print(el2, 'log tau', el2_logtau,"{:e}".format(10.**el2_logtau))
    print(el3, 'log tau', el3_logtau,"{:e}".format(10.**el3_logtau))
    
    log_el1el2_dp= get_el1el2_wrt_time(log_el1el2, time, el1_logtau, el2_logtau)
    log_el3el2_dp= get_el1el2_wrt_time(log_el3el2, time, el3_logtau, el2_logtau)
    return log_el1el2_dp, log_el3el2_dp

def recover_lost_element_number(t_passed, log_elHe_atm, log_m_cvz, el, el_tau):
    """
    Get back how much of whatever element must have fallen out of the convection zone for the given amount of 
    time that has passed since accretion
    
    """
    return

def get_accreted_mass( el, log_elHe,t_passed, teff=5000., logg=8.0,log_q=-5.0, m_wd=0.56, cross_extrap=True,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot):
    """
    Assumes all of the convective mass can be treated as being helium and also that the helium isotope abundances are the same as that found on Earth(?, whatever the default periodic table mean molecular weight is). 
    
    returns the mass of whatever element that was present in the photosphere at the beginning of decline. 
    
    This should then be taken and divided by the mass of that element as a fraction of something like a chondrite 
    or eucrite or whatever to get the total mass of the accreted body.
    """
    el_num=cp.el_nums[el]
    if cross_extrap:
        log_el_tau=itau.extrapolate_tau_x_logg(teff, logg, el,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    else:
        log_el_tau=itau.extrapolate_single_el_tau(teff, el, input_logg=logg,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
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


def easy_dist_decline(wd_row, el1, el2, el3, desired_log_el3el2, n_points=n_points, plot_all=False, start_he=False,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot):
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
    t_decline, dp_el1el2= el1el2_DP_el3el2(teff_dist, log_el1el2, log_el3el2, desired_log_el3el2, el1, el2, el3, logg=logg_dist, cross_extrap=True,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    
    dp_el1el2=dp_el1el2[~np.isinf(dp_el1el2)]
    dp_el1el2=dp_el1el2[~np.isnan(dp_el1el2)]
    
    t_decline=t_decline[~np.isinf(t_decline)]
    t_decline=t_decline[~np.isnan(t_decline)]
    
    print(wd_row['name'], 'mean log('+el1+'/'+el2+')' , np.mean(dp_el1el2), '+/-',np.std(dp_el1el2))
    print(wd_row['name'], 'median log('+el1+'/'+el2+')' , np.median(dp_el1el2), 'up/down',np.percentile(dp_el1el2, 84),np.percentile(dp_el1el2, 16))
    print('mean t_decline', np.mean(t_decline*1e-6),'Myr')
    if plot_all:
        plt.hist(dp_el1el2, bins=101, density=True)
        plt.title(wd_row['name'])
        plt.show()
        plt.hist(t_decline*1e-6, bins=101, density=True)
        plt.xlabel('t_decline (Myr)')
        plt.title(wd_row['name'])
        plt.show()
    else:
        pass
    return

def easy_dist_ssp(wd_row,elements,n_points=n_points, plot_all=False, tau_rand=False, tau_add=0.2,modeler=default_modeler, overshoot=default_overshoot):
    atm_type=wd_row['diff_atm_type']
    string1= elements[0].lower()+'/'+elements[1].lower()
    string2=elements[2].lower()+'/'+elements[1].lower()
    #times= np.arange(0, t_max+t_step, t_step)
    reset_el1el2_err=False
    reset_el3el2_err= False
    
    target_el1el2=wd_row[string1]
    target_el3el2=wd_row[string2]
    logg=wd_row['logg']
    teff=wd_row['teff']
    teff_dist=np.random.normal(loc=teff, scale=wd_row['teff_err'], size=n_points)
    logg_dist= np.random.normal(loc=logg, scale=wd_row['logg_err'], size=n_points)
    if np.abs(wd_row[string1+'_err']) > limit_indicator:
        el1el2_dist=np.random.normal(loc=target_el1el2,scale=0,size=n_points) #the errorbar provided is too large one way or the other and so actually indicates a limit
        reset_el1el2_err=True
    else:
        el1el2_dist=np.random.normal(loc=target_el1el2,scale=wd_row[string1+'_err'],size=n_points)
    if np.abs(wd_row[string2+'_err']) > limit_indicator:
        el3el2_dist=np.random.normal(loc=target_el3el2,scale=0,size=n_points)
        reset_el3el2_err=True
    else:
        el3el2_dist=np.random.normal(loc=target_el3el2,scale=wd_row[string2+'_err'],size=n_points)
    def get_ssp(teff, logg, el1el2, el3el2, elements, plot_all=plot_all,tau_rand=False, tau_add=tau_add,modeler=modeler, atm_type=atm_type, overshoot=overshoot):
        tau_el1=itau.extrapolate_tau_x_logg(teff, logg, elements[0],modeler=modeler, atm_type=atm_type, overshoot=overshoot)
        tau_el2=itau.extrapolate_tau_x_logg(teff, logg, elements[1],modeler=modeler, atm_type=atm_type, overshoot=overshoot)
        tau_el3=itau.extrapolate_tau_x_logg(teff, logg, elements[2],modeler=modeler, atm_type=atm_type, overshoot=overshoot)
        
        tau2_tau1= tau_el2-tau_el1
        tau2_tau3=tau_el2-tau_el3
        if tau_rand:
            #print('tau2_tau1.shape',tau2_tau1.shape)
            tau2_tau1=np.random.normal(loc=tau2_tau1,scale=tau_add)
            tau2_tau3=np.random.normal(loc=tau2_tau3,scale=tau_add)
            print("tau_"+elements[0]+"-tau_"+elements[1], -1*np.mean(tau2_tau1),"+/-",np.std(tau2_tau1))
            print("tau_"+elements[2]+"-tau_"+elements[1], -1*np.mean(tau2_tau3),"+/-",np.std(tau2_tau3))
            #print('tau2_tau1.shape',tau2_tau1.shape)
        else:
            
            pass
        
        
        el1el2=el1el2+tau2_tau1
        el3el2=el3el2+tau2_tau3
        if plot_all:
            plt.hist(tau_el1-tau_el2, alpha=0.5,label='no addition',density=True)
            plt.hist(-1*tau2_tau1, alpha=0.5, label='addition',density=True)
            plt.title(r'$\tau$' + elements[0]+'-'+ elements[1])
            plt.legend()
            plt.show()
            
            plt.hist(tau_el3-tau_el2, label='no addition',density=True)
            plt.hist(-1*tau2_tau3, alpha=0.5, label='addition',density=True)
            plt.title(r'$\tau$' + elements[2]+'-'+ elements[1])
            plt.legend()
            plt.show()
            
        return el1el2, el3el2
    target_ssp_el1el2,target_ssp_el3el2=get_ssp(teff,logg, target_el1el2, target_el3el2,elements,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    el1el2_ssp_dist, el3el2_ssp_dist=get_ssp(teff_dist,logg_dist, el1el2_dist, el3el2_dist,elements,tau_rand=tau_rand,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    print('target_ssp', elements[0],elements[1],target_ssp_el1el2)
    print('target_ssp', elements[2], elements[1], target_ssp_el3el2)
    print('dist ssp',elements[0], elements[1], np.median(el1el2_ssp_dist), np.mean(el1el2_ssp_dist), np.std(el1el2_ssp_dist))
    print('dist ssp',elements[2], elements[1], np.median(el3el2_ssp_dist), np.mean(el3el2_ssp_dist), np.std(el3el2_ssp_dist))
    if plot_all:
        
        plt.hist(el1el2_ssp_dist)
        plt.axvline(target_ssp_el1el2)
        plt.title(elements[0]+'/'+elements[1])
        plt.show()
        get_ssp(teff_dist,logg_dist, el1el2_dist, el3el2_dist,elements,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
        plt.hist(el3el2_ssp_dist)
        plt.axvline(target_ssp_el3el2)
        plt.title(elements[2]+'/'+elements[1])
        plt.show()
        plt.scatter(el3el2_ssp_dist,el1el2_ssp_dist)
        plt.plot(target_ssp_el3el2,target_ssp_el1el2,marker='*', markersize=14)
        plt.show()
    el1el2_err= np.std(el1el2_ssp_dist)
    el3el2_err=np.std(el3el2_ssp_dist)
    if reset_el3el2_err:
        el3el2_err=wd_row[string2+'_err']
    if reset_el1el2_err:
        el1el2_err=wd_row[string1+'_err']
    else:
        pass
    return target_ssp_el1el2, target_ssp_el3el2,el1el2_err, el3el2_err


def get_ssp_accretion_rate(teff, logg, log_elabund, atm_type=default_atm_type, overshoot=default_overshoot, modeler=default_modeler,el_ratio_str='Ca/H',log_q=-5.0, m_wd=0.56, cross_extrap=default_cross_extrap):
    """
    inputs: teff, logg, atm_type, abundance, modeler, el (string indicating what the metal is), overshoot
    
    outputs: accretion rate in g/s 
    
    
    """
    el_pol,el_main=el_ratio_str.split('/') #obtain the different relevant element strings
    m_wd=(m_wd*const.M_sun).to(u.g).value #converting the white dwarf mass to grams from solar mass
    if cross_extrap:
        log_el_tau=itau.extrapolate_tau_x_logg(teff, logg, el_pol,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    else:
        log_el_tau=itau.extrapolate_single_el_tau(teff, el_pol, input_logg=logg,modeler=modeler, atm_type=atm_type, overshoot=overshoot)
    el_tau_seconds=((10.**log_el_tau)*u.yr).to(u.s).value
    log_accretion_rate=log_elabund+np.log10(pt.elements[cp.el_nums[el_pol]].mass/pt.elements[cp.el_nums[el_main]].mass)+log_q+np.log10(m_wd)-np.log10(el_tau_seconds)
    print('log accretion rate (g/s):', log_accretion_rate)
    accretion_rate=10.**log_accretion_rate
    
    return accretion_rate

if __name__ == '__main__':
    #wd_name='WDJ2356-209'
    target_logg=7.98
    target_logg_err=0.07
    target_teff= 4040. #K
    target_teff_err=110.
    
    #WDJ0212-5522
    target_teff=4590.
    target_teff_err=70.
    target_logg=7.97
    target_logg_err=0.02
    target_log_q=-5.5 #guessed value by Ben based on Koester2020 tables
    target_log_cah=-8.1 #plus or minus 0.3
    target_m_wd=0.56 #guessed by Ben based on logg, but should put through Bedard tables to check later
    
    macc= get_ssp_accretion_rate(target_teff, target_logg, target_log_cah, atm_type='H', el_ratio_str='Ca/H', log_q=target_log_q, m_wd=target_m_wd)
    
    
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
    tFeCa, dp_FeCa= el1el2_DP_el3el2(teff_dist,feca_dist, naca_dist, -0.01, 'Li', 'Ca', 'Na', logg=logg_dist, cross_extrap=True,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot)
    #plt.scatter(teff_dist, dp_FeCa)
    #plt.show()
    print("mean", np.mean(dp_FeCa), np.std(dp_FeCa),np.min(dp_FeCa), np.max(dp_FeCa))
    plt.hist(dp_FeCa, bins=101, density=True, label='log(Na/Ca)=-0.01', alpha=0.2)
    tFeCa, dp_FeCa= el1el2_DP_el3el2(teff_dist,feca_dist, naca_dist, -1.1, 'Li', 'Ca', 'Na', logg=logg_dist, cross_extrap=True,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot)
    plt.hist(dp_FeCa, bins=101, density=True, label='log(Na/Ca)=-1.1', alpha=0.2)
    print("mean", np.mean(dp_FeCa), np.std(dp_FeCa),np.min(dp_FeCa),np.max(dp_FeCa))
    plt.legend()
    plt.show()
    
    print(cp.el_nums)
    time_range=np.linspace(0, 40, 20)
    time_range=time_range*1e6

    dp_lica= declining_phase(3830., 1.7,time_range, 'Li', 'Ca',modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot)
    print('DP log(Li/Ca)', dp_lica)

    #plt.plot(np.log10(time_range), dp_lica)
    #plt.plot(time_range, dp_lica)
    #plt.plot(time_range, declining_phase(3830., 1.7,time_range, 'Li', 'Ca', steady_state_start=True), label='Steady state start')
    #plt.legend()
    #plt.show()

    t_NaCa= get_time_since_accretion(3830., 0.0, -1.1, 'Na', 'Ca', steady_state_start=False, cross_extrap=True, logg=7.77,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot)
    print("t_NaCa", t_NaCa, np.log10(t_NaCa))
    print(declining_phase(3830., 1.7,t_NaCa, 'Li', 'Ca', steady_state_start=False,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot))

    #t_NaCa= get_time_since_accretion(3830., 0.0, -1.0, 'Na', 'Ca', logg=7.5)
    print("t_NaCa", t_NaCa, np.log10(t_NaCa))
    print(declining_phase(3830., 1.7,t_NaCa, 'Li', 'Ca', steady_state_start=False, logg=7.5,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot))
    log_mCa= get_accreted_mass('Ca', -9.5, t_NaCa,  teff=3830., logg=8.0, log_q=-4.88, m_wd=0.56,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot)
    print(log_mCa)
    mCa=10.**log_mCa
    print(mCa)
    m_total=mCa/0.073
    print('total accreted mass', m_total, 'kg')
    log_mNa= get_accreted_mass('Na', -9.5, t_NaCa,  teff=3830., logg=8.0, log_q=-4.88, m_wd=0.56,modeler=default_modeler, atm_type=default_atm_type, overshoot=default_overshoot)
    print(log_mNa)
    mNa=10.**log_mNa
    m_total_Na=mNa/0.0033
    print('total accreted mass', m_total_Na, 'kg')
    print("mNa/mCa", mNa/mCa)
    print(np.log10(pt.elements[cp.el_nums['Na']].mass/pt.elements[cp.el_nums['Ca']].mass*mNa/mCa))
