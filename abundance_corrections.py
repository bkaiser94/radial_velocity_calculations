"""
Created by Ben Kaiser (UNC Chapel Hill) 2020-01-19

This should contain all of the functions that perform the diffusion corrections on the abundances


"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column

import interp_tau as itau

def declining_phase(target_teff, log_el1_over_el2, time,el1, el2, logg=8.0, steady_state_start=False):
    """
    provide log_el1_over_el2 as an absolute number ratio not the one normalized to solar abundances.
    
    the output value will be in log10 of absolute number abundances... hopefully
    """
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


def get_time_since_accretion(target_teff, log_atm_ratio, log_desired_ratio, el1, el2,logg=8.0, steady_state_start=False):
    
    el1_logtau= itau.extrapolate_single_el_tau(target_teff, el1, input_logg=logg)
    el2_logtau= itau.extrapolate_single_el_tau(target_teff, el2, input_logg=logg)
    el1_tau=10.**el1_logtau
    el2_tau=10.**el2_logtau
    t_coeffs= (el1_tau*el2_tau)/(el2_tau-el1_tau)
    pollutant_term=log_desired_ratio- log_atm_ratio
    if steady_state_start:
        pollutant_term=pollutant_term+el1_logtau-el2_logtau
    else:
        pass
    pollutant_term=pollutant_term*np.log(10.)
    return t_coeffs*pollutant_term


def LiCa_DP_NaCa(target_teff, log_LiCa, log_NaCa, desired_log_NaCa, logg=8.0, steady_state_start=False):
    t_NaCa= get_time_since_accretion(target_teff, log_NaCa, desired_log_NaCa, "Na", "Ca",  logg=logg, steady_state_start=steady_state_start)
    dp_LiCa= declining_phase(target_teff, log_LiCa,t_NaCa,  'Li', "Ca", logg=logg, steady_state_start=steady_state_start)
    
    return t_NaCa, dp_LiCa


if __name__ == '__main__':
    time_range=np.linspace(0, 40, 20)
    time_range=time_range*1e6

    dp_lica= declining_phase(3830., 1.7,time_range, 'Li', 'Ca')
    print('DP log(Li/Ca)', dp_lica)

    #plt.plot(np.log10(time_range), dp_lica)
    plt.plot(time_range, dp_lica)
    plt.plot(time_range, declining_phase(3830., 1.7,time_range, 'Li', 'Ca', steady_state_start=True), label='Steady state start')
    plt.legend()
    plt.show()

    t_NaCa= get_time_since_accretion(3830., 0.0, -1.0, 'Na', 'Ca')
    print("t_NaCa", t_NaCa, np.log10(t_NaCa))
    print(declining_phase(3830., 1.7,t_NaCa, 'Li', 'Ca', steady_state_start=False))

    t_NaCa= get_time_since_accretion(3830., 0.0, -1.0, 'Na', 'Ca', logg=7.5)
    print("t_NaCa", t_NaCa, np.log10(t_NaCa))
    print(declining_phase(3830., 1.7,t_NaCa, 'Li', 'Ca', steady_state_start=False, logg=7.5))
