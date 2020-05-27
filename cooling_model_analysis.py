"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-10-04

Plot the WD cooling models with Teff and logg values hopefully in an attempt to figure out which one best fits the 
target.



"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from astropy.table import Table, Column
import scipy.interpolate as scinterp


import spec_plot_tools as spt
import cal_params as cp
import hurley_polynomials as hp


cooling_model_file='COModel_ThinH.csv'

wd_file='temp_wd_abundances.csv'


cooling_model_file=cp.ref_dir+'WD_cooling_models/'+cooling_model_file

#default_ms_method='MIST'
#default_ms_method='Fontaine'
#default_ms_method='Hurley'
default_ms_method='Hurley_He'

output_dir= '/Users/BenKaiser/Desktop/'

export_total_ages=False

universe_age= 13.8 #Gyr
percent_range=0.68 #error bar coverage for total age estimate.
#percent_range=0.95
#null_age_val=20. #usually 20
null_age_val=50. #I'm experimenting though for the moment
default_limit_universe=True
default_randomize=True
massloss_default=True
#default_z=0.0001 #allegedly thick disk value
default_z=0.02 #approximately solar


#madeup
#wd_name='made-up'
#target_logg= 8.0
#target_logg_err=0.15
#target_teff= 6000. #K
#target_teff_err=190


##1330
wd_name='SDSSJ1330+6435'
target_logg= 8.26
target_logg_err=0.15
target_teff= 4310. #K
target_teff_err=190

#2356
#wd_name='WDJ2356-209'
#target_logg=7.98
#target_logg_err=0.07
#target_teff= 4040. #K
#target_teff_err=110.

#target_teff=3000.


##1644
#wd_name='GaiaJ1644-0449'
#target_logg=7.77
#target_logg_err= 0.23
#target_teff= 3830.
#target_teff_err= 230.

#J1150
#from Gentile Fusillo et al. 2019
#wd_name='SDSSJ1150+2403'
#target_logg=7.945462
#target_logg_err=2.24414e-01
#target_teff= 3454.718922
#target_teff_err= 126.986438

desired_NaCa= -1.1 #Sioux county meteorite, achondrite

############################3
##############################
prob_char=['+', '-', ' ','.'] #characters to be replaced
rep_char=['p','m','_',''] #characters to use to replace those other characters
def get_output_name(wd_name=wd_name):
    name_string=wd_name
    for prob, rep in zip(prob_char, rep_char):
        name_string=name_string.replace(prob,rep)
    return output_dir+name_string+'_tot_age_MC.csv'




#n=100000
n=int(1e6)


simon_mass= 0.45
simon_mass_err= 0.12

target_logg_dist= np.copy(np.random.normal(loc=target_logg, scale=target_logg_err, size=n))
target_teff_dist= np.copy(np.random.normal(loc=target_teff, scale=target_teff_err, size=n))
simon_mass_dist= np.copy(np.random.normal(loc=simon_mass, scale=simon_mass_err, size=n))
#simon_mass_dist[np.where(simon_mass_dist<0.2)] = np.nan

#target_logg=7.77
#target_teff= 4000.
#given_target_mass= 0.6
given_target_mass= 0.5


interp_kind='cubic'
#interp_kind='linear'

#target_logg=8.26
#target_teff= 4310. #K
##############################
#cummings_m_ranges= [
    #[[-np.inf,0.555],[np.nan,np.nan]],
    #[[0.555,0.717],[0.080,0.489]],
    #[[0.717,0.856],[0.187,0.184]],
    #[[0.856,1.24],[0.107,0.471]],
    #[[1.24,np.inf],[np.nan,np.nan]]
    #]

#setting the progenitor mass to be huge for masses greater than largest allowed so that the MS lifetime is essentially 0.
#setting the masses to produce Nan's for  progenitor mass if Mwd is below the range covered.
#I artificially changed the bounds in the below set so that the min mass is 0.52 This is decidedly outside Jeff's target area
#cummings_m_ranges= [
    #[[-np.inf,0.52],[np.nan,np.nan]],
    #[[0.52,0.717],[0.080,0.489]],
    #[[0.717,0.856],[0.187,0.184]],
    #[[0.856,1.24],[0.107,0.471]],
    #[[1.24,np.inf],[0.,0.]]
    #]
    
##Now that I've fixed the typos, I've determined that the actual minimum WD mass for which the result isn't stupid unreasonable is M_prog=M_wd, which is 0.532 (well slightly less than that, but 3 digits seems like enough)
#cummings_m_ranges= [
    #[[-np.inf,0.532],[np.nan,np.nan]],
    #[[0.532,0.717],[0.080,0.489]],
    #[[0.717,0.856],[0.187,0.184]],
    #[[0.856,1.24],[0.107,0.471]],
    #[[1.24,np.inf],[np.nan,np.nan]]
    #]
    
    
#those same values but with the uncertainties on each one also included. I'm not going to rescale the 
#boundaries to make it continuous for the randomized values as well because that would then require 
#randomizing each one and there's no guarantee that randomly selected values would even by continuous. 

cummings_m_ranges= [
    [[-np.inf,0.532],[np.nan,np.nan],[np.nan,np.nan]],
    [[0.532,0.717],[0.080,0.489],[0.016,0.030]],
    [[0.717,0.856],[0.187,0.184],[0.061,0.199]],
    [[0.856,1.24],[0.107,0.471],[0.016,0.077]],
    [[1.24,np.inf],[np.nan,np.nan],[np.nan,np.nan]]
    ]

def get_progenitor_mass(mass_wd, randomize=False, only_lose_mass=massloss_default):
    def mfunc(mass_wd, coeffs):
        return (mass_wd-coeffs[1])/coeffs[0]
    arrayvalue=True #if input is a float or not
    try:
        output_masses= np.ones(mass_wd.shape)
    except AttributeError:
        arrayvalue=False
        pass        
    for element in cummings_m_ranges:
        massrange= element[0]            
        if arrayvalue:
            inplay= np.where((mass_wd > massrange[0]) & (mass_wd <= massrange[1]))
            if randomize:
                random_element=[np.random.normal(loc=element[1][0], scale=element[2][0],size=mass_wd[inplay].shape),np.random.normal(loc=element[1][1], scale=element[2][1],size=mass_wd[inplay].shape)]
                #print("random_element", random_element)
                output_masses[inplay]= mfunc(mass_wd[inplay], random_element)
            else:
                output_masses[inplay]= mfunc(mass_wd[inplay], element[1])
            
        else:
            if((mass_wd > massrange[0])&(mass_wd < massrange[1])):
                output_masses= mfunc(mass_wd, element[1])
    if arrayvalue:
        if only_lose_mass:
            gained_mass=np.where((output_masses-mass_wd)<0)
            print("shape of those with calculated progenitors", output_masses.shape)
            output_masses[gained_mass]=np.nan
            print("shape of those that actually lost mass", output_masses[~np.isnan(output_masses)].shape)
        else:
            pass
    else:
        pass
    return output_masses

def get_ms_lifetime(mass_wd, method=default_ms_method, z=default_z, randomize=False):
    bins=np.linspace(0,10,1000)
    plt.hist(mass_wd, label='Original WD masses', bins=bins)
    prog_mass= get_progenitor_mass(mass_wd, randomize=randomize)
    plt.hist(prog_mass[~np.isnan(prog_mass)], label='progenitor masses remaining', bins=bins, alpha=0.5)
    plt.hist(mass_wd[~np.isnan(prog_mass)], label='WD masses remaining', bins=bins, alpha=0.5)
    plt.legend()
    plt.show()
    if method=='Fontaine':
        return 10*prog_mass**(-2.5)
    elif method=='MIST':
        return 61*prog_mass**(-2.5)
    elif method=='Hurley':
        return hp.get_t_ms(prog_mass, z=z)
    elif method=='Hurley_He':
        return hp.get_t_he(prog_mass, z)
    



#def get_ms_lifetime(mass_wd):
    #return 10*(8 *np.log(mass_wd/0.4))**(-2.5)
def get_old_ms_lifetime(mass_wd):
    return 10*(8 *np.log(mass_wd/0.4))**(-2.5)

def operate_on_dist(dist1, dist2, function):
    """
    dist1 is the distribution of the first input to 'function'
    
    dist2 is the distribution of the second input to 'function'
    
    function is the most likely scinterp.interp2d() function output or object or whatever that you've defined that you'd 
    like to feed distributions through and be able to track which inputs yielded them.
    
    For whatever reason, scinterp.interp2d() -generated objects seem to completely randomize the indices of the 
    outputs relative to the inputs. I have no idea why. It also wants to take the 2 input (10,) shape arrays and make 
    a (10,10) output, whose indices have no correlation to the indices of either input (10,) array. I know, that 
    seems ridiculous, but it's what I've been experiencing, so I just made a damn for-loop.
    
    """
    output_dist= []
    for el1, el2 in zip(dist1, dist2):
        out_el= function(el1, el2)
        output_dist.append(out_el)
    output_dist=np.array(output_dist).T[0]
    return output_dist

def clean_and_trim_age(total_age_dist, limit_universe=default_limit_universe, wd_mass_dist=[]):
    
    clean_total_age_dist= np.copy(total_age_dist)
    clean_total_age_dist[np.isnan(clean_total_age_dist)]=null_age_val
    clean_total_age_dist[np.where(clean_total_age_dist>null_age_val)] = null_age_val #set
    
    trimmed_total_age_dist=np.copy(total_age_dist)
    trimmed_wd_mass=np.copy(wd_mass_dist)
    try:
        trimmed_wd_mass= trimmed_wd_mass[~np.isnan(trimmed_total_age_dist)]
    except IndexError:
        #no wd_mass_dist_provided
        pass
    trimmed_total_age_dist=trimmed_total_age_dist[~np.isnan(trimmed_total_age_dist)]
    if limit_universe:
        try:
            trimmed_wd_mass= trimmed_wd_mass[trimmed_total_age_dist<universe_age]
            trimmed_total_age_dist=trimmed_total_age_dist[trimmed_total_age_dist<universe_age]
            return trimmed_total_age_dist, trimmed_wd_mass
        except IndexError:
            #no wd_mass_dist_provided
            pass
        #trimmed_total_age_dist=trimmed_total_age_dist[trimmed_total_age_dist<universe_age]
    else:
        try:
            print(wd_mass_dist)
            return trimmed_total_age_dist, trimmed_wd_mass
        except IndexError:
            pass
    return trimmed_total_age_dist
##########################3

cooling_table= Table.read(cooling_model_file)

#cooling_table.pprint()

#loggteff_to_m = scinterp.interp2d(cooling_table['Teff'], cooling_table['logg'], cooling_table['Mass'], kind=interp_kind)
loggteff_to_m_interp = scinterp.SmoothBivariateSpline(cooling_table['Teff'], cooling_table['logg'], cooling_table['Mass'])

teffm_to_age_interp= scinterp.interp2d(cooling_table['Teff'], cooling_table['Mass'], cooling_table['Age'], kind=interp_kind)
#teffm_to_age_interp= scinterp.SmoothBivariateSpline(cooling_table['Teff'], cooling_table['Mass'], cooling_table['Age'])

#loggteff_to_logTc_interp = scinterp.interp2d([cooling_table['Teff']], [cooling_table['logg']], cooling_table['Log(Tc)'], kind=interp_kind)
loggteff_to_logTc_interp = scinterp.SmoothBivariateSpline([cooling_table['Teff']], [cooling_table['logg']], cooling_table['Log(Tc)'])

def loggteff_to_m(teff, logg):
    return loggteff_to_m_interp(teff,logg)[0]

def teffm_to_age(teff, m):
    #return teffm_to_age_interp(teff, m)[0]
    output_vals=teffm_to_age_interp(teff, m)
    try:
        print(output_vals.shape[1])
        return output_vals.T[0]
    except IndexError:
        return output_vals

def loggteff_to_logTc(teff, logg):
    return loggteff_to_logTc_interp(teff, logg)[0]

#loggteff_to_logTc = scinterp.interp2d([cooling_table['Teff']], [cooling_table['logg']], cooling_table['Log(Tc)'], kind=interp_kind)

target_mass=  loggteff_to_m(target_teff, target_logg)
print('Target mass:', target_mass)

#target_mass_dist=loggteff_to_m(target_teff_dist, target_logg_dist)
target_mass_dist=operate_on_dist(target_teff_dist, target_logg_dist,loggteff_to_m)
print(target_mass_dist.shape)
print(target_logg_dist.shape)
print(target_teff_dist.shape)
print('Target mean mass:', np.mean(target_mass_dist), '+/-', np.std(target_mass_dist))


#calc_pro_masses= get_progenitor_mass(target_mass_dist)
test_masses= np.linspace(0.2, 1.3, 100)

#plt.scatter(target_mass_dist, calc_pro_masses, label='WD distribution')
#plt.plot(test_masses, get_progenitor_mass(test_masses), label='Cummings IFMR')
#plt.legend(loc='best')
#plt.xlabel('M_wd (M_sol)')
#plt.ylabel('M_prog. (M_sol)')
#plt.show()

#print(np.nanmax(calc_pro_masses))
#calc_pro_masses[np.isnan(calc_pro_masses)]=10.
#calc_pro_masses[np.isinf(calc_pro_masses)]=10.

#plt.hist(calc_pro_masses)
#plt.xlabel('M_progenitor')
#plt.show()

print(target_teff_dist.shape)
#age_bins=np.arange(0,21,0.25)
age_bins=np.arange(0,null_age_val+0.25, 0.25)
#print(np.where(target_mass_dist== loggteff_to_m(target_teff_dist[5], target_logg_dist[5])))
#print('comparison', target_mass_dist[5], loggteff_to_m(target_teff_dist[5], target_logg_dist[5]))
print('target_mass_dist.shape', target_mass_dist.shape)
mean_mass= np.nanmean(target_mass_dist)
std_mass= np.std(target_mass_dist)
median_mass= np.nanmedian(target_mass_dist)
print('Mass:', mean_mass, '+/-', std_mass, 'or', median_mass)


#plt.hist(target_mass_dist, bins=50, label='logg-> M')
#plt.hist(simon_mass_dist, bins=50,label='Simon M', alpha=0.5)
#plt.axvline(x=np.median(target_mass_dist), color='r', linestyle='--', label='median(logg -> M): '+str(np.round(np.median(target_mass_dist),2)))
#plt.axvline(x=np.median(simon_mass_dist), color='k', linestyle= '--', label= 'median(Simon M): '+ str(np.round(np.median(simon_mass_dist),2) ))
#plt.legend(loc='best')
#plt.xlabel('Mass')
#plt.ylabel('N')
#plt.show()





target_age= teffm_to_age(target_teff, target_mass)
print("target cooling age", target_age*1e-9, 'Gyr')

target_age_dist=operate_on_dist(target_teff_dist, target_mass_dist, teffm_to_age)*1e-9 #Gyr units
#target_age_dist=operate_on_dist(target_teff_dist, simon_mass_dist, teffm_to_age)*1e-9 #Gyr units
simon_age_dist=operate_on_dist(target_teff_dist, simon_mass_dist, teffm_to_age)*1e-9 #Gyr units

#ms_age_dist=get_ms_lifetime(target_mass_dist)
ms_age_dist=get_ms_lifetime(target_mass_dist, randomize=default_randomize)
#lowz_ms_age_dist=get_ms_lifetime(target_mass_dist, z=0.0001)
#ms_age_dist=get_ms_lifetime(simon_mass_dist)

total_age_dist= target_age_dist+ms_age_dist
#lowz_total_age_dist= target_age_dist+lowz_ms_age_dist

clean_total_age_dist= np.copy(total_age_dist)
#clean_total_age_dist[np.isnan(clean_total_age_dist)]=20.
#clean_total_age_dist[np.where(clean_total_age_dist> 20.)] = 20. #setting a max
clean_total_age_dist[np.isnan(clean_total_age_dist)]=null_age_val
clean_total_age_dist[np.where(clean_total_age_dist> null_age_val)] = null_age_val#setting a max


#cleanz_total_age_dist= np.copy(lowz_total_age_dist)
#cleanz_total_age_dist[np.isnan(cleanz_total_age_dist)]=20.
#cleanz_total_age_dist[np.where(cleanz_total_age_dist> 20.)] = 20. #setting a max


#################
################


#trimmed_total_age_dist=np.copy(total_age_dist)
#trimmed_total_age_dist=trimmed_total_age_dist[~np.isnan(trimmed_total_age_dist)]
#trimmed_total_age_dist=trimmed_total_age_dist[trimmed_total_age_dist<universe_age]

#trimmed_total_age_dist=clean_and_trim_age(total_age_dist)
trimmed_total_age_dist, trimmed_mass_dist=clean_and_trim_age(total_age_dist,wd_mass_dist= target_mass_dist)
throway_trimming, trimmed_teff_dist=clean_and_trim_age(total_age_dist,wd_mass_dist= target_teff_dist)
#trimmed_lowz_total_age_dist=clean_and_trim_age(lowz_total_age_dist)
throway_trimming, trimmed_cooling_dist=clean_and_trim_age(total_age_dist,wd_mass_dist= target_age_dist)
#trimmed_lowz_total_age_dist=clean_and_trim_age(lowz_total_age_dist)

#ln_trim_ages= np.log(trimmed_total_age_dist)

#plt.hist(ln_trim_ages)
#plt.xlabel('ln(total ages)')
#plt.show()

print('trimmed_total_age_dist.shape', trimmed_total_age_dist.shape)
print(wd_name)
print('relative remaining fraction', np.float_(trimmed_total_age_dist.shape[0])/total_age_dist.shape[0])

trim_vals, trim_edges, trim_patches= plt.hist(trimmed_total_age_dist, bins=np.arange(0,null_age_val+0.1, 0.1), label='total ages limited to universe', normed=True, alpha=0.2)
sub_edges=trim_edges[:-1] #remove the last edge to make the length the same as the probability values
sort_order= np.argsort(-1*trim_vals) #by multiplying by a negative you make the largest values the smallest effectively. Thanks stackoverflow!
bin_widths=np.median(trim_edges-np.roll(trim_edges,1))
sort_trim_vals= trim_vals[sort_order]
sort_trim_probs=sort_trim_vals*bin_widths #probability contribution of each bin
sort_edges=sub_edges[sort_order]
max_arg=np.argmax(trim_vals)
max_val=sub_edges[max_arg]+(bin_widths*0.5)
cumprob= np.cumsum(sort_trim_probs)
inbounds=np.where(cumprob-percent_range < 0)
print("max percentage", np.max(cumprob[inbounds]))
print("total percentage in distribution", np.max(cumprob))


print('bin_widths', bin_widths)
print('max val', max_val)


if export_total_ages:
    print('\n\n************\n')

    output_name= get_output_name()
    print('Saving', output_name)
    np.savetxt(output_name, trimmed_total_age_dist)
    print('Saved')

    print('\n************\n\n')
else:
    print('Skipping saving of MC total ages because export_total_ages=',export_total_ages)


lowbound=np.nanmin(sort_edges[inbounds])
highbound= np.nanmax(sort_edges[inbounds])+bin_widths
bounding_range=np.where((trimmed_total_age_dist < highbound) & (trimmed_total_age_dist > lowbound))
bounded_masses= trimmed_mass_dist[bounding_range]
bounded_cooling=trimmed_cooling_dist[bounding_range]
ml_mass_vals= trimmed_mass_dist[np.where((trimmed_total_age_dist > sub_edges[max_arg])&(trimmed_total_age_dist < sub_edges[max_arg]+bin_widths))]
ml_cooling_vals=trimmed_cooling_dist[np.where((trimmed_total_age_dist > sub_edges[max_arg])&(trimmed_total_age_dist < sub_edges[max_arg]+bin_widths))]

print('min WD mass included in age of universe:', np.min(trimmed_mass_dist))
minarg= np.argmin(trimmed_mass_dist)
print('Teff of that min WD mass:', trimmed_teff_dist[minarg])
print('total age of that min WD mass:', trimmed_total_age_dist[minarg])
print("\n***********\n")
print("low total age:", lowbound)
print("high total age:", highbound)
print("most likely total age:", max_val)
print("Minimum WD mass in that age range", np.min(bounded_masses))
print("Maximum WD mass in that age range", np.max(bounded_masses))
print("Mean WD mass for M/L", np.mean(ml_mass_vals))
print(np.min(ml_mass_vals), np.max(ml_mass_vals))

print("Mean M/L cooling vals", np.mean(ml_cooling_vals))
print(np.min(ml_cooling_vals), np.max(ml_cooling_vals))
print("\n***********\n")
plt.show()

these_bins= np.arange(0,1.4, 0.025)
plt.hist(bounded_masses, alpha =0.2, label='bounded_masses', bins=these_bins, normed=True)
plt.hist(target_mass_dist, alpha=0.2, label='target_mass_dist', bins=these_bins, normed=True)
plt.hist(trimmed_mass_dist, alpha=0.2, label='trimmed_mass_dist', bins=these_bins, normed=True)
plt.legend()
plt.title(wd_name)

plt.xlabel('Mass (M_sol)')
plt.show()
plt.plot(sort_trim_vals)
plt.show()

plt.hist(trimmed_teff_dist, normed=True)
plt.xlabel('Teff of trimmed masses')
plt.show()



print('Minimum total age: ', np.nanmin(trimmed_total_age_dist), 'Gyr')
print('Median total age: ', np.nanmedian(trimmed_total_age_dist), 'Gyr')


#prog_bins=np.arange(0.0,12, 0.05)
##prog_bins=50
#plt.hist(clean_and_trim_age(get_progenitor_mass(target_mass_dist),limit_universe=False)[0], alpha=0.2, bins=prog_bins, normed=True, label='prog. mass for all logg in defined IFMR')
#plt.hist(clean_and_trim_age(get_progenitor_mass(bounded_masses), limit_universe=False)[0], alpha=0.2, bins=prog_bins, normed=True, label='prog. mass for WD masses in M/L total age bin')
#plt.hist(clean_and_trim_age(get_progenitor_mass(trimmed_mass_dist), limit_universe=False)[0], alpha=0.2, bins=prog_bins, normed=True, label='prog. mass for total ages < age of universe')
#plt.legend()
#plt.title(wd_name)
#plt.xlabel('Prog. Mass (M_sol)')
#plt.show()


test_wd_masses=np.linspace(0.25, 1.5, 1000)
plt.plot(test_wd_masses, get_progenitor_mass(test_wd_masses),marker='o')
#plt.plot(test_wd_masses, clean_and_trim_age(get_progenitor_mass(test_wd_masses),limit_universe=False))
plt.ylabel('Progenitor Mass (M_sol)')
plt.xlabel('WD Mass (M_sol)')
plt.title('InverseCummings et al. 2018 IFMR with downside extended down to 0.532 from default 0.55')
plt.show()

plt.plot(get_progenitor_mass(test_wd_masses),test_wd_masses, marker='o')
#plt.plot(test_wd_masses, clean_and_trim_age(get_progenitor_mass(test_wd_masses),limit_universe=False))
plt.xlabel('Progenitor Mass (M_sol)')
plt.ylabel('WD Mass (M_sol)')
plt.title('Inverted Inverse Cummings et al. 2018 IFMR with downside extended down to 0.532 from default 0.55')
plt.show()

############################
###########################

#total_simon_dist= get_ms_lifetime(simon_mass_dist)+simon_age_dist
#total_simon_dist[np.where(total_simon_dist > null_age_val)]=null_age_val
#total_simon_dist[np.isnan(total_simon_dist )]=null_age_val

mean_age= np.nanmean(target_age_dist)
std_age= np.std(target_age_dist)
mean_total_age=np.nanmean(total_age_dist)
std_total_age=np.nanstd(total_age_dist)

med_total_age=np.nanmedian(clean_total_age_dist)
upper_total_age=np.nanpercentile(clean_total_age_dist, 84)
lower_total_age= np.nanpercentile(clean_total_age_dist, 16)

print('med total age', med_total_age, 'up to', upper_total_age, 'or down to', lower_total_age)
#print('\n99\%\ chance that total age > ', np.nanpercentile(clean_total_age_dist, 1),'Gyr\n')
#print('\n99\%\ chance that total age < ', np.nanpercentile(clean_total_age_dist, 99),'Gyr\n')
print('\n99\%\ chance that total age > ', np.nanpercentile(trimmed_total_age_dist, 1),'Gyr\n')
print('\n99\%\ chance that total age < ', np.nanpercentile(trimmed_total_age_dist, 99),'Gyr\n')
print('mean cooling age', mean_age, '+/-', std_age)
print('mean total age', mean_total_age, '+/-', std_total_age)
print(np.nanpercentile(clean_total_age_dist,16),np.nanmedian(clean_total_age_dist), np.nanpercentile(clean_total_age_dist,84))
#print(np.nanmedian(total_simon_dist), np.nanpercentile(total_simon_dist,16), np.nanpercentile(total_simon_dist,84))
plt.hist(target_age_dist, bins=np.arange(0,null_age_val+0.25, 0.25),label='cooling ages', normed=True)
#plt.hist(ms_age_dist[~np.isnan(ms_age_dist)], bins=50, alpha=0.4, label='ms ages')
plt.hist(clean_total_age_dist, bins=np.arange(0,null_age_val+0.25, 0.25), alpha=0.5,label= 'logg-> M -> Total Ages', normed=True)
plt.hist(trimmed_total_age_dist, bins=np.arange(0,null_age_val+0.1, 0.1), label='total ages limited to universe', normed=True, alpha=0.2)
#plt.hist(total_simon_dist, bins=np.arange(0,21, 0.25), alpha=0.5, label='Simon M -> total ages', normed=True)
#plt.hist(trimmed_lowz_total_age_dist, bins=np.arange(0,null_age_val+0.1, 0.1), label='lowz total ages limited to universe', normed=True, alpha=0.2)
plt.xlabel('Age (Gyr)')
plt.axvline(x=med_total_age, linestyle='--', color='k')
plt.axvline(x=upper_total_age, linestyle='--', color='k')
plt.axvline(x=lower_total_age, linestyle='--', color='k')
plt.axvline(x=max_val, linestyle='--', color='r', label='M/L estimates')
plt.axvline(x=lowbound, linestyle='--', color='r')
plt.axvline(x=highbound, linestyle='--', color='r')
plt.legend(loc='best')
plt.title(wd_name)
#plt.yscale('log')
plt.show()

log_total= np.log10(total_age_dist)
clean_log_dist= np.copy(log_total)
clean_log_dist[np.isnan(clean_log_dist)]=null_age_val

#clean_log_dist[np.where(clean_total_age_dist> 20.)] = 20.
#plt.hist(clean_log_dist, bins=200, normed=True, color='g', alpha=0.5, label='Total Ages')
plt.hist(clean_log_dist, bins=np.linspace(0,null_age_val+1,1000), normed=True, color='g', alpha=0.5, label='Total Ages')
plt.xlabel('log10(age(Gyr))')
#plt.yscale('log')
plt.show()

target_age_gmass= teffm_to_age(target_teff, given_target_mass)
print('Target age:', target_age)
print('Target age assuming mass=', given_target_mass, ':', target_age_gmass)
target_logTc=  loggteff_to_logTc(target_teff, target_logg)
target_Tc= 10.** target_logTc

print('Target core temperature:', target_Tc)


#loggteff_to_age= scinterp.interp2d(cooling_table['Teff'], cooling_table['logg'], cooling_table['Age'], kind=interp_kind)
loggteff_to_age_interp= scinterp.SmoothBivariateSpline(cooling_table['Teff'], cooling_table['logg'], cooling_table['Age'])

def loggteff_to_age(teff, logg):
    return loggteff_to_age_interp(teff, logg)[0]

target_age2= loggteff_to_age(target_teff, target_logg)
ms_lifetime= get_ms_lifetime(target_mass)

print("Target age from logg and teff:", target_age2)
print("MS lifetime:", ms_lifetime, "Gyr")
print("Total age from logg and teff:", ms_lifetime+(target_age2*1e-9))
#print("Total age from given mass:", get_ms_lifetime(given_target_mass)+(target_age_gmass*1e-9), 'Gyr')

approx_inds= np.where((cooling_table['Teff']< 4000) & (cooling_table['Teff']> 3500))
approx_masses= cooling_table['Mass'][approx_inds]
approx_ages=cooling_table['Age'][approx_inds]*1e-9

#wd_mass_vals= np.linspace(0.2, 1.3, 100)
#wd_mass_vals= np.linspace(0.2, 1.3, 1000)
wd_mass_vals=10.**np.linspace(np.log10(0.2),np.log10(1.3), 1000)
cooling_ages=teffm_to_age(target_teff, wd_mass_vals)*1e-9
print('cooling_ages.shape', cooling_ages.shape)
ms_ages= get_ms_lifetime(wd_mass_vals)
#lowz_ms_ages=get_ms_lifetime(wd_mass_vals, z=0.0001)
print(ms_ages.shape)

total_ages= cooling_ages+ms_ages
#lowz_total_ages= cooling_ages+lowz_ms_ages
print('total_ages.shape', total_ages.shape)
print((teffm_to_age(target_teff, wd_mass_vals)*1e-9).shape)
#print(get_ms_lifetime(wd_mass_vals).shape)
#plt.plot(wd_mass_vals, get_ms_lifetime(wd_mass_vals))
#plt.axvline(x=0.5, linestyle='--', color='k')
#plt.axhline(y=10, linestyle='--', color='k')
#plt.axvline(x=spt.naninfmax(cummings_m_ranges), color='r', linestyle='--', label=r'IFMR $M_{WD}$ Range')
#plt.axvline(x=spt.naninfmin(cummings_m_ranges), color='r', linestyle='--')
plt.axhline(y=13.8, color='k', label='Age of the Universe', linestyle= '--' )
plt.plot(wd_mass_vals, total_ages, label='Total Age')
#plt.plot(wd_mass_vals, lowz_total_ages, label='Total Age Z='+str(0.0001))
plt.plot(wd_mass_vals, cooling_ages, label='WD Cooling Age')
#plt.plot(wd_mass_vals, ms_ages, label='Z='+str(default_z)+'MS lifetime from' +default_ms_method)
plt.plot(wd_mass_vals, ms_ages, label=r'Progenitor $t_{BGB} + t_{He}$')
#plt.plot(wd_mass_vals, lowz_ms_ages, label='Z='+str(0.0001)+' MS lifetime from' +default_ms_method)
#plt.plot(wd_mass_vals, get_ms_lifetime(wd_mass_vals, method='Fontaine'), label='Fontaine')
#plt.scatter(approx_masses, approx_ages, color='r', label='Grid vals with Teff ~3800K')
plt.hist(target_mass_dist, normed=True, label=wd_name+' MC Masses', color='k')
#plt.scatter(0.56,  teffm_to_age(target_teff, 0.56)*1e-9+get_ms_lifetime(0.56), label='M=0.56 at teff'+str(target_teff))
plt.xlabel(r'$M_{wd}$ $(M_{\odot})$')
#plt.ylabel('MS lifetime (Gyr)')
plt.ylabel('Age (Gyr)')
plt.legend()
#plt.yscale('log')
plt.ylim(0,15)
#plt.title(wd_name)
plt.show()

approx_inds= np.where((cooling_table['Teff']< 4000) & (cooling_table['Teff']> 3500))
approx_masses= cooling_table['Mass'][approx_inds]
approx_loggs= cooling_table['logg'][approx_inds]
approx_ages=cooling_table['Age'][approx_inds]*1e-9
approx_teffs= cooling_table['Teff'][approx_inds]
#other_inds= np.where((approx_loggs <  8.1 ) & (approx_loggs > 7.5))
#other_masses=approx_masses[other_inds]
#other_teffs= approx_teffs[other_inds]
#other_masses= approx_masses[other_inds]
#other_ages=approx_ages[other_inds]
plt.scatter(approx_loggs, approx_masses, label='cooling models')
test_loggs= np.linspace(7.0, 9.5, 100)
test_masses= loggteff_to_m(target_teff, test_loggs)
plt.plot(test_loggs, test_masses, label='for Teff = '+str(target_teff))
plt.legend(loc='best')
plt.xlabel('log g')
plt.ylabel('M_wd (M_sol)')
plt.show()

plt.scatter(approx_loggs, approx_ages, label='cooling models')
test_ages= loggteff_to_age(target_teff, test_loggs)*1e-9
other_test_ages= teffm_to_age(target_teff, test_masses)*1e-9
plt.plot(test_loggs, test_ages, label='logg-> Age for Teff = '+str(target_teff))
plt.plot(test_loggs, other_test_ages, label='logg -> M -> Age')
plt.legend(loc='best')
plt.xlabel('log g')
plt.ylabel('Age (Gyr)')
plt.show()


plt.scatter(cooling_table['Teff'], cooling_table['logg'], label='cooling models')
plt.plot(target_teff, target_logg, marker='*', color='r', markersize=12)
plt.xlabel('Teff (K)')
plt.ylabel('log(g)')
plt.title(cooling_model_file.split('/')[-1])
plt.show()


plt.scatter(cooling_table['Teff'], cooling_table['Mass'], label='cooling models')
plt.plot(target_teff, target_mass, marker='*', color='r', markersize=12)
plt.xlabel('Teff (K)')
plt.ylabel('Mass')
plt.title(cooling_model_file.split('/')[-1])

plt.show()


given_mass= given_target_mass
given_inds= np.where(cooling_table["Mass"]==given_mass)
plt.scatter(cooling_table['Teff'][given_inds], cooling_table['Age'][given_inds], label='table vals for '+ str(given_mass))
#print((cooling_table['Teff'][given_inds]).shape, teffm_to_age(cooling_table['Teff'][given_inds], given_mass))
plt.scatter(cooling_table['Teff'][given_inds], teffm_to_age(cooling_table['Teff'][given_inds], given_mass), label='interpolated vals', color='r')
plt.plot(target_teff, target_age, marker='*', color='r', markersize=12, label='Target')
plt.plot(target_teff, target_age2, marker='*', color='g', markersize=12, label='Target (from teff+logg)')
plt.legend()
plt.xlabel('Teff')
plt.xlim(0,70000)
plt.ylabel('Age')
plt.show()


mass_vals= np.linspace(0.2,1.2, 100)
plt.plot(mass_vals, teffm_to_age(target_teff, mass_vals), label='interpolated vals for teff='+str(target_teff))
plt.plot(target_mass, target_age, marker='*', color='r', markersize=12, label='Target')
plt.legend()
plt.xlabel('mass')
plt.ylabel('age')
plt.show()

teff_vals= np.linspace(3000., 30000., 100)
allowed_inds= np.where(cooling_table['Mass']==given_target_mass)
print('shapes again', teff_vals.shape, target_logg, loggteff_to_logTc(teff_vals, target_logg).shape)
target_logg_array= np.ones(teff_vals.shape[0])*target_logg
plt.plot(teff_vals, loggteff_to_logTc(teff_vals, target_logg_array), label='log(Tc)')
plt.scatter(cooling_table['Teff'][allowed_inds], cooling_table['Log(Tc)'][allowed_inds], label='cooling table direct values for '+str(given_target_mass)+ 'M')
plt.plot(target_teff, target_logTc, marker='*', label='target')
plt.xlabel('Teff')

plt.ylabel('log(Tc)')
plt.legend()
plt.show()


logg_vals= np.linspace(7., 8.5, 100)
allowed_inds= np.where(cooling_table['Mass']==given_target_mass)
plt.plot(logg_vals, loggteff_to_logTc(target_teff, logg_vals), label='log(Tc)')
plt.scatter(cooling_table['logg'][allowed_inds], cooling_table['Log(Tc)'][allowed_inds], label='cooling table direct values for '+str(given_target_mass)+ 'M')
plt.plot(target_logg, target_logTc, marker='*', label='target')
plt.xlabel('logg')

plt.ylabel('log(Tc)')
plt.legend()
plt.show()


print('t_ms for 0.56 M_wd', get_ms_lifetime(0.56,method='Hurley'))





