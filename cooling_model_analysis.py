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

import time


import spec_plot_tools as spt
import cal_params as cp
import hurley_polynomials as hp


#print(sys.version.split())
python_version_long=sys.version.split()[0]
python_digits=python_version_long.split('.')
python_version_number=float('.'.join([python_digits[0],python_digits[1]]))
#print(python_version_number)

#wd_file='temp_wd_abundances.csv'
wd_file='20210527_all_wd_abundances_new_Blouin.csv'

atm_type='thinH'
#atm_type='thickH'

#cooling_modeler='fontaine2001'
cooling_modeler='bedard2020'

#cooling_model_file=cp.ref_dir+'WD_cooling_models/'+cooling_model_file

#default_ms_method='MIST'
#default_ms_method='Fontaine'
#default_ms_method='Hurley'
default_ms_method='Hurley_He'

output_dir= '/Users/BenKaiser/Desktop/'
output_dir=output_dir+'wd_70_30_total_age_partial_contributions_each_atm_type/'

export_total_ages=False

#universe_age= 13.8 #Gyr
universe_age=20. #Gyr, artifical age of the universe to use as upper limit to get around Malmquist bias
percent_range=0.68 #error bar coverage for total age estimate.
#percent_range=0.95
#null_age_val=20. #usually 20
null_age_val=50. #I'm experimenting though for the moment
default_limit_universe=True
default_randomize=True
massloss_default=True
#default_z=0.0001 #allegedly thick disk value
default_z=0.02 #approximately solar
#default_z=2.0 #super enriched


##madeup
#wd_name='made-up'
#target_logg= 8.0
#target_logg_err=0.15
#target_teff= 3000. #K
#target_teff_err=180.

##
##1644
#wd_name='WD J1644-0449'
#target_logg=7.77
#target_logg_err= 0.23
#target_teff= 3830.
#target_teff_err= 230.


##1644 w/ eDR3 approximation applied
#wd_name='WD J1644-0449'
#target_logg=7.85
#target_logg_err= 0.23
#target_teff= 3830.
#target_teff_err= 230.


##1330
#wd_name='SDSSJ1330+6435'
#target_logg= 8.26
#target_logg_err=0.15
#target_teff= 4310. #K
#target_teff_err=190.

####2356
#wd_name='WDJ2356-209'
#target_logg=7.98
#target_logg_err=0.07
#target_teff= 4040. #K
#target_teff_err=110.

###J1636
wd_name='SDSSJ1636+1619'
target_logg= 8.10
target_logg_err=0.06
target_teff= 4410. #K
target_teff_err=200.

#J1636 from Table2 extended Blouin Magnesium 2020
wd_name='SDSSJ1636+1619'
target_logg= 8.096
target_logg_err=0.059
target_teff= 4410. #K
target_teff_err=200.



##J2317 Hollands et al 2021 parameters
#wd_name='WDJ2317+1830'
#target_logg=8.64
#target_logg_err=0.03
#target_teff= 4210. #K
#target_teff_err=50.

#####J2317 Simon's new parameters
#wd_name='WDJ2317+1830'
#target_logg=8.74
#target_logg_err=0.06
#target_teff= 4430. #K
#target_teff_err=120.

#J1824 Hollands et al 2021 parameters
#wd_name='WDJ1824+1213'
#target_logg=7.41
#target_logg_err=0.07
#target_teff= 3350. #K
#target_teff_err=50.

#####J1824 Simon's new parameters
#wd_name='WDJ1824+1213'
#target_logg=7.53
#target_logg_err=0.09
#target_teff= 3540. #K
#target_teff_err=90.

##LHS2534 Hollands et al 2021 parameters
#wd_name='LHS2534'
#target_logg=7.97
#target_logg_err=0.04
#target_teff= 4780. #K
#target_teff_err=50.

##LHS2534 Simon's new parameters
#wd_name='LHS2534'
#target_logg=8.101
#target_logg_err=0.075
#target_teff= 5020. #K
#target_teff_err=100.



### WD J1922+0233 Tremblay et al. 2020 parameters
#wd_name='WDJ1922+0233'
#target_logg= 9.1
#target_logg_err=0.02
#target_teff=5800.
#target_teff_err=390.



#J1150
#from Gentile Fusillo et al. 2019
#wd_name='SDSSJ1150+2403'
#target_logg=7.945462
#target_logg_err=2.24414e-01
#target_teff= 3454.718922
#target_teff_err= 126.986438

#J2339-0424 from Klein et al. 2021
#wd_name='GALEX J2339-0424'
#target_logg=7.93
#target_logg_err=0.09
#target_teff= 13735
#target_teff_err= 500

##GD378 from Klein et al. 2021
#wd_name='GD 378'
#target_logg=7.93
#target_logg_err=0.06
#target_teff= 15620
#target_teff_err= 500



#Made up white dwarf for total age estimates
#wd_name='Imaginary'
#target_logg=8.74
#target_logg_err=0.06
#target_teff= 7800. #K
#target_teff_err=300.

### WD 1856+534 Xu et al. 2021 (modeled by Blouin; white dwarf with eclipsing gas giant candidate)
#wd_name='WD 1856+534'
#target_logg=7.995
#target_logg_err=0.065
#target_teff=4860.
#target_teff_err=80.

##WD 1425+540 Bonsor et al. 2021 (polluted white dwarf with wide binary companion),
#wd_name='WD1425+540'
#target_logg=8.04
#target_logg_err=0.05 #arbitrarily guessed
#target_teff= 14213. #K
#target_teff_err=300. #arbitarily guessed

desired_NaCa= -1.1 #Sioux county meteorite, achondrite


########################

def get_wd_parameters(wd_name, wd_file=wd_file):
    
    
    
    return

#####################

#cooling_model_file='COModel_ThinH.csv'
#cooling_model_file='COModel_ThickH.csv'
#cooling_model_file='bedard2020_seq_thickH.csv'
#cooling_model_file='bedard2020_seq_thinH.csv'



#interp_kind='cubic'
interp_kind='quintic'

if cooling_modeler== 'bedard2020':
    if atm_type=='thinH':
        cooling_model_file='bedard2020_seq_thinH.csv'
        interp_kind='cubic'
        print("\n\n*****\n\nbedard2020 with thinH atmosphere selected.\nSetting interp_kind='cubic'\nBecause I know it works correctly-ish.")
    elif atm_type=='thickH':
        cooling_model_file='bedard2020_seq_thickH.csv'
elif cooling_modeler=='fontaine2001':
    if atm_type=='thinH':
        cooling_model_file='COModel_ThinH.csv'
        print("\n\n*****\n\nfontaine2001 with thinH atmosphere selected.\nSetting interp_kind='cubic'\nBecause I know it works correctly-ish.")
        interp_kind='cubic'

    elif atm_type=='thickH':
        cooling_model_file='COModel_ThickH.csv'
        print("\n\n fontaine2001 with thickH atmosphere doesn't actually work correctly for some reason with any interpolation, so... we're gonna stop here. Pick a different setting. We also have the bedard2020 cooling_modeler option, which I have hopefully gotten running by the time you see this message...\n\n")
        sys.exit()
else:
    print('\n\n**********\n\nInvalid cooling_modeler and atm_type selected:')
    print('cooling_modeler:',cooling_modeler,'atm_type:',atm_type,'\n********\n\n')

cooling_model_file=cp.ref_dir+'WD_cooling_models/'+cooling_model_file

############################3
##############################
prob_char=['+', '-', ' ','.'] #characters to be replaced
rep_char=['p','m','_',''] #characters to use to replace those other characters
def get_output_name(wd_name=wd_name):
    name_string=wd_name
    for prob, rep in zip(prob_char, rep_char):
        name_string=name_string.replace(prob,rep)
    time_string= str(time.time()).split('.')[0]
    #return output_dir+name_string+'_tot_age_MC.csv'
    #return output_dir+time_string+'_'+name_string+'_tot_age_MC.csv'
    return output_dir+time_string+'_'+name_string+'_'+cooling_modeler+'_'+atm_type+'_tot_age_MC.csv'





#n=100000
#n=int(1e6)
#n=int(4e6)
#n=1000
#n=10
#n=int(6e6)#6 million for the 60% M_H>10^-10 from Cunningham et al. (2020)
n=int(2.5e6)#2.5 million for the 25% 10^-14< M_H<10^-10 from Cunningham et al. (2020)


#bin_widths=0.01
#age_bin_widths=0.01
#bin_widths=1.0 
#age_bin_widths=1.0
bin_widths=0.1 #default width of age bin widths for histogram.
age_bin_widths=0.1
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

test_mass=0.5


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
#randomizing each one and there's no guarantee that randomly selected values would even be continuous. 
######2021-05-28 This was the version (below this comment) that was in place for Kaiser et al. 2020/2021 (the Science paper on the discovery of lithium)
cummings_m_ranges= [
    [[-np.inf,0.532],[np.nan,np.nan],[np.nan,np.nan]],
    [[0.532,0.717],[0.080,0.489],[0.016,0.030]],
    [[0.717,0.856],[0.187,0.184],[0.061,0.199]],
    [[0.856,1.24],[0.107,0.471],[0.016,0.077]],
    [[1.24,np.inf],[np.nan,np.nan],[np.nan,np.nan]]
    ]


####2021-05-28 I'm going to get experimental here. Technically other checks outside the IFMR lower boundary should remove unphysical solutions, so I'm going to lower the IFMR lower boundary to a mass that would normally be mass-gaining for MS to WD, but remember the IFMR has randomized boundaries, so I'm going to try letting just the testing mass loss requirement perform the lower truncation in coordination with confinement to whatever total age.


#cummings_m_ranges= [
    #[[-np.inf,0.0001],[np.nan,np.nan],[np.nan,np.nan]],
    #[[0.0001,0.717],[0.080,0.489],[0.016,0.030]],
    #[[0.717,0.856],[0.187,0.184],[0.061,0.199]],
    #[[0.856,1.24],[0.107,0.471],[0.016,0.077]],
    #[[1.24,np.inf],[np.nan,np.nan],[np.nan,np.nan]]
    #]

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

mteff_to_age_interp=scinterp.interp2d(cooling_table['Mass'], cooling_table['Teff'], cooling_table['Age'], kind=interp_kind)

#loggteff_to_logTc_interp = scinterp.interp2d([cooling_table['Teff']], [cooling_table['logg']], cooling_table['Log(Tc)'], kind=interp_kind)
loggteff_to_logTc_interp = scinterp.SmoothBivariateSpline([cooling_table['Teff']], [cooling_table['logg']], cooling_table['Log(Tc)'])

def loggteff_to_m(teff, logg):
    return loggteff_to_m_interp(teff,logg)[0]

def teffm_to_age(teff, m):
    #return teffm_to_age_interp(teff, m)[0]
    output_vals=teffm_to_age_interp(teff, m)
    #print('teffm_to_age',output_vals)
    #print('teff_to_age.shape',output_vals.shape)
    try:
        print(output_vals.shape[1])
        return output_vals.T[0]
        #return output_vals[0]
    except IndexError:
        return output_vals
    
def mteff_to_age(m,teff):
    #return teffm_to_age_interp(teff, m)[0]
    output_vals=mteff_to_age_interp(m,teff)
    #print('mteff_to_age',output_vals)
    #print('mteff_to_age.shape',output_vals.shape)
    try:
        print(output_vals.shape[1])
        #return output_vals.T[0]
        return output_vals[0]
    except IndexError:
        return output_vals

def loggteff_to_logTc(teff, logg):
    return loggteff_to_logTc_interp(teff, logg)[0]

#loggteff_to_logTc = scinterp.interp2d([cooling_table['Teff']], [cooling_table['logg']], cooling_table['Log(Tc)'], kind=interp_kind)




##########################
#Testing how good the interpolation is to the model grid for ages and other stuff.
print('\n\n#############doing test checks (interp to grid points)\n\n')
#size_point=(np.log10(cooling_table['Age']+1)-np.min(np.log10(cooling_table['Age']+1)))+4
#print(np.max(size_point), np.min(size_point))
#print(np.max(cooling_table['Age']),np.min(cooling_table['Age']))
#plt.plot(cooling_table['Teff'],size_point)
#plt.show()
#plt.scatter(cooling_table['Teff'],cooling_table['Mass'],c=size_point)
#plt.contour(cooling_table['Teff'],cooling_table['Mass'],twodim_age,levels=10)
#plt.ylabel('Mass')
#plt.xlabel('Teff')
#plt.show()
#test_mass=0.5
mass_subset_inds=np.where(cooling_table['Mass']==test_mass)
print('Teff shape',cooling_table['Teff'][mass_subset_inds].shape)
#interped_ages=teffm_to_age(cooling_table['Teff'][mass_subset_inds],cooling_table['Mass'][mass_subset_inds])
interped_ages=operate_on_dist(cooling_table['Teff'][mass_subset_inds],cooling_table['Mass'][mass_subset_inds],teffm_to_age)
#interped_ages=teffm_to_age(cooling_table['Teff'][mass_subset_inds], [test_mass,test_mass])
#for i in range(0,interped_ages.shape[0]):
    #print(cooling_table['Teff'][mass_subset_inds][i],interped_ages[i],cooling_table['Age'][mass_subset_inds][i])

plt.scatter(cooling_table['Teff'][mass_subset_inds], cooling_table['Age'][mass_subset_inds],label='grid points')
plt.scatter(cooling_table['Teff'][mass_subset_inds], interped_ages, label='interpolated ages for mass '+str(test_mass))
#plt.scatter(cooling_table['Teff'][mass_subset_inds],mteff_to_age(cooling_table['Mass'][mass_subset_inds],cooling_table['Teff'][mass_subset_inds]),label='interpolated with flipped axes...better be the same as the other one', alpha=0.2)
plt.legend()
plt.ylabel('Cooling Age')
plt.xlabel('Teff')
plt.show()


subsubset=np.where(cooling_table['Teff'][mass_subset_inds]<=11000.)
teff_subset=cooling_table['Teff'][mass_subset_inds][subsubset]
interped_ages_subset=interped_ages[subsubset]
cooling_subset=cooling_table['Age'][mass_subset_inds][subsubset]

plt.plot(cooling_table['Teff'][mass_subset_inds], (interped_ages-cooling_table['Age'][mass_subset_inds])/cooling_table['Age'][mass_subset_inds])
plt.ylabel("Fractional Error of Interpolated Cooling Age compared to Grid Point")
plt.xlabel("Teff")
plt.title('Mass '+str(test_mass)+ ' Fractional Error of Cooling Ages for ' +interp_kind+' interpolation')
#plt.xlim([0,10000])
plt.show()

plt.plot(teff_subset, (interped_ages_subset-cooling_subset)/cooling_subset)
plt.ylabel("Fractional Error of Interpolated Cooling Age compared to Grid Point")
plt.xlabel("Teff")
plt.title('Mass '+str(test_mass)+ ' Fractional Error of Cooling Ages for ' +interp_kind+' interpolation')
plt.show()



plt.plot(cooling_table['Teff'][mass_subset_inds], (interped_ages-cooling_table['Age'][mass_subset_inds])*1e-9)
plt.ylabel("Interpolated Cooling Age minus  Grid Point (Gyr)")
plt.xlabel("Teff")
plt.title('Mass = '+str(test_mass)+ '  Error of Cooling Ages for ' +interp_kind+' interpolation')
#plt.xlim([0,10000])
plt.show()

plt.plot(teff_subset, (interped_ages_subset-cooling_subset)*1e-9)
plt.ylabel("Interpolated Cooling Age minus  Grid Point (Gyr)")
plt.xlabel("Teff")
plt.title('Mass = '+str(test_mass)+ '  Error of Cooling Ages for ' +interp_kind+' interpolation')
#plt.xlim([0,10000])
plt.show()

unique_masses=np.unique(cooling_table['Mass'])
#for test_mass in np.arange(0.2,1.3,0.05):
color_list=['r','g','b','orange','k','pink','purple','brown','olive','cyan','yellow','dodgerblue','tomato','maroon','darkolivegreen','goldenrod','lime','burlywood','aquamarine','k','k','k','k','k','k']

#sub_masses=np.where((unique_masses>=0.4) & (unique_masses <=1.2))
#for test_mass,plot_color  in zip(unique_masses[sub_masses],color_list):
for test_mass,plot_color  in zip(unique_masses,color_list):
    print(test_mass)
    mass_subset_inds=np.where(cooling_table['Mass']==test_mass)
    #print(mass_subset_inds)
    print('Teff shape',cooling_table['Teff'][mass_subset_inds].shape)
    #interped_ages=teffm_to_age(cooling_table['Teff'][mass_subset_inds],cooling_table['Mass'][mass_subset_inds])
    interped_ages=operate_on_dist(cooling_table['Teff'][mass_subset_inds],cooling_table['Mass'][mass_subset_inds],teffm_to_age)
    #interped_ages=teffm_to_age(cooling_table['Teff'][mass_subset_inds], [test_mass,test_mass])
    #for i in range(0,interped_ages.shape[0]):
        #print(cooling_table['Teff'][mass_subset_inds][i],interped_ages[i],cooling_table['Age'][mass_subset_inds][i])

    subsubset=np.where(cooling_table['Teff'][mass_subset_inds]<=11000.)
    teff_subset=np.copy(cooling_table['Teff'][mass_subset_inds][subsubset])
    interped_ages_subset=np.copy(interped_ages[subsubset])
    cooling_subset=np.copy(cooling_table['Age'][mass_subset_inds][subsubset])
    plt.plot(teff_subset, (interped_ages_subset-cooling_subset)*1e-9, label='Mass='+str(test_mass), color=plot_color)
    #plt.plot(teff_subset, cooling_subset*1e-9,linestyle=':', color=plot_color)
    #plt.plot(teff_subset, interped_ages_subset*1e-9,label='Mass='+str(test_mass),color=plot_color)
plt.ylabel("Delta Cooling Age (Interpolated -  Grid Point) (Gyr)")
#plt.ylabel("Cooling Age Grid (Gyr)")
plt.xlabel("Teff (K)")
#plt.title('All masses Error of Cooling Ages for ' +interp_kind+' interpolation')
plt.title('All masses Error of Cooling Ages for ' +interp_kind+' '+cooling_modeler+' ' + atm_type)
#plt.xlim([0,10000])
plt.legend()
plt.grid()
plt.show()


print('\n\n\n########end test checks##########\n\n\n')
#########################

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

#try:
if python_version_number <3.:
    trim_vals, trim_edges, trim_patches= plt.hist(trimmed_total_age_dist, bins=np.arange(0,null_age_val+bin_widths, bin_widths), label='total ages limited to universe', normed=True, alpha=0.2)
#except AttributeError:
elif python_version_number>=3.:
    #plt.show()
    trim_vals, trim_edges, trim_patches= plt.hist(trimmed_total_age_dist, bins=np.arange(0,null_age_val+bin_widths, bin_widths), label='total ages limited to universe', density=True, alpha=0.2)
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
try:
    print("max percentage", np.max(cumprob[inbounds]))
    print("total percentage in distribution", np.max(cumprob))
except ValueError:
    pass


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
print(wd_name)
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

print("\n\n")
print("Median Total Age:", np.median(trimmed_total_age_dist))
print("16th percentile:", np.nanpercentile(trimmed_total_age_dist,16), "difference w/ median:", np.median(trimmed_total_age_dist)-np.nanpercentile(trimmed_total_age_dist,16))
print("84th percentile:", np.nanpercentile(trimmed_total_age_dist,84),"difference w/ median:", np.nanpercentile(trimmed_total_age_dist,84)- np.median(trimmed_total_age_dist))

print("\n\n")
plt.show()

these_bins= np.arange(0,1.4, 0.025)
try:
    plt.hist(bounded_masses, alpha =0.2, label='bounded_masses', bins=these_bins, normed=True)
    plt.hist(target_mass_dist, alpha=0.2, label='target_mass_dist', bins=these_bins, normed=True)
    plt.hist(trimmed_mass_dist, alpha=0.2, label='trimmed_mass_dist', bins=these_bins, normed=True)
except AttributeError:
    plt.show()
    plt.hist(bounded_masses, alpha =0.2, label='bounded_masses', bins=these_bins, density=True)
    plt.hist(target_mass_dist, alpha=0.2, label='target_mass_dist', bins=these_bins, density=True)
    plt.hist(trimmed_mass_dist, alpha=0.2, label='trimmed_mass_dist', bins=these_bins, density=True)
plt.legend()
plt.title(wd_name)

plt.xlabel('Mass (M_sol)')
plt.show()
plt.plot(sort_trim_vals)
plt.show()

try:
    plt.hist(trimmed_teff_dist, normed=True)
except AttributeError:
    plt.show()
    plt.hist(trimmed_teff_dist, density=True)
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

print('\n\n\n#####\n\nclear out the plotting\n\n\n\n\n')
plt.show()

#try:
if python_version_number < 3.:
    plt.hist(target_age_dist, bins=np.arange(0,null_age_val+0.25, 0.25),label='cooling ages', normed=True)
    #plt.hist(ms_age_dist[~np.isnan(ms_age_dist)], bins=50, alpha=0.4, label='ms ages')
    plt.hist(clean_total_age_dist, bins=np.arange(0,null_age_val+0.25, 0.25), alpha=0.5,label= 'logg-> M -> Total Ages', normed=True)
    plt.hist(trimmed_total_age_dist, bins=np.arange(0,null_age_val+bin_widths, bin_widths), label='total ages limited to universe', normed=True, alpha=0.2)
    #plt.hist(total_simon_dist, bins=np.arange(0,21, 0.25), alpha=0.5, label='Simon M -> total ages', normed=True)
    #plt.hist(trimmed_lowz_total_age_dist, bins=np.arange(0,null_age_val+0.1, 0.1), label='lowz total ages limited to universe', normed=True, alpha=0.2)
#except AttributeError:
elif python_version_number>=3.:
    plt.show()
    plt.hist(target_age_dist, bins=np.arange(0,null_age_val+0.25, 0.25),label='cooling ages', density=True)
    #plt.hist(ms_age_dist[~np.isnan(ms_age_dist)], bins=50, alpha=0.4, label='ms ages')
    plt.hist(clean_total_age_dist, bins=np.arange(0,null_age_val+0.25, 0.25), alpha=0.5,label= 'logg-> M -> Total Ages', density=True)
    plt.hist(trimmed_total_age_dist, bins=np.arange(0,null_age_val+bin_widths, bin_widths), label='total ages limited to universe', density=True, alpha=0.2)
    #plt.hist(total_simon_dist, bins=np.arange(0,21, 0.25), alpha=0.5, label='Simon M -> total ages', normed=True)
    #plt.hist(trimmed_lowz_total_age_dist, bins=np.arange(0,null_age_val+0.1, 0.1), label='lowz total ages limited to universe', normed=True, alpha=0.2)
else:
    print('python_version_number',python_version_number, "somehow couldn't be <3. or >=3., so something is very wrong")
    sys.exit()
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
#try:
if python_version_number < 3.:
    plt.hist(clean_log_dist, bins=np.linspace(0,null_age_val+1,1000), normed=True, color='g', alpha=0.5, label='Total Ages')
#except AttributeError:
elif python_version_number>=3.:
    plt.hist(clean_log_dist, bins=np.linspace(0,null_age_val+1,1000), density=True, color='g', alpha=0.5, label='Total Ages')
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

#approx_inds= np.where((cooling_table['Teff']< 4000) & (cooling_table['Teff']> 3500))
#approx_masses= cooling_table['Mass'][approx_inds]
#approx_ages=cooling_table['Age'][approx_inds]*1e-9

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
#try:
if python_version_number < 3.:
    plt.hist(target_mass_dist, normed=True, label=wd_name+' MC Masses', color='k')
#except:
elif python_version_number>=3.:
    mass_bins=np.arange(0.1,1.4,0.01)
    plt.hist(target_mass_dist, density=True, stacked=True,label=wd_name+' MC Masses', color='k',bins=mass_bins)
#plt.scatter(0.56,  teffm_to_age(target_teff, 0.56)*1e-9+get_ms_lifetime(0.56), label='M=0.56 at teff'+str(target_teff))
plt.xlabel(r'$M_{wd}$ $(M_{\odot})$')
#plt.ylabel('MS lifetime (Gyr)')
plt.ylabel('Age (Gyr)')
plt.legend()
#plt.yscale('log')
plt.ylim(0,15)
#plt.title(wd_name)
plt.show()

approx_inds= np.where((cooling_table['Teff']< (target_teff+250)) & (cooling_table['Teff']> (target_teff-250)))
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

if atm_type=='thinH':
    alt_cooling_model_file='COModel_ThinH.csv'
elif atm_type=='thickH':
        alt_cooling_model_file='COModel_ThickH.csv'
alt_cooling_model_file= cp.ref_dir+'WD_cooling_models/'+alt_cooling_model_file
alt_cooling_table=Table.read(alt_cooling_model_file)


bed_inds=np.where(cooling_table['Mass'] ==0.5)
font_inds=np.where(alt_cooling_table['Mass']==0.5)

plt.scatter(cooling_table['Teff'][bed_inds], cooling_table['Age'][bed_inds]*1e-9, label='cooling models'+cooling_model_file.split('/')[-1],s=6,alpha=0.5)
plt.scatter(alt_cooling_table['Teff'][font_inds],alt_cooling_table['Age'][font_inds]*1e-9,label='Fontaine2001 '+atm_type,s=6,alpha=0.5)
#plt.plot(target_teff, target_logg, marker='*', color='r', markersize=12)
plt.xlabel('Teff (K)')
plt.ylabel('Age (Gyr)')
plt.legend()
plt.title(cooling_model_file.split('/')[-1]+' for 0.5 M_sol WD')
plt.show()

plt.scatter(cooling_table['Teff'], cooling_table['logg'], label='cooling models'+cooling_model_file.split('/')[-1],s=6,alpha=0.5)
plt.scatter(alt_cooling_table['Teff'],alt_cooling_table['logg'],label='Fontaine2001 '+atm_type,s=6,alpha=0.5)
plt.plot(target_teff, target_logg, marker='*', color='r', markersize=12)
plt.xlabel('Teff (K)')
plt.ylabel('log(g)')
plt.legend()
plt.title(cooling_model_file.split('/')[-1])
plt.show()

plt.scatter(np.log10(cooling_table['Teff']), np.log10(cooling_table['Lum']), label='cooling models '+cooling_model_file.split('/')[-1],s=6,alpha=0.5)
plt.scatter(np.log10(alt_cooling_table['Teff']),np.log10(alt_cooling_table['Lum']),label='Fontaine2001 '+atm_type,s=6,alpha=0.5)
#plt.plot(target_teff, target_logg, marker='*', color='r', markersize=12)
plt.xlabel('log10(Teff (K))')
plt.ylabel('log10(Luminosity)')
plt.legend()
plt.title(cooling_model_file.split('/')[-1])
plt.show()

plt.scatter(cooling_table['Teff'], cooling_table['rad'], label='cooling models'+cooling_model_file.split('/')[-1],s=6,alpha=0.5)
plt.scatter(alt_cooling_table['Teff'],alt_cooling_table['rad'],label='Fontaine2001 '+atm_type,s=6,alpha=0.5)
#plt.plot(target_teff, target_logg, marker='*', color='r', markersize=12)
plt.xlabel('Teff (K)')
plt.ylabel('Radius')
plt.legend()
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
#plt.scatter(cooling_table['Teff'][given_inds], teffm_to_age(cooling_table['Teff'][given_inds], given_mass).T, label='interpolated vals', color='r')
plt.scatter(cooling_table['Teff'][given_inds], operate_on_dist(cooling_table['Teff'][given_inds], np.full((cooling_table['Teff'][given_inds]).shape, given_mass),teffm_to_age).T, label='interpolated vals Transposed', color='r')
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
plt.plot(teff_vals, loggteff_to_logTc(teff_vals, target_logg_array).T, label='log(Tc) transposed')
plt.scatter(cooling_table['Teff'][allowed_inds], cooling_table['Log(Tc)'][allowed_inds], label='cooling table direct values for '+str(given_target_mass)+ 'M')
plt.plot(target_teff, target_logTc, marker='*', label='target')
plt.xlabel('Teff')

plt.ylabel('log(Tc)')
plt.legend()
plt.show()


logg_vals= np.linspace(7., 8.5, 100)
allowed_inds= np.where(cooling_table['Mass']==given_target_mass)
plt.plot(logg_vals, loggteff_to_logTc(target_teff, logg_vals), label='log(Tc)')
plt.plot(logg_vals, loggteff_to_logTc(target_teff, logg_vals).T, label='transposed log(Tc)')

plt.scatter(cooling_table['logg'][allowed_inds], cooling_table['Log(Tc)'][allowed_inds], label='cooling table direct values for '+str(given_target_mass)+ 'M')
plt.plot(target_logg, target_logTc, marker='*', label='target')
plt.xlabel('logg')

plt.ylabel('log(Tc)')
plt.legend()
plt.show()


print('t_ms for 0.56 M_wd', get_ms_lifetime(0.56,method='Hurley'))





