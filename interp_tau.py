"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-12-12


Import the diffusion timescales and other arrays we need from evolution_pulled.py


"""
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
import time

from evolution_pulled import ttDB, ggDB, tauDB, elemnameDB

start= time.time()


color_wheel=['b', 'r', 'g', 'k']
teff_max=10000. #max teff to use in the extrapolation of tau values.


tauDB_array= np.array(tauDB)
ttDB_array= np.array(ttDB)
ggDB_array=np.array(ggDB)
DBels= np.array(elemnameDB)

teffs= np.ones(tauDB_array.shape)
#teffs=ttDB_array*teffs
loggs= np.ones(tauDB_array.shape)
#loggs=ggDB_array*loggs
print(tauDB_array.shape)
tauDB_array=np.transpose(tauDB_array, (1,2,0))
print(tauDB_array.shape)

li_index= np.where(DBels=='Li')[0][0]
ca_index=np.where(DBels=='Ca')[0][0]
na_index=np.where(DBels=='Na')[0][0]
fe_index=np.where(DBels=='Fe')[0][0]
mg_index=np.where(DBels=='Mg')[0][0]
k_index=np.where(DBels=='K ')[0][0]

for grav,subtau in zip(ggDB_array,tauDB_array):
    li_taus=  subtau[li_index][:]
    ca_taus= subtau[ca_index][:]
    na_taus=subtau[na_index][:]
    fe_taus=subtau[fe_index][:]
    k_taus=subtau[k_index][:]
    print(li_taus.shape)
    print(subtau[:][0].shape)
    print(subtau.shape)
    #plt.plot(np.log10(ttDB_array),li_taus, label='Li logg='+str(grav), marker='s')
    #plt.plot(np.log10(ttDB_array),ca_taus, label='Ca logg='+str(grav), marker='x')
    #plt.plot(np.log10(ttDB_array),na_taus, label='Na logg='+str(grav), marker='^')
    #plt.plot(np.log10(ttDB_array),k_taus, label='K logg='+str(grav), marker='+')
plt.legend(loc='best')
plt.xlabel('log10(Teff)')
plt.ylabel(r'log10($\tau$)')
plt.title("'DB' diffusion timescales from MWDD/Fontaine et al. 2015")
plt.show()


for grav,subtau, color  in zip(ggDB_array,tauDB_array, color_wheel):
    li_taus=  subtau[li_index][:]
    ca_taus= subtau[ca_index][:]
    na_taus=subtau[na_index][:]
    fe_taus=subtau[fe_index][:]
    mg_taus=subtau[mg_index][:]
    k_taus=subtau[k_index][:]
    #print(li_taus.shape)
    #print(subtau[:][0].shape)
    #print(subtau.shape)
    #plt.plot(ttDB_array,li_taus-ca_taus, label='Li-Ca logg='+str(grav), marker='s', color=color)
    #plt.plot(ttDB_array,na_taus-ca_taus, label='Na-Ca logg='+str(grav), marker='^', color=color)
    ##plt.plot(ttDB_array,mg_taus-ca_taus, label='Mg-Ca logg='+str(grav), marker='o', color=color)
    #plt.plot(ttDB_array,fe_taus-ca_taus, label='Fe-Ca logg='+str(grav), marker='<', color=color)
    #plt.plot(ttDB_array,k_taus-ca_taus, label='K-Ca logg='+str(grav), marker='+', color=color)
    ##plt.plot(ttDB_array,na_taus-fe_taus, label='Fe-Ca logg='+str(grav), marker='^', color=color)
    
    
    #plt.plot(np.log10(ttDB_array),fe_taus-ca_taus, label='Fe-Ca logg='+str(grav), marker='<', color=color)
    #plt.plot(np.log10(ttDB_array),k_taus-ca_taus, label='K-Ca logg='+str(grav), marker='+', color=color)
#plt.legend(loc='best')
#plt.xlabel('Teff')
#plt.xscale('log')
#plt.ylabel(r'$\tau_{X}-\tau_{Ca}$')
#plt.title("'DB' diffusion timescales from MWDD/Fontaine et al. 2015")
#plt.show()

#print(tauDB_array.shape)
#print(tauDB_array[0][0][:].shape)
#print(tauDB_array[0][:,0].shape)
#plt.plot(tauDB_array[0][:,0])
#plt.plot(tauDB_array[3][:,-1])
#plt.show()

#print('li_taus.shape', li_taus.shape)
#print('tauDB_array.shape',tauDB_array.shape)

#print('ttDB_array.shape', ttDB_array.shape)
#print('ggDB_array.shape', ggDB_array.shape)

def get_logg_index(input_logg):
    return np.where(ggDB_array==input_logg)

def make_DB_tau_table(input_logg):
    #print('logg:', input_logg)
    subtau=tauDB_array[get_logg_index(input_logg)][0]
    #print('subtau.shape', subtau.shape)
    li_taus=  subtau[li_index][:]
    ca_taus= subtau[ca_index][:]
    na_taus=subtau[na_index][:]
    fe_taus=subtau[fe_index][:]
    mg_taus=subtau[mg_index][:]
    k_taus=subtau[k_index][:]
    full_array=np.vstack([ttDB_array, li_taus, ca_taus, na_taus, fe_taus, mg_taus, k_taus]).T
    #print('full_array.shape', full_array.shape)
    db_table=Table(full_array, names=('teff', 'Li', 'Ca', 'Na', 'Fe', 'Mg', 'K'))
    return db_table


#tau_table= make_DB_tau_table(8.0)
#tau_table.pprint()

#for  el in ['Li', 'Ca', 'Na', 'Fe', 'Mg', 'K']:
    ##plt.plot(np.log10(tau_table['teff']),tau_table[el], label=el)
    #plt.plot(tau_table['teff'],tau_table[el], label=el)
#plt.legend()
#plt.ylabel('tau')
#plt.xlim(6000, teff_max)
##plt.xlim(np.log10(6000.),np.log10( teff_max))

#plt.show()

def extrapolate_tau_single_logg( element,  input_logg= 8.0, teff_max=teff_max, plot_all=False):
    db_table= make_DB_tau_table(input_logg)
    valid_teff_inds= np.where(db_table['teff']<teff_max)
    valid_teffs= db_table[valid_teff_inds]['teff']
    valid_taus= db_table[valid_teff_inds][element]
    log_valid_teffs= np.log10(valid_teffs)
    poly_coeffs= np.polyfit(log_valid_teffs, valid_taus, 1) #it's going to only be a linear fit to the log-log space, tau is already logged
    if plot_all:
            plt.scatter(log_valid_teffs, valid_taus, label=element+ ' DB grid, logg='+str(input_logg))
    return poly_coeffs

def make_extrap_cross_logg(element, teff_max=teff_max):
    logg_coeffs=[]
    for grid_logg in ggDB_array:
        coeff_single=extrapolate_tau_single_logg(element, input_logg=grid_logg, teff_max=teff_max)
        logg_coeffs.append(coeff_single)
    return logg_coeffs

def extrapolate_single_el_tau(target_teff, element,  input_logg= 8.0, teff_max=teff_max, plot_all=False):
    #db_table= make_DB_tau_table(input_logg)
    #valid_teff_inds= np.where(db_table['teff']<teff_max)
    #valid_teffs= db_table[valid_teff_inds]['teff']
    #valid_taus= db_table[valid_teff_inds][element]
    #log_valid_teffs= np.log10(valid_teffs)
    #poly_coeffs= np.polyfit(log_valid_teffs, valid_taus, 1) #it's going to only be a linear fit to the log-log space, tau is already logged
    poly_coeffs=extrapolate_tau_single_logg(element,  input_logg= input_logg, teff_max=teff_max, plot_all=plot_all)
    target_tau= np.polyval(poly_coeffs, np.log10(target_teff))
    test_teffs=np.log10(np.linspace(3000,teff_max, 100))
    if plot_all:
        plt.plot(test_teffs, np.polyval(poly_coeffs, test_teffs), label='test extrapolation'+str(input_logg))
        #plt.scatter(log_valid_teffs, valid_taus, label='grid')
        plt.scatter(np.log10(target_teff), target_tau, label='target', color='r')
        plt.legend(loc='best')
        plt.xlabel('log10(Teff)')
        plt.ylabel('Tau')
        plt.title(element)
        plt.show()
    else:
        pass
    #plt.show()
    return target_tau
    #need the teff values to also be logged... 
    
def extrapolate_tau_x_logg(target_teff, target_logg, element):
    input_vals_array=True
    poly_coeffs= make_extrap_cross_logg(element)
    try:
        print(target_teff[0], target_logg[0])
    except TypeError:
        input_vals_array=False
    except IndexError:
        input_vals_array=False
        #for some crazy reason attempting to index a float value in the terminal python script is a 
        #TypeError, but doing the exact same thing when calling from a jupyter notebook is an IndexError. 
        #Can't make this shit up. Maybe they're slightly different python version?
        #print("Not an array input, which is ok.")
    if input_vals_array:
        tau_list=[]
        #loop_count=0
        for teff, logg in zip(target_teff, target_logg):
            #if loop_count%1000==0:
                #print('loop round done', loop_count)
            #loop_count+=1
            logg_diffs= np.abs(logg-ggDB_array)
            small_args=np.argsort(logg_diffs)[:2]
            bound_taus=[]
            for arg in small_args:
               
                tau_teff= np.polyval(poly_coeffs[arg], np.log10(teff))
                bound_taus.append(tau_teff)
            cross_coeffs=np.polyfit(ggDB_array[small_args], bound_taus, 1)
            single_tau= np.polyval(cross_coeffs, logg)
            tau_list.append(single_tau)
        return np.array(tau_list)
    else:
        logg_diffs= np.abs(target_logg-ggDB_array)
        small_args=np.argsort(logg_diffs)[:2]
        bound_taus=[]
        for arg in small_args:
            
            tau_teff= np.polyval(poly_coeffs[arg], np.log10(target_teff))
            bound_taus.append(tau_teff)
        cross_coeffs=np.polyfit(ggDB_array[small_args], bound_taus, 1)
        single_tau= np.polyval(cross_coeffs, target_logg)
        return single_tau
    return

if __name__ == '__main__':
    wd_name='GaiaJ1644-0449'
    target_logg=7.77
    target_logg_err= 0.23
    target_teff= 3830.
    target_teff_err= 230.
    
    
    #wd_name='WDJ2356-209'
    #target_logg=7.98
    #target_logg_err=0.07
    #target_teff= 4040. #K
    #target_teff_err=110.
    
    #wd_name='SDSSJ1330+6435'
    #target_logg= 8.26
    #target_logg_err=0.15
    #target_teff= 4310. #K
    #target_teff_err=190
    
    ######from Coutu et al. 2019 retrieved from MWDD
    #wd_name='SDSSJ1636+1619'
    #target_logg= 8.08
    #target_logg_err=0.07 #not provided in MWDD, but I assumed it was similar to error for J2356
    ##target_teff= 4416. #K #Coutu values
    ##target_teff_err=66.
    #target_teff= 4410. #K #Hollands values
    #target_teff_err=150.
    
    n_points=1e5
    logg_dist=np.random.normal(loc=target_logg, scale=target_logg_err, size=n_points)
    teff_dist= np.random.normal(loc=target_teff, scale=target_teff_err, size=n_points)
    tau_li_dist=extrapolate_tau_x_logg(teff_dist, logg_dist, 'Li')
    tau_ca_dist=extrapolate_tau_x_logg(teff_dist, logg_dist, 'Ca')
    tau_na_dist=extrapolate_tau_x_logg(teff_dist, logg_dist, 'Na')
    tau_fe_dist=extrapolate_tau_x_logg(teff_dist, logg_dist, 'Fe')
    tau_k_dist=extrapolate_tau_x_logg(teff_dist, logg_dist, 'K')
    stop= time.time()
    print(stop-start)
    print((stop-start)/60.)
    print('Na', extrapolate_tau_x_logg(target_teff, target_logg, 'Na'))
    print('Na', extrapolate_single_el_tau(target_teff, 'Na', plot_all=True))
    print('Ca', extrapolate_tau_x_logg(target_teff, target_logg, 'Ca'))
    print('K', extrapolate_tau_x_logg(target_teff, target_logg, 'K'))
    print('K', extrapolate_single_el_tau(target_teff, 'K', plot_all=True))
    print('Ca', extrapolate_single_el_tau(target_teff, 'Ca', plot_all=True))
    print('Li', extrapolate_tau_x_logg(target_teff, target_logg, 'Li'))
    print('Li', extrapolate_single_el_tau(target_teff, 'Li', plot_all=True))
    print('Fe', extrapolate_tau_x_logg(target_teff, target_logg, 'Fe'))
    plt.show()
    #bins=np.arange(4.5,7., 0.05)
    #bins=np.arange(4.0,9.0, 0.05)
    bins=np.arange(3.,9.0, 0.05)
    #plt.hist(tau_li_dist, label='Li', alpha=0.2, bins=bins)
    plt.hist(tau_ca_dist, label='Ca', alpha=0.2, bins=bins)
    #plt.hist(tau_na_dist, label='Na', alpha=0.2, bins=bins)
    plt.hist(tau_fe_dist, label='Fe', alpha=0.2,bins=bins)
    plt.hist(tau_k_dist, label='K', alpha=0.2,bins=bins)
    plt.legend()
    plt.title(wd_name)
    plt.xlabel(r'$\tau_{el}$')
    plt.show()
    plt.hist(logg_dist)
    plt.axvline(x=7.5, linestyle='--')
    plt.axvline(x=9.0, linestyle='--')
    plt.show()
    plt.scatter(teff_dist, tau_li_dist)
    plt.scatter(teff_dist, tau_ca_dist)
    plt.show()
    plt.scatter(logg_dist, tau_li_dist)
    plt.scatter(logg_dist, tau_ca_dist)
    plt.show()
    plt.scatter(logg_dist, np.array(tau_li_dist)-np.array(tau_ca_dist))
    plt.title('tau_Li - tau_ca '+wd_name)
    plt.show()
    plt.scatter(logg_dist, np.array(tau_na_dist)-np.array(tau_ca_dist))
    plt.title('tau_na - tau_ca '+wd_name)
    plt.xlabel('log g')
    plt.show()
    plt.scatter(logg_dist, np.array(tau_k_dist)-np.array(tau_ca_dist))
    plt.title('tau_k - tau_ca '+wd_name)
    plt.xlabel('log g')
    plt.show()
    plt.hist(np.array(tau_li_dist)-np.array(tau_ca_dist))
    plt.show()
    target_tau_li=extrapolate_single_el_tau(3830., 'Li')
    target_tau_li=extrapolate_single_el_tau(3830., 'Li', input_logg=7.5)
    target_tau_ca=extrapolate_single_el_tau(3830., 'Ca')
    print('tau_li-tau_ca', target_tau_li-target_tau_ca)
    plt.show()
    for grav in ggDB_array:
        extrapolate_single_el_tau(3830., 'Li', input_logg=grav)
    plt.show()
    for grav in ggDB_array:
        extrapolate_single_el_tau(3830., 'Ca', input_logg=grav)
    plt.show()

    for grav in ggDB_array:
        extrapolate_single_el_tau(3830., 'Na', input_logg=grav)
    plt.show()

    dist_tau_li= extrapolate_single_el_tau(np.random.normal(loc=3830., scale=230., size=1000.),'Li')
    plt.hist(dist_tau_li)
    plt.show()
