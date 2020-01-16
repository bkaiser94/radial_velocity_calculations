"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-12-12


Import the diffusion timescales and other arrays we need from evolution_pulled.py


"""
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

from evolution_pulled import ttDB, ggDB, tauDB, elemnameDB


color_wheel=['b', 'r', 'g', 'k']
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

#for grav,subtau in zip(ggDB_array,tauDB_array):
    #li_taus=  subtau[li_index][:]
    #ca_taus= subtau[ca_index][:]
    #na_taus=subtau[na_index][:]
    #print(li_taus.shape)
    #print(subtau[:][0].shape)
    #print(subtau.shape)
    #plt.plot(np.log10(ttDB_array),li_taus, label='Li logg='+str(grav), marker='o')
    #plt.plot(np.log10(ttDB_array),ca_taus, label='Ca logg='+str(grav), marker='o')
    #plt.plot(np.log10(ttDB_array),na_taus, label='Na logg='+str(grav), marker='o')
#plt.legend(loc='best')
#plt.xlabel('log10(Teff)')
#plt.ylabel('tau_Li')
#plt.title("'DB' diffusion timescales from MWDD/Fontaine et al. 2015")
#plt.show()


for grav,subtau, color  in zip(ggDB_array,tauDB_array, color_wheel):
    li_taus=  subtau[li_index][:]
    ca_taus= subtau[ca_index][:]
    na_taus=subtau[na_index][:]
    fe_taus=subtau[fe_index][:]
    mg_taus=subtau[mg_index][:]
    print(li_taus.shape)
    print(subtau[:][0].shape)
    print(subtau.shape)
    plt.plot(ttDB_array,li_taus-ca_taus, label='Li-Ca logg='+str(grav), marker='s', color=color)
    plt.plot(ttDB_array,na_taus-ca_taus, label='Na-Ca logg='+str(grav), marker='^', color=color)
    plt.plot(ttDB_array,mg_taus-ca_taus, label='Mg-Ca logg='+str(grav), marker='o', color=color)
    #plt.plot(ttDB_array,na_taus-fe_taus, label='Fe-Ca logg='+str(grav), marker='^', color=color)
plt.legend(loc='best')
plt.xlabel('Teff')
#plt.xscale('log')
plt.ylabel(r'$\tau_{X}-\tau_{Ca}$')
plt.title("'DB' diffusion timescales from MWDD/Fontaine et al. 2015")
plt.show()

print(tauDB_array.shape)
print(tauDB_array[0][0][:].shape)
print(tauDB_array[0][:,0].shape)
plt.plot(tauDB_array[0][:,0])
plt.plot(tauDB_array[3][:,-1])
plt.show()
