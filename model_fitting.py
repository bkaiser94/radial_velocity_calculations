"""
this script should open a model file (or all of them I suppose more accurately, and step through them.

"""
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
import cosmics
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const


import wdatmos

target_list_name = 'listFWCTB'
target_list = np.genfromtxt(target_list_name, dtype = 'str')

scaling_range = [4600,4650]

target_file = target_list[0]
print target_file
i= fits.open(target_file)
header = fits.getheader(target_file)
target_waves= i[0].data
target_flux = i[1].data

def chi_squared(observed, actual):
    return (observed - actual)**2/actual

#David's instructions for loading the model
wd=wdatmos.wdmodel(filename='ELM.hdf5')
teff = 8000
logg = 6.25
####3

teff_array = np.arange(6000, 15000, 250)
logg_array = np.arange(3.75, 6.5, 0.25)
print wd.Teffs[0]
print wd.loggs
model=wd(Teff= teff, logg =logg)
plt.plot(model['w'], model['flux'])
plt.show()
#print model
#for teff in teff_array:
    #for logg in logg_array:
        #model = wd(Teff = teff , logg = logg)
        ##print model['w'][0]
        #if model != None:
            #model_num +=1


#####

model_waves = model['w']
model_flux = model['flux'] #since we'll be arbitrarily-ish scaling this it won't work.

target_vals = np.where((target_waves > scaling_range[0]))
scale_factor =np.mean( target_flux[scaling_range[0]:scaling_range[1]])/np.mean(model_flux[scaling_range[0]:scaling_range[1]])
print scale_factor
scale_model_flux = model_flux* scale_factor

print scale_model_flux.mean()
print target_flux.mean()
plt.plot(model_waves, scale_model_flux, label = 'model'+str(teff) + ' ' + str(logg))
plt.plot(target_waves, target_flux, label = 'Target')
plt.legend()
plt.xlabel('Angstroms')
plt.ylabel('Flux in cgs 10**-16')
plt.show()
