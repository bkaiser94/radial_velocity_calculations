"""
Created by Ben Kaiser 2019-06-19 (UNC-Chapel Hill)

Presumably this will be a one-time use script that will convert Stritzinger et al. 2005's extinction curve to be an
astropy table that has headers, which will make it less opaque as to its meaning and make it safer to read-in in the 
future without causing confusion.

"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, QTable


import cal_params as cp

filename= cp.ref_dir+'extinction/Stritzinger_2005_extinction_curve.txt'
output_filename=cp.ref_dir+'extinction/Stritzinger_2005_extinction_curve.csv'
all_array = np.genfromtxt(filename)

#print(all_array)
#print(all_array.shape)

waves= all_array[0]
extinction_mag= all_array[1]

plt.plot(waves, extinction_mag)
plt.xlabel(r'Wavelength($\AA$)')
plt.ylabel('Extinction (mag/airmass)')
plt.show()

output_table= Table(all_array.T, names=('lambda', 'extinction'))

output_table.pprint()




output_table.write(output_filename, format='ascii.csv', overwrite=True)
