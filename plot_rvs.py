"""
I'm changing the radial velocity plot  to read in an astropy file and plot that instead.
"""

import sys
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, Column

plt.rc('font', size =18)
plt.rc('lines', markersize=8)
plotfile = glob(sys.argv[1])[0]
#all_array = np.genfromtxt(plotfile[0], names = True, delimiter = ',')
#print all_array.dtype.names
#bmjd_array = all_array['BMJD_TDB']
#H_delta= all_array['H_delta']
#H_gamma = all_array['H_gamma']
#H_beta = all_array['H_beta']
#H_delta_s = all_array["H_delta_s"]
#H_gamma_s = all_array['H_gamma_s']
#H_beta_s = all_array['H_beta_s']
#print all_array
##mean_rv = np.mean(all_array[:,1 :], axis = 1)
##std_dev = np.std(all_array[:,1 :], axis =1)
##print mean_rv
##plt.plot(bmjd_array, H_delta, label = r"H-$\delta$", linestyle = 'none', marker = '*')
##plt.plot(bmjd_array, H_gamma, label = r"H-$\gamma$", linestyle = 'none', marker = '*')
##plt.plot(bmjd_array, H_beta, label = r"H-$\beta$", linestyle = 'none', marker = '*')
#plt.errorbar(bmjd_array, H_delta, H_delta_s, label = r"H-$\delta$", linestyle = 'none', marker = '*')
#plt.errorbar(bmjd_array, H_gamma, H_gamma_s, label = r"H-$\gamma$", linestyle = 'none', marker = '*')
#plt.errorbar(bmjd_array, H_beta, H_beta_s, label = r"H-$\beta$", linestyle = 'none', marker = '*')
##plt.plot(mjd_array, mean_rv, label = 'Mean RV')
##plt.errorbar(bmjd_array, mean_rv, std_dev, label = 'Mean RV')
#plt.title('')
#plt.ylabel('Radial Velocity (km/s)')
#plt.xlabel('BMJD')
#plt.legend()
#plt.show()


input_table= Table.read(plotfile, format='ascii.csv')
input_table.pprint()
plt.errorbar(input_table['bmjd_tdb'], input_table['rv'], yerr=input_table['rv_error'], marker='o', linestyle='none')
plt.xlabel('BMDJ_TDB')
plt.ylabel('Radial Velocity (km/s)')
plt.show()
