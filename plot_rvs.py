import sys
from glob import glob
import matplotlib.pyplot as plt
import numpy as np

plotfile = glob(sys.argv[1])
all_array = np.genfromtxt(plotfile[0], names = True, delimiter = ',')
print all_array.dtype.names
bmjd_array = all_array['BMJD_TDB']
H_delta= all_array['H_delta']
H_gamma = all_array['H_gamma']
H_beta = all_array['H_beta']

print all_array
#mean_rv = np.mean(all_array[:,1 :], axis = 1)
#std_dev = np.std(all_array[:,1 :], axis =1)
#print mean_rv
plt.plot(bmjd_array, H_delta, label = r"H-$\delta$", linestyle = 'none', marker = '*')
plt.plot(bmjd_array, H_gamma, label = r"H-$\gamma$", linestyle = 'none', marker = '*')
plt.plot(bmjd_array, H_beta, label = r"H-$\beta$", linestyle = 'none', marker = '*')
#plt.plot(mjd_array, mean_rv, label = 'Mean RV')
#plt.errorbar(bmjd_array, mean_rv, std_dev, label = 'Mean RV')
plt.title('')
plt.ylabel('Radial Velocity (km/s)')
plt.xlabel('BMJD')
plt.legend()
plt.show()
