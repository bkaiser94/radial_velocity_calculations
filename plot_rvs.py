import sys
from glob import glob
import matplotlib.pyplot as plt
import numpy as np

plotfile = glob(sys.argv[1])
all_array = np.genfromtxt(plotfile[0]).T
mjd_array = all_array[0]
H_8 = all_array[1]
H_9 = all_array[2]
H_10 = all_array[3]

print all_array
mean_rv = np.mean(all_array[1:, :], axis = 0)
std_dev = np.std(all_array[1:, :], axis =0)
print mean_rv
plt.plot(mjd_array, H_8, label = "H8", linestyle = 'none', marker = '*')
plt.plot(mjd_array, H_9, label = "H9", linestyle = 'none', marker = '*')
plt.plot(mjd_array, H_10, label = "H10", linestyle = 'none', marker = '*')
#plt.plot(mjd_array, mean_rv, label = 'Mean RV')
plt.errorbar(mjd_array, mean_rv, std_dev, label = 'Mean RV')
plt.title('')
plt.ylabel('Radial Velocity (km/s)')
plt.xlabel('MJD')
plt.legend()
plt.show()
