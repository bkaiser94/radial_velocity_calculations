"""
Created by Ben Kaiser (UNC-Chapel Hill) 2022-03-18


"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const
from astropy import convolution as conv
from astropy.table import Table, Column
import scipy.interpolate as scinterp
import time

input_filename='Moss_2022_WDMS_age_metallicity.csv'

full_table=Table.read(input_filename)

plt.errorbar(full_table['total_age'],full_table['fe/h'],xerr=full_table['total_age_err'],yerr=full_table['fe/h_err'],linestyle='None',marker='o')
plt.ylabel('[Fe/H]')
plt.xlabel('Total Age (Gyr)')
plt.title('Moss et al. 2022 subsample of 38 WD-MS binaries with MS spectroscopic [Fe/H]')
#plt.axhline(y=0, linestyle=':', color='k')
plt.grid()
plt.show()

def plot_vals(string1, string2):
    
    
    plt.title('Moss et al. 2022 subsample of 38 WD-MS binaries with MS spectroscopic [Fe/H]')
    plt.grid()
    plt.show()
    return
