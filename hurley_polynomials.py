"""
Created by Ben Kaiser (UNC-Chapel Hill) 2020-01-06

This is essentially a Python implementation of some of the routines/equations presented in Hurley, Pols, and Tout 2000 with an emphasis (and probably entirety of implementation) focused on being used for white dwarf progenitor MS lifetimes.


"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy import constants as const
from astropy.time import Time
from astropy.table import Table
import scipy.interpolate as scinterp



import cal_params as cp
from a_coeffs import a_coeffs
