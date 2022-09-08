"""
Created by Ben Kaiser (UNC-Chapel Hill) 2022-09-08

Read in the "recommended" Milky Way dataset from the SAGA database of extremely metal poor stars and 
produce plots of various abundances using subsets of the sample.

This script will mostly be called by other scripts most likely in order to produce the backgrounds of plots.

Database URL: http://sagadatabase.jp

I'm going to write this to work in Python3 as its main goal. If it happens to be compatible with Python 2.7 
that's a happy accident.


"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.table import join as ATjoin
from astropy.table import hstack as AThstack
import periodictable as pt
from astropy import units as u
from astropy import constants as const

import cal_params as cp
