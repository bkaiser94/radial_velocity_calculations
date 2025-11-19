"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-07-08.

This is basically the mirror of saga_plotting.py, but this one will work with data from the JINAbase of metal poor stars.




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


jina_file='JINAbase_20250708_metalpoor.txt'


lodders_file='Lodders2020_solarsystem_abundances.csv'


jina_file=cp.abundance_dir+jina_file
lodders_file=cp.abundance_dir+lodders_file


jina_table=Table.read(jina_file,format='ascii')


jina_table.pprint()
















