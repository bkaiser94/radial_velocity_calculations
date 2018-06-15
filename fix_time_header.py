"""
This fixes the GPS-synched time header (well technically the year of it) because it was set to 2019 for some reason when the data was taken in 2018. This edits the original files in listSpec, so it should not be run in the RAW directory (not that anything ever really should be).

This is technically actually step1 in the step-by-step process, but it's not really a reduction step. It's a revision thing.

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


header_name = 'OPENDATE'
speclistname = 'listSpec'
correct_year = 2018

speclist= np.genfromtxt(speclistname, dtype = 'str')


for filename in speclist:
    i = fits.open(filename)
    header= i[0].header
    header_time = header['OPENDATE']
    print header_time
    split_header = header_time.split('-')
    year = int(split_header[0])
    if year != correct_year:
        year = correct_year
        split_header[0]= str(year)
        print '-'.join(split_header)
        new_header = '-'.join(split_header)
        header['OPENDATE']= new_header
        header.append(card= ('typocorr', 'True', 'Fixed incorrect year in OPENDATE'))
        i.writeto(filename, overwrite = True)
    else:
        print "correct time originally"
