
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
