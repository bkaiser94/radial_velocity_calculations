"""
This fixes the GPS-synched time header (well technically the year of it) because it was set to 2019 for some reason when the data was taken in 2018. This edits the original files in listSpec, so it should not be run in the RAW directory (not that anything ever really should be).

This is technically actually step1 in the step-by-step process, but it's not really a reduction step. It's a revision thing.

We're having time issues again in 2024 because they disconnected SOAR from the internet, so I'm updating this to python3 (adding parentheses for print statements) and changing the correct year.

"""
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const


header_name = 'OPENDATE'
speclistname = 'listSpec'
#correct_year = 2018
correct_year=2024
header_list=[
    'DATE-OBS',
    'DATE',
    'OPENDATE'
    ]


speclist= np.genfromtxt(speclistname, dtype = 'str')


for filename in speclist:
    i = fits.open(filename)
    header= i[0].header
    typocorr_added=False
    for header_name in header_list:
        header_time = header[header_name]
        print(header_name, header_time)
        split_header = header_time.split('-')
        year = int(split_header[0])
        if year != correct_year:
            if ((year %4==0)^(correct_year%4==0)):
                if year%4==0:
                    day=int(split_header[2])
                    day=day+1
                    split_header[2]=str(day)
                    print('year was leap year')
                    print("didn't get around to actually fixing this because it would require me to fix months too most likely")
                elif current_year%4==0:
                    day=int(split_header[2])
                    day=day-1
                    split_header[2]=str(day)
                    print('current_year was leap year')
                    print("didn't get around to actually fixing this because it would require me to fix months too most likely")
                    pass
            year = correct_year
            split_header[0]= str(year)
            print('-'.join(split_header))
            new_header = '-'.join(split_header)
            header[header_name]= new_header
            if typocorr_added:
                pass
            else:
                header.append(card= ('typocorr', 'True', 'Fixed incorrect year in headers, day may be off'))
                typocorr_added=True
            i.writeto(filename, overwrite = True)
        else:
            print("correct time originally")
