"""
This fixes the GPS-synched time header (well technically the year of it) because it was set to 2019 for some reason when the data was taken in 2018. This edits the original files in listSpec, so it should not be run in the RAW directory (not that anything ever really should be).

This is technically actually step1 in the step-by-step process, but it's not really a reduction step. It's a revision thing.

We're having time issues again in 2024 because they disconnected SOAR from the internet, so I'm updating this to python3 (adding 
parentheses for print statements) and changing the correct year.

2025-03-01 it happened again, so I'm going to have this just run on all the fits files in a directory since they're all screwed up. It should also 
work on all of the time headers instead of just the default one... I should also add a header that says that I fixed the time headers so there's a 
record. Ah I already did that before. Yay past Ben! Oh, I also already compiled all of the time headers instead of just the GPS time header. This 
is fun discovering features I want that I previously wanted and already added. Sad I didn't document it up here though to save myself the time. 
The globbing of all fits files in the directory didn't exist yet though, so I am adding that.
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
correct_year=2025
header_list=[
    'DATE-OBS',
    'DATE',
    'OPENDATE'
    ]

#the dict is going to indicate which index will hold the year in the split list that results from the time string. I believe they're split on hyphens.
#Well, wait, all of these headers have the same header format, at least in so far as they all have the year first... why was I casting the 2 index to an int?
#Ah for the leap year...
header_dict={
    'DATE-OBS':0,
    'DATE':0,
    'OPENDATE':0
    
    }


#speclist= np.genfromtxt(speclistname, dtype = 'str') #removing this on 2025-03-01
speclist=glob('*fits')
speclist=sorted(speclist)

for filename in speclist:
    print(filename)
    i = fits.open(filename)
    header= i[0].header
    typocorr_added=False
    for header_name in header_list:
        try:
            header_time = header[header_name]
            print(header_name, header_time)
            split_header = header_time.split('-')
            year = int(split_header[0])
            if year != correct_year:
                if ((year %4==0)^(correct_year%4==0)):
                    if year%4==0:
                        try:
                            day=int(split_header[2])
                            day=day+1
                            split_header[2]=str(day)
                        except ValueError as error:
                            #this will handle the Date-Obs header... and I just realized this is actually wrong if it's any day before February 29th... you don't add a day to the leap year day if it's... Wait why did this even run? It isn't a leap year!
                            print("ValueError:",error)
                            
                        print('year was leap year')
                        print("didn't get around to actually fixing this because it would require me to fix months too most likely")
                    elif correct_year%4==0:
                        day=int(split_header[2])
                        day=day-1
                        split_header[2]=str(day)
                        print('correct_year was leap year')
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
        except KeyError as error:
            print("KeyError:",error,filename)
            print("Probably a Zero frame because those don't open the shutter so they're missing GPS-synched headers")
