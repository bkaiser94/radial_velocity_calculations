"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-04-01

This should take the spectroscopic FITS files and rename them using the common name in 
the MORDOR survey CSV files.

When we undertook the MORDOR Survey starting back in 2019, we were pretty fast and loose with the naming. And I was apparently dyslexic. As a consequence the majority of the files have names like GaiaJ1644m0449_400m1.fits (in the case of WD J1644-0449). However, these aren't even self-consistent, for example capitalization etc.

Therefore, it makes sense to make this file, which will take the "name" of each object in the MORDOR survey file and then rename each file to have that name base in there and to have conventional capitalization. It will also change the "400m1" and "400m2" columns of the table to match the new file names for each object.

I should probably output that new table to yet another new CSV file instead of overwriting the existing one because being careful.



"""

from __future__ import print_function
from glob import glob
import sys
from astropy.table import Table, Column
import os
import time


sys.path.append('../')
sys.path.append('/Users/BenKaiser/Desktop/radial_velocity_calculations/')

dest_path='/Users/BenKaiser/Desktop/MORDOR_Survey_paper/MORDOR_Survey_forPaper/Goodman_spectra_copy/'
#based on the internet, all I need to do this is basically os.rename(original_name, new_name)

#need the CSV file with the filenames and such inside it

#apparently I need to first actually go down into the subfiles for each spectral type and bring those back up into the main CSV to get the filenames.

input_file='full_MORDOR_survey_1681141339_fullphotometry_fitparams.csv'

input_table=Table.read(input_file)

