"""
take in the MJD times for the different RV's and reoutput the same file with the BMJD tdb times instead.
"""

import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
from astropy.units import cds
cds.enable()


ra = float(sys.argv[1]) #values in decimal degrees
dec = float(sys.argv[2])
cerro_pachon_location = coord.EarthLocation.from_geodetic(lat =(-30, 14, 16.41), lon = (-70, 44, 01.11), height = 2748* u.m)


target_coord = coord.SkyCoord(ra, dec, unit= (u.deg, u.deg), frame= 'icrs')
input_file= 'radial_velocities_combinedb.txt'
output_file= 'radial_velocities_combinedb_bmjd.txt'
all_array = np.genfromtxt(input_file).T

def to_barycenter(input_times):
    bary_corr =input_times.tdb.light_travel_time(target_coord)
    return (input_times.tdb+ bary_corr.tdb).mjd




mjd_array = all_array[0]
og_mjd = np.copy(mjd_array)
obs_times = Time(mjd_array, format = 'mjd', scale = 'utc', location = cerro_pachon_location)

obs_times_bary = to_barycenter(obs_times)
mjd_array = obs_times_bary

header = "BMJD_TDB,H_delta,H_gamma,H_beta"
np.savetxt(output_file, all_array.T, delimiter = ',', header = header)
