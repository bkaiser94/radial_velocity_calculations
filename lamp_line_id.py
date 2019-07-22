"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-07-22


This is supposed to take a line list (I'm using one from NIST at the moment, but there's no reason it has to be
from there other than consistency of headers), and it uses that line list on a lamp image that is wavelength calibrated already using existing line lists and a different lamp image (one that doesn't have the new lines to be identified or actually might as a check).

This doesn't do any actual polynomial fitting. It should literally produce an intensity plot, and that's it.





"""

def barycentric_vel_corr(header, wavelengths):
    input_year = header['OPENDATE'] #gps-synched date
    input_hours = header['OPENTIME'] #gps-synched time
    exp_time= header['EXPTIME']*u.s
    input_times = input_year+'T'+input_hours #formatting correctly
    obs_time = Time(input_times, format = 'isot', scale = 'utc',location = cerro_pachon_location)
    obs_time= obs_time+exp_time/2.
    ra = header['RA']
    dec = header['DEC']
    radec = coords.SkyCoord(ra, dec, frame = 'icrs', unit= (u.hourangle, u.deg))
    bary_corr = radec.radial_velocity_correction(obstime= obs_time, location = cerro_pachon_location)
    bary_corr = bary_corr.to(u.km/u.s)
    lambda_rest = (wavelengths*(u.Angstrom))*const.c.to(u.km/u.s)/(-1*bary_corr+const.c.to(u.km/u.s))
    lambda_rest = lambda_rest.value
