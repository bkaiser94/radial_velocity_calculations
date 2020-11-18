"""
Created by Ben Kaiser 2019-05-02 (UNC-Chapel Hill)

Calls cal_params.py and then using a header as input, will return whatever desired parameter specific to the 
called needs...hopefully, as always.


"""
from __future__ import print_function

from astropy.io import fits
from astropy import coordinates as coords
from astropy import units as u

import cal_params as cp


def get_cal_params(header):
    if (('goodman' in header['INSTRUME'].lower()) or (header['INSTRUME'] =='')):
        setup_dict= cp.cal_params[header['GRATING']][header['CAM_TARG']][header['GRT_TARG']][header['INSTCONF']]
        if header['INSTRUME']=='':
            print('\n\n**********************\n"INSTRUME" header is blank "" so we are assuming it is a Goodman observation. If this is not the case, you better stop the script now because it is about to make a bunch of assumptions on that basis.\n\n**********************\n')
        else:
            pass
    elif 'gmos-n' in header['INSTRUME'].lower():
        print('GMOS-N spectrum detected')
        setup_dict=cp.gemini_params
    return setup_dict

def identify_standard(header):
    target_loc= coords.SkyCoord(header['RA'], header['DEC'], frame='icrs', unit=(u.hourangle, u.deg))
    for standard_name in cp.standard_dict:
        ra= cp.standard_dict[standard_name]['ra']
        dec= cp.standard_dict[standard_name]['dec']
        standard_loc= coords.SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg))
        distance= standard_loc.separation(target_loc)
        if distance.to(u.arcmin).value < cp.distance_threshold:
            print('Standard matched: ', standard_name)
            return standard_name
        else:
            pass
    print('No standards matched. Check that the file is in fact a standard, or maybe change cal_params.distance_threshold to be higher.')
    return
