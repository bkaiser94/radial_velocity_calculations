"""
Created by Ben Kaiser 2019-05-02 (UNC-Chapel Hill)

Calls cal_params.py and then using a header as input, will return whatever desired parameter specific to the 
called needs...hopefully, as always.


"""
from astropy.io import fits


import cal_params as cp


def get_cal_params(header):
    setup_dict= cp.cal_params[header['GRATING']][header['CAM_TARG']][header['GRT_TARG']][header['INSTCONF']]
    return setup_dict
