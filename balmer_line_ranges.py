"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-10-09.

This file is effectively a text file, but is a python file, so I don't have to mess with importing things.

It contains the wavelength ranges for each Balmer line as recorded in Josh Fuchs' fitspec.py code. I have collected them into a list for easier accessibility however.

Presumably these values are based of those from Liebert et al 2005 (which didn't publish the actual numbers), or 


"""


import numpy as np


#The ranges that are supposed to be used to fit for each balmer line. I let the range
#go all the way out to h-alpha so I'd have all of the values there.

balmer_fit_ranges=[
    [3757.,3785.],
    [3785.,3815.],
    [3815.,3855.],
    [3859.,3925.],
    [3925.,4030.],
    [4031.,4191.],
    [4200.,4510.],
    [4680.,5040.],
    [6380.,6760.]]

balmer_norm_ranges=[
    [3500.,3782.],
    [4191.,4220.],
    [4460.,4721.],
    [5001.,6413.],
    [6713.,7000.]]
