"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-10-09.

This file is effectively a text file, but is a python file, so I don't have to mess with importing things.

It contains the wavelength ranges for each Balmer line as recorded in Josh Fuchs' fitspec.py code. I have collected them into a list for easier accessibility however.

Presumably these values are based of those from Liebert et al 2005 (which didn't publish the actual numbers), or 


"""


#import numpy as np


#The ranges that are supposed to be used to fit for each balmer line. I let the range
#go all the way out to h-alpha so I'd have all of the values there.

def make_inside_out(input_list, min_val, max_val):
    """
    Takes a list of lists (usually something like a mask_list from the other scripts) and then makes the selected regions become the outer boundaries.
    
    I.E. make_inside_out([[10,12],[14,18]], 8, 30)
    returns
    [[8,10],[12,14],[18,30]]
    which can then be used as a mask_list to do other stuff
    """
    new_list=[[]]
    for index in range(0,len(input_list)+1):
        if index==0:
            new_list.append([min_val,input_list[index][0]])
        elif index==len(input_list):
            new_list.append([input_list[index-1][1],max_val])
        else:
            new_list.append([input_list[index-1][1],input_list[index][0]])
        #print new_list
    new_list= new_list[1:]
    #print new_list
    return new_list
        
        
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

io_balmer_fit_ranges=make_inside_out(balmer_fit_ranges, 3600, 5200)
#print "io_balmer_fit_ranges:", io_balmer_fit_ranges

#balmer_norm_ranges=[
    #[3500.,3782.],
    #[4191.,4220.],
    #[4460.,4721.],
    #[5001.,6413.],
    #[6713.,7000.]]

balmer_norm_ranges=[
    [3500.,3782.],
    [4025.,4035.],
    [4191.,4220.],
    [4460.,4721.],
    [5001.,6413.],
    [6713.,7000.]] #new version as of 20190301. Added a continuum region in he 4000 area
io_balmer_norm_ranges=make_inside_out(balmer_norm_ranges, 3600, 5200)
#balmer_norm_masks=[
    #[3782., 4191.],
    #[4220., 4460.],
    #[4721., 5001.]]


balmer_norm_masks=[
    [3782., 4191.],
    [4220., 4460.],
    [4721., 5001.]]

continuum_list = [[3812,3819],
                  [3858,3864],
                  [3924,3932],
                  [4012,4022],
                  [4034,4040],
                  [4172, 4260],
                  [4422,4460],
                  [4500,4780],
                  [4950,5010],
                  [5035, 5200]]#Best one there is. shortened red side

io_continuum_list= make_inside_out(continuum_list, 3600, 5200)
