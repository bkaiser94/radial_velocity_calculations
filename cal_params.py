"""
Created by Ben Kaiser 2019-05-02 (UNC-Chapel Hill)

This file should contain all of the trim values,  for the different gratings, camera angles , and grating angles (and possibly binnings).

It should also eventually have the different line lists for the different setups and lamps being on.

'GRATING':{
    CAM_TARG:{
        GRT_TARG:{
            'INSTCONF':{
                'linelistname': linelistname,
                'trimregions':{
                    'x':[xstart, xend]
                    'y':[ystart, yend]
                    }
                }
            }
        }
    }
    }


lamp line pixel coordinates correspond to the trimmed images that are 2x2 binned. I need to add a feature in the future that maps them to whatever binning is used.

ra and dec should be from gaia in the future, but I won't have the automatic identification implemented yet anyway.

For the spectrophotometric standards dictionary:

name_all_lower_case:{
    ra: ,
    dec: ,
    standard_file: ,
    balmer_mask: ,
    other_masks: [[firstblue, firstred],[secondblue, secondred]}
    


"""

blue_default_trim={'x':[9,2055],
                   'y':[1,199]}

red_default_trim={'x':[26,2071],
                  'y':[1,199]}


#930-ZZcetiblue
#400 M1
#400 M2
#need to add 930-ZZcetired



cal_params={
    'SYZY_930':{
        24.00000:{
            13.00000:{
                'Blue':{
                    'linelistname':'JJ_FeAr_lines.txt',
                    'offset':90,
                    'trimregions': blue_default_trim},
                'Red':{
                    'linelistname':'',
                    'offset':90,
                    'trimregions':red_default_trim}
                }
            }
        },
    'SYZY_400':{
        11.60000:{
            5.80000:{
                'Blue':{
                    'linelistname':'400M1_HgAr.txt',
                    'offset':0,
                    'trimregions':{
                        'x':[380,blue_default_trim['x'][1]],
                        'y':blue_default_trim['y']}
                    },
                'Red':{
                    'linelistname':'JJ_FeAr_lines.txt',
                    'offset':0,
                    'trimregions':{
                        'x':[380,red_default_trim['x'][1]],
                        'y':red_default_trim['y']}
                    }
                }
            },
        16.10000:{
            7.50000:{
                'Blue':{
                    'linelistname':'400M2_HgAr.txt',
                    'offset':0,
                    'trimregions':blue_default_trim},
                'Red':{
                    'linelistname':'400M2_HgAr.txt',
                    'offset':0,
                    'trimregions':red_default_trim}
                }
            }
        }
    }


#######################
####Dictionary of spectrophotometric standard-related information#####

balmer_lines_wide=[
    [3792.92, 3811.62],
    [3823.59, 3853.88],
    [3867.34,3915.21],
    [3939.52, 4029.45],
    [4046.53, 4189.13],
    [4251.3, 4470.2],
    [4661.77, 4994.76],
    [6469.59, 6703.61]
    ] #from EG274

#lines identified by Ben Kaiser looking at GD153 and EG274 in 400M2 setup
#telluric_lines=[[6803.1, 6976.18],
                #[7528.0,7760.0],
                #[7160.0,7375.0],
                #[8099.0,8405.0]
                #]

#lines identified in Moehler et al. 2014
telluric_lines=[
    [5855.,5992.],
    [6261.,6349.],
    [6438.,6600.],
    [6821.,7094.],
    [7127.,7434.],
    [7562.,7731.],
    [7801.,8613.],
    [8798.,10338.],
    [10500.,20000.]]


#old flux cal references
#standard_dict={
    #'eg274':{
        #'ra':'0',
        #'dec':'0',
        #'filename':'fhamuy/feg274.dat',
        #'sens_filename':'EG274_sensitivity_curve.txt',
        #'balmer_masks':balmer_lines_wide,
        #'other_masks': telluric_lines}
    #}
        
standard_dict={
    'eg274':{
        'ra':'0',
        'dec':'0',
        'filename':'xshooter_standards/fEg274.dat',
        'sens_filename':'EG274_sensitivity_curve.txt',
        'balmer_masks':balmer_lines_wide,
        'other_masks': telluric_lines},
    'gd153':{
        'ra':'0',
        'dec':'0',
        'filename':'xshooter_standards/fGD153.dat',
        'sens_filename':'GD153_sensitivity_curve.txt',
        'balmer_masks':balmer_lines_wide,
        'other_masks': telluric_lines}
    }




