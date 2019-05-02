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
                    'linelistname':'',
                    'trimregions': blue_default_trim},
                'Red':{
                    'linelistname':'',
                    'trimregions':red_default_trim}
                }
            }
        },
    'SYZY_400':{
        11.60000:{
            5.80000:{
                'Blue':{
                    'linelistname':'',
                    'trimregions':{
                        'x':[380,blue_default_trim['x'][1]],
                        'y':blue_default_trim['y']}
                    },
                'Red':{
                    'linelistname':'',
                    'trimregions':{
                        'x':[380,red_default_trim['x'][1]],
                        'y':red_default_trim['y']}
                    }
                }
            },
        16.10000:{
            7.50000:{
                'Blue':{
                    'linelistname':'',
                    'trimregions':blue_default_trim},
                'Red':{
                    'linelistname':'',
                    'trimregions':red_default_trim}
                }
            }
        }
    }
