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
ref_dir= '/Users/BenKaiser/Desktop/Goodman_ref_files/'
wave_sol_dir= 'wavelength_solns/'
tell_dir='telluric_libs/'
abundance_dir='abundances/'

standard_dir=ref_dir+'standards/'
line_list_dir= ref_dir+'line_lists/'
tell_dir= ref_dir+tell_dir
abundance_dir= ref_dir+abundance_dir

airline_name='OstMart1992_airglow.txt'

blue_default_trim={'x':[9,2055],
                   'y':[1,199]}

red_default_trim={'x':[26,2071],
                  'y':[1,199]}

goodman_unbinned_pixscale= 0.15 #arcsec/pixel
#930-ZZcetiblue
#400 M1
#400 M2
#need to add 930-ZZcetired


#Goodman cal_params

cal_params={
    'SYZY_930':{
        24.00000:{
            13.00000:{
                'Blue':{
                    'linelistname':'JJ_FeAr_lines.txt',
                    'offset':90,
                    'trimregions': blue_default_trim,
                    'air_corr':True,
                    'setupname':'930blue'},
                'Red':{
                    'linelistname':'JJ_FeAr_lines.txt',
                    'offset':90,
                    'trimregions':red_default_trim,
                    'air_corr': True,
                    'setupname':'930blue'}
                }
            }
        },
    'SYZY_400':{
        11.60000:{
            5.80000:{
                'Blue':{
                    'linelistname':'400M1_HgAr.txt',
                    #'linelistname':'400m1_HgArNe_calc20190814.txt',
                    'offset':0,
                    'trimregions':{
                        'x':[380,blue_default_trim['x'][1]],
                        'y':blue_default_trim['y']},
                    'air_corr': True,
                    'setupname':'400m1'
                    },
                'Red':{
                    'linelistname':'400m1_HgArNe_calc20190814.txt',
                    'offset':0,
                    'trimregions':{
                        'x':[380,red_default_trim['x'][1]],
                        'y':red_default_trim['y']},
                    'air_corr': True,
                    'setupname':'400m1'
                    }
                }
            },
        16.10000:{
            7.50000:{
                'Blue':{
                    'linelistname':'400M2_HgAr.txt',
                    'offset':0,
                    'trimregions':blue_default_trim,
                    'air_corr':True,
                    'setupname':'400m2'},
                'Red':{
                    'linelistname':'400m2_HgArNe_calc20190814.txt',
                    'offset':0,
                    'trimregions':red_default_trim,
                    'air_corr':True,
                    'setupname':'400m2'}
                }
            }
        },
    '400_SYZY':{
        11.60000:{
            5.80000:{
                'Blue':{
                    'linelistname':'400M1_HgAr.txt',
                    'offset':0,
                    'trimregions':{
                        'x':[380,blue_default_trim['x'][1]],
                        'y':blue_default_trim['y']},
                    'air_corr': True,
                    'setupname':'400m1'
                    },
                'Red':{
                    'linelistname':'400m1_HgArNe_calc20190814.txt',
                    'offset':0,
                    'trimregions':{
                        'x':[380,red_default_trim['x'][1]],
                        'y':red_default_trim['y']},
                    'air_corr': True,
                    'setupname':'400m1'
                    }
                }
            },
        16.10000:{
            7.50000:{
                'Blue':{
                    'linelistname':'400M2_HgAr.txt',
                    'offset':0,
                    'trimregions':blue_default_trim,
                    'air_corr':True,
                    'setupname':'400m2'},
                'Red':{
                    'linelistname':'400m2_HgArNe_calc20190814.txt',
                    'offset':0,
                    'trimregions':red_default_trim,
                    'air_corr':True,
                    'setupname':'400m2'}
                }
            }
        }
    }
                
gemini_params={
    'linelistname':'no_linelist',
    'offset':0,
    'trimregions':[],
    'air_corr': True,
    'setupname':'Gemini'}


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

#Gemini-North Hamamatsu in the setup centered on 6500 angstroms
gem_gaps=[
    [5942.,5986.],
    [7023.,7080.]]

#systemic velocities taken from Simbad. See Page 37 of General Clemens VII for further info.
#all velocities in km/s
#Feige110 didn't have one listed at the top of Simbad; I assume there is one, but I haven't gotten it yet.
#GD71 = 12.0 +/- 3.4
#LTT3218 29.3+/- 2.9
#all others don't have uncertainties listed except EG274 which is +/- 0.1 km/s
#coordinates are approximately the Gaia-projected coordinates for 2019

distance_threshold= 1. #distance in arcminutes that an observation can be away from a standard to count as that standard

standard_dict={
    'eg274':{
        'ra':'16:23:33.96',
        'dec':' -39:13:46.4',
        'sys_vel':13.0,
        'filename':'xshooter_standards/fEg274.dat',
        'sens_filename':'EG274_sensitivity_curve.txt',
        'balmer_masks':balmer_lines_wide,
        'other_masks': telluric_lines},
    'gd153':{
        'ra':'12:57:02.26',
        'dec':' +22:01:48.8',
        'sys_vel': 67.0,
        'filename':'xshooter_standards/fGD153.dat',
        'sens_filename':'GD153_sensitivity_curve.txt',
        'balmer_masks':balmer_lines_wide,
        'other_masks': telluric_lines},
    'feige110':{
        'ra':'23:19:58.39',
        'dec':' -05:09:56.1',
        'sys_vel': 0.0,
        'filename':'xshooter_standards/fFeige110.dat',
        'sens_filename':'Feige110_sensitivity_curve.txt',
        'balmer_masks':balmer_lines_wide,
        'other_masks': telluric_lines},
    'ltt7987':{
        'ra':'20:10:56.35',
        'dec':' -30:13:11.4',
        'sys_vel':75.40,
        'filename':'xshooter_standards/fLTT7987.dat',
        'sens_filename':'LTT7987_sensitivity_curve.txt',
        'balmer_masks':balmer_lines_wide,
        'other_masks': telluric_lines},
    'ltt3218':{
        'ra':'08:41:32.43',
        'dec':'-32:56:32.92',
        'sys_vel':29.3,
        'filename':'xshooter_standards/fLTT3218.dat',
        'sens_filename':'LTT3218_sensitivity_curve.txt',
        'balmer_masks':balmer_lines_wide,
        'other_masks': telluric_lines},
     'gd71':{
        'ra':'08:41:32.43',
        'dec':'-32:56:32.92',
        'sys_vel':12.0,
        'filename':'xshooter_standards/fGD71.dat',
        'sens_filename':'GD71_sensitivity_curve.txt',
        'balmer_masks':balmer_lines_wide,
        'other_masks': telluric_lines},
     'eg131':{
         'ra': '19:20:34.86',
         'dec':'-07:40:02.57',
         'sys_vel':0.,
         'filename':'gemini_north_standards/eg131_flambda.dat',
         'sens_filename': 'EG131_sensitivity_curve.txt',
         'balmer_masks':balmer_lines_wide,
         'other_masks':telluric_lines+gem_gaps}
    }
     
#EG 131 does not need Balmer masks because it is a DC. Actually no it's not; it does need the Balmer masks. I don't know where we concluded that but now that I've looked it up in Simbad to find the RV, which I realized in hindsight should have actually been 0 for a DC since no absorption, it turns out EG 131 is a DBZQA. Which is pretty much every letter a dwarf can have, which would make it nearly the opposite of a DC. I don't know what was going on with me and Erik that we concluded it was a DC. ----- Actually ACTUALLY maybe this is still wrong.... the only match was actually EGGR 131 which I think is different from EG131... to be continued. End conclusion, it's a DC at low enough resolution... because it actually has absorption. It's more or less the same as my other objects though, so there's that.

#dictionary to be used for plotting atmospheric emission lines so elements are a consistent color marker
airline_color={
    'Na': 'red',
    'KI':'red',
    'Li':'red',
    '[O':'blue',
    'OH':'green',
    'O_':'magenta'
    }

#dictionary for plotting line identifications; this one should probably replace airline_color as the default dict for colors

line_color_dict={
    'Li': 'green',
    'Na': 'red',
    'K': 'blue',
    'Rb': 'magenta',
    'Ca': 'cyan'
        }

#dictionary of guesses for the box1d model used to generate the shape for background skylines
slit_airline_p0={
    'amplitude': 4.0,
    'x_0': 10.0,
    'width': 9.0
    }

slit_airline_ampbounds=(0.2,100)

cont_box_p0={
    'amplitude':1,
    'width':5000,
    'x_0': 100
    }

cont_box_bounds={
    'amplitude': (0,30),
    'width': (4000,6000),
    'x_0': (0,3000)
    }

soar_diameter= 4.1 #SOAR diameter in meters

#####

flux_cal_dict={
    'model_poly_degree':{
        '930blue':7,
        '400m1':5,
        '400m2':5,
        'Gemini':5},
    'obs_poly_degree':{
        '930blue':9,
        '400m1':7,
        '400m2':7,
        'Gemini':7},
    'sens_fit_method':{
        '930blue':'poly/poly',
        '400m1':'empirical',
        '400m2':'poly/poly',
        'Gemini':'poly/interp'}
    }

line_id_dict={
    'alkali':'alkali_lines_vac.csv',
    'ca':'Ca_1_2_lines_vac.csv'
    }


#header names to be read-out and then put into flux_calibration.py-produced sens_curves and the like
#don't include 'AIRMASS' and the method to get the MJD value because of the MJD part
in_headers = [
   'WIDTH',
   'SEE_FWHM',
   'ENVHUM',
   'ENVPRE',
   'ENVTEM',
   'ENVWIN',
   'ENVDIR'
   ]
#the names that are to be used in the sens_curve files that are written
out_headers=[
    'Ext_width',
    'See_FWHM',
    'Hum',
    'Pres',
    'Temp',
    'Wind_sp',
    'Wind_dir'
    ]


el_names = [
 "null",
 "H",
 "He",
 "Li",
 "Be",
 "B",
 "C",
 "N",
 "O",
 "F",
 "Ne",
 "Na",
 "Mg",
 "Al",
 "Si",
 "P",
 "S",
 "Cl",
 "Ar",
 "K",
 "Ca",
 "Sc",
 "Ti",
 "V",
 "Cr",
 "Mn",
 "Fe",
 "Co",
 "Ni",
 "Cu",
     ]

#el_nums=[]
#for el_num, el_name in enumerate(el_names):
    
el_nums=dict([(el[1],el[0]) for el in enumerate(el_names)])


