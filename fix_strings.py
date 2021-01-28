"""
Created by Ben Kaiser (UNC-Chapel Hill) 2020-11-27

This file holds the fix_display_string() function separately from spec_plot_tools.py because the en dash in the 
string requires some sort of declaration of encoding of the file at the beginning if it's going to work in Python 2. 
The previous solution makes all Python 2 scripts choke at import of spec_plot_tools.py which is obviously 
untenable. I will replace the spt prefix to the calls for the function in the Jupyter Notebooks and add an import call 
for this script.

"""


def fix_display_string(input_string):
    """
    INPUT: input_string, some string that needs to be examined for if it needs to be replaced with something else.
    
    
    This functions primary utility is to replace all occurrences of "GaiaJ1644-0449" with the new name "WD J1644–0449". However this can be changed very easily if we find out that name isn't acceptable either. I went ahead and had it fix the space in the other name
    """
    bad_string_list=['GaiaJ1644-0449',
                     "SDSSJ1330+6435",
                     "WDJ2356-209",
                     "Gaia J1644-0449",
                     "WD J2356-209"]
    good_string_list=["WD J1644–0449",
                      "SDSS J1330+6435",
                      "WD J2356–209",
                      "WD J1644–0449",
                      "WD J2356–209"]
    for bad_string, good_string in zip(bad_string_list, good_string_list):
        print('bad_string:', bad_string)
        print('input_string:', input_string)
        output_string=input_string.replace(bad_string, good_string)
        print('output_string:', output_string)
        if input_string!= output_string:
            return output_string
        else:
            pass
    
    #that should essentially try every bad_string inside the input_string and if it is not present, nothing will happen.
    return output_string
