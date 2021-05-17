"""
Created by Ben Kaiser (UNC-Chapel Hill) 2021-05-17

Take the white dwarf cooling model files from Bedard et al. 2020 (https://www.astro.umontreal.ca/
~bergeron/CoolingModels/), and compile them into single files that include all of the relevant information to 
be able to interpolate between them.

The files are structured in a strange way by default in which there are essentially lines with missing columns, 
namely the Mod column. Not to mention that it includes a #, so that will be shooting us in the foot when we 
try to read that in most likely. I'm not sure how to handle that one... I hope it doesn't comment out the whole 
header portion but I bet it will now...



"""
