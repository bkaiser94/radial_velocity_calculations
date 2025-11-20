"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-11-20


I need this to open all of the compressed FITS files Jon gave me and then resave them as uncompressed FITS files because that's probably what's breaking my code.

Maybe I'll have this one trim the images down to size at the same time while I'm at it.



Copilot assisted.

"""

from astropy.io import fits
from glob import glob


input_list=glob("*.fits.fz")
trace_coordinate=500 #trace coordinate in y values
desired_height=200


def decompress_file(input_file):
    # open compressed FITS
    with fits.open(input_file) as hdu:
        hdu.info()
        # primary header
        header = hdu[1].header
        # image or table data (replace 1 with the HDU index you need)
        #hdu.compressed._tiled_compression.decompress_image_data_section
        data = hdu[1].data
        #print(header)
        print(type(data), data.shape)
        output_name_parts=input_file.split('.')
        output_name=".".join(output_name_parts[:-1]) #Drop the .fz from the name
        if data.shape[1]==2071:
            #Apparently there are quite a few biases kicking around in here that shouldn't be. so this will limit to things binned 2x2 hopefully. Especially a problem with biases. They're just all in there together it seems
            sub_array=data[trace_coordinate-desired_height//2:trace_coordinate+desired_height//2,:]
            print('sub_array.shape',sub_array.shape)
            #print(header)
            #for thing in header:
                #print(thing)
            #print('\n\nNext Header\n\n')
            header.append(card=('downtrim',True,'trimmed CCD segment in postprocess'))
            #print(header)
            #for thing in header:
                #print(thing)
            hdu_out= fits.PrimaryHDU(sub_array, header=header)
            hdu_out.writeto(output_name, overwrite=True)
            print('image saved.',output_name)
        else:
            pass

for input_file in input_list:
    decompress_file(input_file)
