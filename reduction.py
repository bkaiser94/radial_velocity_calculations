"""

"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob



zerolistname= 'listZero'
flatlistname = 'listFlat'

flatlist= np.genfromtxt(flatlistname,dtype = 'str' )
zerolist = np.genfromtxt(zerolistname, dtype ='str')

print flatlist

########

def make_image_stack(imagelist):
    """
    
    """
    images = []
    for img in imagelist:
        filename = glob(img)[0]
        i= fits.open(img)
        img_data= i[0].data
        images.append(img_data)
    return np.array(images)



###############


flat_stack = make_image_stack(flatlist)
flat_med = np.nanmedian(flat_stack, axis=0)

plt.imshow(flat_med)
plt.show()
