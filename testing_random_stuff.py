from matplotlib import pyplot as plt
import numpy as np


xvals=[1013.,1018.,1022., 1022.5,1018.5,1013.5,1023.5,1028.,1032.,1032.5,1028.5,1024.]

yvals=[136.5,132., 129., 77.5, 75., 71., 146., 143., 140., 71.5, 69., 66.5]

radius=34.4
xcirc=np.linspace(998.-radius,998.+radius,1000)
yupcirc=np.sqrt(radius**2-(xcirc-998.)**2)
ydowncirc=yupcirc*-1+104.
yupcirc=yupcirc+104.


radius2=np.sqrt((1028.5-998.)**2+(69.-104.)**2)
xcirc2=np.linspace(998.-radius2,998.+radius2,1000)
yupcirc2=np.sqrt(radius2**2-(xcirc2-998.)**2)
ydowncirc2=yupcirc2*-1+104.
yupcirc2=yupcirc2+104.

plt.scatter(xvals,yvals,s=4)

plt.scatter(998.,104.,color='r',marker='o',s=14)
plt.plot(xcirc,yupcirc,label='upcirc r='+str(radius)[:4])
plt.plot(xcirc,ydowncirc,label='downcirc r='+str(radius)[:4])
plt.plot(xcirc2,yupcirc2,label='upcirc r='+str(radius2)[:4])
plt.plot(xcirc2,ydowncirc2,label='downcirc r='+str(radius2)[:4])
plt.legend()

plt.xlim(900,1100)
plt.ylim(0,200)

plt.show()
