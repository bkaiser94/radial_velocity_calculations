"""
Created by Ben Kaiser (UNC-Chapel Hill) 2021-05-27


2-D interpolation is getting absurd in my other files. The outputs' orders make no sense, nor do the actual interpolation functions half the damn time. I'm going to try the exact code used in the examples to see if this is nonsense.




"""
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

x = np.arange(-5.01, 5.01, 0.25)
y = np.arange(-5.01, 5.01, 0.25)
xx, yy = np.meshgrid(x, y)
z = np.sin(xx**2+yy**2)
f = interpolate.interp2d(x, y, z, kind='cubic')

xnew = np.arange(-5.01, 5.01, 1e-2)
ynew = np.arange(-5.01, 5.01, 1e-2)
znew = f(xnew, ynew)
plt.plot(x, z[0, :], 'ro-', xnew, znew[0, :], 'b-')
plt.show()

plt.contour(xnew,ynew, znew, levels=10)
plt.show()
