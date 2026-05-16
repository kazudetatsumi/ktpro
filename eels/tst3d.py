from scipy import misc
dem=misc.imread('dem.tif')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

print dem.shape

import numpy as np
ny, nx = dem.shape
print dem

x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
xv, yv = np.meshgrid(x, y)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
dem3d=ax.plot_surface(xv,yv,dem)
plt.show()


