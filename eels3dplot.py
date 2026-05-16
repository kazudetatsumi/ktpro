#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
#from PIL import Image as im
from mpl_toolkits.mplot3d import Axes3D
from scipy import misc

dem = misc.imread('l1c1r1_combined_with_temb40p9.tif')
dem = dem[50:80,210:-390]
ny, nx = dem.shape
print dem.shape
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
xv, yv = np.meshgrid(x, y)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
dem3d=ax.plot_surface(xv,yv,dem, cmap='jet')
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

yjump=15000
eini=636
efin=eini+0.05*nx
x = np.linspace(eini, efin, nx)
fig = plt.figure(figsize=(6,6))
plt.plot(x,dem[0,:])
plt.plot(x,dem[5,:]+yjump*1)
plt.plot(x,dem[10,:]+yjump*2)
plt.plot(x,dem[15,:]+yjump*3)
plt.plot(x,dem[20,:]+yjump*4)
plt.plot(x,dem[25,:]+yjump*5)
plt.plot(x,dem[29,:]+yjump*6)
plt.yticks([])


plt.show()
#from scipy.stats import multivariate_normal


#im = im.open('l1c1r1_combined_with_temb40p9.tif')
#imarray = np.array(im)

#imsize = imarray.shape
#x1 = np.linspace(-5, 5, imsize[0])
#x2 = np.linspace(-5, 5, imsize[1])
#X1, X2 = np.meshgrid(x1, x2)
#

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#surf = ax.plot_surface(X1, X2, imarray)
#fig.show()

#m = 2 #dimension
#mean = np.zeros(m)
#sigma = np.eye(m)
#
#N = 1000
#x1 = np.linspace(-5, 5, N)
#x2 = np.linspace(-5, 5, N)
#
#X1, X2 = np.meshgrid(x1, x2)
#X = np.c_[np.ravel(X1), np.ravel(X2)]
#
#Y_plot = multivariate_normal.pdf(x=X, mean=mean, cov=sigma)
#Y_plot = Y_plot.reshape(X1.shape)
#print Y_plot.shape
#
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#surf = ax.plot_surface(X1, X2, Y_plot, cmap='bwr', linewidth=0)
#fig.colorbar(surf)
#ax.set_title("Surface Plot")
#plt.show()

