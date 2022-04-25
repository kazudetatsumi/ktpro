#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack

x = np.linspace(-1., 1., 64)
y = np.linspace(-1., 1., 64)
z = np.zeros((64, 64))
for ix, _x in enumerate(x):
    for iy, _y in enumerate(y):
        if _x >= 0. and _y >= 0. and _x + _y < 1.:
            z[ix, iy] = 1. - _x - _y
        elif _x >= 0. and _y < 0. and _x - _y < 1.:
            z[ix, iy] = 1. - _x + _y
        elif _x < 0. and _y >= 0. and _y - _x < 1.:
            z[ix, iy] = 1. + _x - _y
        elif _x < 0. and _y < 0. and _x + _y > -1.:
            z[ix, iy] = 1. + _x + _y
        else:
            z[ix, iy] = 0.0

#X, Y = np.meshgrid(x, y)
z_ft = fftpack.fft2(z, axes=(0,1))
#n = 64
##w = 0.25 / (x[1] - x[0]) 
#w = 2.0 
#f = np.linspace(0, n-1, n) / n
#f = np.concatenate((-f[0: int(n / 2 + 1)],
#                    f[1: int(n / 2 - 1 + 1)][::-1]))
#print(f)
#F, G = np.meshgrid(f, f)
#
#K = np.exp(-0.5 * (w * 2 * np.pi * F) ** 2) * np.exp(-0.5 * (w * 2 * np.pi * G) ** 2)
#
#z2 = np.real(np.fft.ifftn(z*K))

t = np.linspace(-10, 10, 64)
bump = np.exp(-10.1*t**2)
bump /= np.trapz(bump) # normalize the integral to 1
X, Y = np.meshgrid(t, t)

# make a 2-D kernel out of it
kernel = bump[:, np.newaxis] * bump[np.newaxis, :]
kernel_ft = fftpack.fft2(kernel, shape=z.shape[:], axes=(0, 1))


# the 'newaxis' is to match to color direction
z2_ft = kernel_ft[:, :, np.newaxis] * z_ft[:, :, np.newaxis]
z2 = fftpack.ifft2(z2_ft, axes=(0, 1)).real
z2 = np.clip(z2, 0, 1)


print(z2.shape)
#plt.pcolor(X, Y, np.abs(kernel_ft))
plt.imshow(z2)
plt.show()
