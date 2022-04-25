#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

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

X, Y = np.meshgrid(x, y)
fftz = np.fft.fftn(z)
n = 64
#w = 0.25 / (x[1] - x[0]) 
w = 2.0 
f = np.linspace(0, n-1, n) / n
f = np.concatenate((-f[0: np.int(n / 2 + 1)],
                    f[1: np.int(n / 2 - 1 + 1)][::-1]))
print(f)
F, G = np.meshgrid(f, f)

K = np.exp(-0.5 * (w * 2 * np.pi * F) ** 2) * np.exp(-0.5 * (w * 2 * np.pi * G) ** 2)

z2 = np.real(np.fft.ifftn(z*K))

plt.pcolor(X, Y, z2)
plt.show()
