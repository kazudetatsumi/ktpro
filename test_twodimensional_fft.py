#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
np.set_printoptions(linewidth=120)

# https://qiita.com/kaityo256/items/64a54bb2e2c477cc6fa1
# https://numpy.org/doc/stable/reference/generated/numpy.fft.fftfreq.html
# https://stackoverflow.com/questions/14582543/extracting-frequencies-from-multidimensional-fft
n = 100

# a test image of a pyramidal intensity distribution
tx = np.linspace(-1., 1., n)
ty = np.linspace(-1., 1., n)
tz = np.zeros((n, n))
for ix, _x in enumerate(tx):
    for iy, _y in enumerate(ty):
        if _x >= 0. and _y >= 0. and _x + _y < 1.:
            tz[ix, iy] = 1. - _x - _y
        elif _x >= 0. and _y < 0. and _x - _y < 1.:
            tz[ix, iy] = 1. - _x + _y
        elif _x < 0. and _y >= 0. and _y - _x < 1.:
            tz[ix, iy] = 1. + _x - _y
        elif _x < 0. and _y < 0. and _x + _y > -1.:
            tz[ix, iy] = 1. + _x + _y
        else:
            tz[ix, iy] = 0.0

#f = np.linspace(0, n-1, n) / n
#f = np.concatenate((-f[0: int(n / 2 + 1)],
#                    f[1: int(n / 2 - 1 + 1)][::-1]))
f = np.fft.fftfreq(n, 1)

#print(f.shape)

#print(f+np.fft.fftfreq(n))

x = np.linspace(0, n, n) - n/2.

s = 4.
z = np.exp(-x**2/(2.0*s**2))/((2*np.pi)**0.5*s)



#plt.plot(x,z)
#plt.scatter(f, np.abs(np.fft.fft(z)))
#plt.scatter(f, np.exp(-(2.*np.pi*f)**2*s**2/2.))
#plt.show()


# 2D gaussian in Real space
xx, yy = np.meshgrid(x, x)
zz = np.exp(-(xx**2+yy**2)/(2.0*s**2))/(2*np.pi*s)
# 2D gaussian in Fourier space
fzz = np.fft.fft2(zz)/4.
zz2 = np.fft.ifft2(fzz)

ftz = np.fft.fft2(tz)*fzz


# 2D gaussian in Fourier space, analytically calculated
fxx, fyy = np.meshgrid(f, f)
fxx = fxx*2.*np.pi
fyy = fyy*2.*np.pi
theo_fzz = np.exp(-(fxx**2+fyy**2)*s**2/2.)

# bluring image by using 2D gaussian in Fourier space, analytically calculated
# theo_fzz = np.fft.fftshift(np.exp(-(fxx**2+fyy**2)*s**2/2.))
theo_zz2 = np.fft.ifft2(theo_fzz)
theo_ftz = np.fft.fft2(tz)*theo_fzz
theo_tz2 = np.fft.ifft2(theo_ftz)

tz2 = np.fft.ifft2(ftz)
#plt.pcolor(xx, yy, np.abs(zz))
#plt.pcolor(fxx, fyy, np.abs(fzz))
#plt.pcolor(fxx, fyy, np.abs(fzz))
#plt.imshow(np.abs(fzz))
fig = plt.figure(figsize=(8, 16))
ax = fig.add_subplot(3, 2, 1)
plt.imshow(np.clip(zz, 0, 1))
plt.title('zz')
#plt.imshow(np.clip(theo_fzz, 0, 1))
#plt.plot(theo_fzz[125,:])
#plt.plot(np.abs(fzz[125,:]))
ax = fig.add_subplot(3, 2, 2)
#plt.imshow(np.clip(np.abs(fzz), 0,1))
plt.imshow(np.clip(np.real(zz2), 0, 1))
plt.title('zz2')
ax = fig.add_subplot(3, 2, 3)
plt.imshow(np.clip(np.abs(theo_fzz), 0, 1))
plt.title('theo_fzz')
ax = fig.add_subplot(3, 2, 4)
plt.imshow(np.clip(np.fft.ifftshift(np.abs(theo_zz2)), 0, 1))
plt.title('theo_zz2')
ax = fig.add_subplot(3, 2, 5)
plt.imshow(np.clip(tz, 0, 1))
plt.title('tz')
ax = fig.add_subplot(3, 2, 6)
plt.imshow(np.clip(np.abs(theo_tz2), 0, 1))
plt.title('theo_tz2')

plt.show()


