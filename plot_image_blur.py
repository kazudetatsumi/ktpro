"""
=======================================================
Simple image blur by convolution with a Gaussian kernel
=======================================================

Blur an an image (:download:`../../../../data/elephant.png`) using a
Gaussian kernel.

Convolution is easy to perform with FFT: convolving two signals boils
down to multiplying their FFTs (and performing an inverse FFT).

"""

import numpy as np
from scipy import fftpack
import matplotlib.pyplot as plt

#####################################################################
# The original image
#####################################################################

# read image
img = plt.imread('sphx_glr_plot_image_blur_001.png')
#plt.figure()
#plt.imshow(img)
nsize = 100
x = np.linspace(-1., 1., nsize)
y = np.linspace(-1., 1., nsize)
z = np.zeros((nsize, nsize))
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

#img = z
#plt.imshow(img)
#plt.show()
#####################################################################
# Prepare an Gaussian convolution kernel
#####################################################################

# First a 1-D  Gaussian
t = np.linspace(-1., 1., nsize)
bump = np.exp(-25.*t**2)
bump /= np.trapz(bump) # normalize the integral to 1

# make a 2-D kernel out of it
kernel = bump[:, np.newaxis] * bump[np.newaxis, :]

plt.imshow(kernel)
plt.show()
#####################################################################
# Implement convolution via FFT
#####################################################################

# Padded fourier transform, with the same shape as the image
# We use :func:`scipy.signal.fftpack.fft2` to have a 2D FFT
#kernel_ft = fftpack.fft2(kernel, shape=img.shape[:-1], axes=(0, 1))
kernel_ft = np.fft.fft2(kernel)

n = nsize
f = np.linspace(0, n-1, n) / n
f = np.concatenate((-f[0: int(n / 2 + 1)],
                    f[1: int(n / 2 - 1 + 1)][::-1]))

test_kft = np.zeros((nsize,nsize))
for ix, _fx in enumerate(f):
    for iy, _fy in enumerate(f):
        test_kft[ix, iy] = np.exp(-0.01*(_fx**2 + _fy**2))

# convolve
img_ft = fftpack.fft2(img, axes=(0, 1))
#plt.imshow(np.log(np.abs(img_ft)))
#plt.show()
# the 'newaxis' is to match to color direction
#img2_ft = kernel_ft[:, :, np.newaxis] * img_ft
img2_ft = kernel_ft
#img2_ft = test_kft * img_ft
#img2 = np.real(fftpack.ifft2(img2_ft, axes=(0, 1)))
img2 = np.real(np.fft.ifft2(img2_ft))
#img2 = np.fft.ifftshift(np.real(img2))

# clip values to range
img2 = np.clip(img2, 0, 1)

# plot output
#plt.figure()
plt.imshow(img2[0:nsize,0:nsize])
plt.show()

#####################################################################
# Further exercise (only if you are familiar with this stuff):
#
# A "wrapped border" appears in the upper left and top edges of the
# image. This is because the padding is not done correctly, and does
# not take the kernel size into account (so the convolution "flows out
# of bounds of the image").  Try to remove this artifact.


#####################################################################
# A function to do it: :func:`scipy.signal.fftconvolve`
#####################################################################
#
# The above exercise was only for didactic reasons: there exists a
# function in scipy that will do this for us, and probably do a better
# job: :func:`scipy.signal.fftconvolve`

#from scipy import signal
# mode='same' is there to enforce the same output shape as input arrays
# (ie avoid border effects)
#img3 = signal.fftconvolve(img, kernel[:, :, np.newaxis], mode='same')
#plt.figure()
#plt.imshow(img3)

#####################################################################
# Note that we still have a decay to zero at the border of the image.
# Using :func:`scipy.ndimage.gaussian_filter` would get rid of this
# artifact


#plt.show()

