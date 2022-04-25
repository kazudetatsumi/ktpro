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
x = np.linspace(-1., 1., 65) 
y = np.linspace(-1., 1., 65) 
z = np.zeros((65, 65))
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

img = z
plt.imshow(img)
plt.show()
#####################################################################
# Prepare an Gaussian convolution kernel
#####################################################################

# First a 1-D  Gaussian
t = np.linspace(-10, 10, 65)
bump = np.exp(-40.1*t**2)
bump /= np.trapz(bump) # normalize the integral to 1

# make a 2-D kernel out of it
kernel = bump[:, np.newaxis] * bump[np.newaxis, :]

#####################################################################
# Implement convolution via FFT
#####################################################################

# Padded fourier transform, with the same shape as the image
# We use :func:`scipy.signal.fftpack.fft2` to have a 2D FFT
kernel_ft = fftpack.fft2(kernel, axes=(0, 1))

# convolve
img_ft = fftpack.fft2(img, axes=(0, 1))
# the 'newaxis' is to match to color direction
img2_ft = kernel_ft * img_ft
img2 = fftpack.ifft2(img2_ft, axes=(0, 1))
img2 = np.fft.ifftshift(np.real(img2))

# clip values to range
img2 = np.clip(img2, 0, 1)
print(np.argmax(img))
print(np.argmax(img2))

# plot output
plt.figure()
plt.imshow(img2)

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


plt.show()

