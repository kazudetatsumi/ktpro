#!/usr/bin/env python
import scipy.signal as ss
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-0.05, 0.15, 8000)
y1 = np.exp(-(x/0.01)**2)*-1.
y2 = np.exp(-np.abs(x-0.04)/0.01)
y3 = ss.convolve(y1, y2, mode='same', method='fft')*0.015



plt.plot(x, y1)
plt.plot(x, y2)
plt.plot(x, y3)
plt.show()
