#!/usr/bin/env python
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt

sig = np.repeat([0., 1., 0.], 500000)
win = signal.windows.hann(250000)
filtered = signal.convolve(sig, win, mode='same', method='direct') / sum(win)
#filtered = signal.convolve(sig, win, mode='same', method='fft') / sum(win)


#fig, (ax_orig, ax_win, ax_filt) = plt.subplots(3, 1, sharex=True)
#ax_orig.plot(sig)
#ax_orig.set_title('Original pulse')
#ax_orig.margins(0, 0.1)
#ax_win.plot(win)
#ax_win.set_title('Filter impulse response')
#ax_win.margins(0, 0.1)
#ax_filt.plot(filtered)
#ax_filt.set_title('Filtered signal')
#ax_filt.margins(0, 0.1)
#fig.tight_layout()
#plt.show()
