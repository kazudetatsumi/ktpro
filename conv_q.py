#!/usr/bin/env python
# convolute a gaussian with INS 4D phantom data.
import numpy as np
import h5py 
from scipy.ndimage import gaussian_filter
from matplotlib import pyplot as plt

def read_h5py(outfile, term):
    f=h5py.File(outfile, 'r')
    return np.array(f[term])

def run():
    filleddata = "/home/kazu/desktop/200312/for_cu_new/filled/2h/eliminated_data.hdf5"
    data4 = read_h5py(filleddata, "data4")
    data4_blur = np.zeros(data4.shape)
    for idx3 in range(0, data4.shape[3]):
        data4_blur[:,:,:,idx3] = gaussian_filter(data4[:,:,:,idx3], sigma=[4.0/2.3548,2.0/2.3548,1.00/2.3548])
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(1, 2, 1)
    ax.pcolor(data4[:, 10, 16, :], cmap='jet')
    ax = fig.add_subplot(1, 2, 2)
    ax.pcolor(data4_blur[:, 10, 16, :], cmap='jet')
    plt.show()


run()
