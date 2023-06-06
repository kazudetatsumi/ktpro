#!/usr/bin/env python
import numpy as np
import pickle
import sys
from matplotlib import pyplot as plt

if len(sys.argv) >= 2:
    prefix = sys.argv[1]
else:
    prefix = "."

def stats(outpklfile):
    with open(outpklfile, 'rb') as f:
        out = pickle.load(f)
    nsizes = np.arange(10, 3000, 1)
    ubs = np.zeros((nsizes.shape[0]))
    lbs = np.zeros((nsizes.shape[0]))
    aves = np.zeros((nsizes.shape[0]))
    stds = np.zeros((nsizes.shape[0]))
    for it, nsize in enumerate(nsizes):
        outall = out['out'][0:nsize, :]
        orderidx1 = np.argsort(outall[:, 1])
        lbs[it] = outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.16))], 1]
        ubs[it] = outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.84))], 1]
        aves[it] = np.average(outall[:, 1])
        stds[it] = np.std(outall[:, 1])
    return lbs, ubs, aves, stds
    #plt.plot(nsizes, ubs)
    #plt.plot(nsizes, lbs)
    #plt.plot(nsizes, aves)
    #plt.plot(nsizes, ubs - lbs)
    #plt.plot(nsizes, stds*2.)
    #plt.ylim(0.0012, 0.0016)
    #plt.show()


def run():
    lbsh, ubsh, avesh, stdsh = stats(prefix + '/outkde.pkl')
    print(lbsh[2989], ubsh[2989])
    #lbsk, ubsk, avesk, stdsk = stats(prefix + '/outkde.pkl')
    #plt.plot(ubsh - lbsh, label='ubsh-lbsh')
    #plt.plot(2.*stdsh, label='2stdh')
    #plt.plot(ubsk - lbsk, label='ubsk-lbsk')
    #plt.plot(2.*stdsk, label='2stdk')
    #plt.xlabel('Number of resampled spectra')
    #plt.ylabel('Uncertainties, 67% ub-lb or 2sigma (meV)')
    #plt.legend()
    #plt.show()


run()
