#!/usr/bin/env python
import numpy as np
import pickle
import sys
from matplotlib import pyplot as plt

if len(sys.argv) >= 2:
    prefix = sys.argv[1]
else:
    prefix = "."

plt.rcParams['font.size'] = 14
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


def run(outfile):
    lbsh, ubsh, avesh, stdsh = stats(prefix + "/" + outfile)
    print((ubsh[2989]-lbsh[2989])/2., stdsh[2989])
    lbsk, ubsk, avesk, stdsk = stats(prefix + '/outkde.pkl')
    plt.plot((ubsh - lbsh)*1000., label='$CI_h$')
    plt.plot(2.*stdsh*1000., label='2$\sigma_h$')
    plt.plot((ubsk - lbsk)*1000., label='$CI_{KDE}$')
    plt.plot(2.*stdsk*1000., label='2$\sigma_{KDE}$')
    plt.xlabel('Number of resampled spectra')
    plt.ylabel('67% CI or 2$\sigma$ of $\Gamma$ ($\mu$eV)')
    #plt.ylim(1, 2.)
    plt.xlim(0, 3000.)
    #plt.yticks([0., 0.5, 1.0, 1.5, 2.0])
    plt.xticks([0, 1000, 2000, 3000])
    plt.tick_params(direction='in', right=True, top=True, labelbottom=True)
    plt.legend()
    plt.show()


def mstats():
    print("68%_ci/2 std")
    for qidx in range(0, 11):
        outfile = "outhist.pkl." + str(qidx)
        run(outfile)

outfile = "outhist.pkl"
run(outfile)

#mstats()
