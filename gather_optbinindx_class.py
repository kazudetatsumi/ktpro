#!/usr/bin/env python
import numpy as np
import subprocess
import h5py
import os
import sys
import pickle
from matplotlib import pyplot as plt
import matplotlib.image as mpimg

class gather_optbinidx:
    def __init__(self, num_try, workdir, timestr, reflog, savefile, dq, eu):
        self.num_try = num_try
        self.workdir = workdir
        self.timestr = timestr
        self.reflog = reflog
        self.savefile = savefile
        self.dq = dq
        self.eu = eu
        print("num_try:", num_try)

    def getbinidx(self, logfile):
        #print(logfile)
        line = subprocess.getoutput("tail --line=1 " + logfile)
        if "(" in line:
            values = line.split("(")
            values2 = values[1].split(")")
            values3 = np.array(values2[0].split(","), dtype='int32')
        elif "minloc" in line:
            values = line.split()
            values3 = np.array(values[2:6][::-1], dtype='int32')
        else:
            line2 = subprocess.getoutput("grep \"with the\" " + logfile)
            values = line2.split("(")
            values2 = values[1].split(")")
            values3 = np.array(values2[0].split(","), dtype='int32')
        return values3

    def getscalarfromhdf(self, hdffile, prop):
        f = h5py.File(hdffile, 'r')
        return f[prop].value

    def init_diffbin(self):
        #print(self.reflog)
        refbinidx = self.getbinidx(self.reflog)
        print(refbinidx)

        binindx = np.zeros((self.num_try, 4), dtype='int32')
        maskfrac = np.zeros((self.num_try))
        radius = np.zeros((self.num_try))
        for tryidx in range(0, self.num_try):
            path = self.workdir + "try" + str(tryidx) + "_" + self.timestr
            logfile = path + "std-" + str(tryidx) + ".log"
            binindx[tryidx, :] = self.getbinidx(logfile)
            hdffile = path + "condition.hdf5"
            maskfrac[tryidx] = self.getscalarfromhdf(hdffile, 'frac_add_mask')
            radius[tryidx] = self.getscalarfromhdf(hdffile, 'radius')
        difbinidx = np.max(np.abs(binindx - refbinidx), axis=1)

        dataset = {}
        dataset['maskfrac'] = maskfrac
        dataset['difbinidx'] = difbinidx
        dataset['radius'] = radius

        with open(self.savefile, 'wb') as f:
            pickle.dump(dataset, f, -1)
        print("Done!")

    def create_fig(self, dataname):
        fig = plt.figure(figsize=(16, 8))
        fig.suptitle("optimal bin-widths maximum difference caused" +
                     " by an additional cyrindrical mask \n for " + dataname)

    def get_max_maskfrac(self):
        if not os.path.exists(self.savefile):
            self.init_diffbin()
        with open(self.savefile, 'rb') as f:
            dataset = pickle.load(f)
        self.difbinidx = dataset['difbinidx']
        self.maskfrac = dataset['maskfrac']
        self.radius = dataset['radius']
        if np.max(self.difbinidx) >= 2:
            self.pngdatano = np.where((self.difbinidx >= 2) &
                                      (self.maskfrac ==
                                       np.min(
                                           self.maskfrac[self.difbinidx >= 2]))
                                      )[0][0]
            self.max_maskfrac = self.maskfrac[self.pngdatano]

    def process(self):
        self.get_max_maskfrac()
        ax = plt.subplot2grid((2, 2), (0, 1), rowspan=2)
        if np.max(self.difbinidx) >= 2:
            print("pngdatano:", self.pngdatano)
            print("that maskfrac:", self.maskfrac[self.pngdatano])
            pngfile = self.workdir + "/data4_" + self.timestr[:-2] + "_" +\
                str(self.pngdatano)+".png"
            print(pngfile)
            if os.path.exists(pngfile):
                plt.imshow(mpimg.imread(pngfile))
            else:
                pngfile = self.workdir + "/../data4_" + self.timestr[:-2] +\
                          "_" + str(self.pngdatano)+".png"
                if os.path.exists(pngfile):
                    plt.imshow(mpimg.imread(pngfile))
                else:
                    sys.exit("No adequate png file exists!")

        ax = plt.subplot2grid((2, 2), (0, 0))
        ax.scatter(self.radius*self.dq*2.0, self.difbinidx, marker='x', c='k',
                   clip_on=True)
        if np.max(self.difbinidx) >= 2:
            ax.scatter(self.radius[self.pngdatano]*self.dq*2.0,
                       self.difbinidx[self.pngdatano],
                       marker='x', facecolors='gray', edgecolor='gray',
                       clip_on=True)

        ax.set_yticks(range(np.min(self.difbinidx), np.max(self.difbinidx)+1))
        ax.set_xlabel('mask size (rlu or ' + str(self.eu) + ' meV)')
        ax.set_ylabel('maximum difference in bin-width (step)')
        ax = plt.subplot2grid((2, 2), (1, 0))
        ax.scatter(self.maskfrac*100.0, self.difbinidx, marker='x', c='k',
                   clip_on=True)
        if np.max(self.difbinidx) >= 2:
            ax.scatter(self.maskfrac[self.pngdatano]*100.0,
                       self.difbinidx[self.pngdatano],
                       marker='x',  facecolors='gray', edgecolor='gray',
                       clip_on=True)
        ax.set_yticks(range(np.min(self.difbinidx), np.max(self.difbinidx)+1))
        ax.set_xlabel('mask volume fraction (%)')
        ax.set_ylabel('maximum difference in bin-width (step)')
        plt.show()


def samplerun():
    num_try = 800
    timestr = "14m/"
    workdir = "/home/kazu/desktop/200522/Ei42/veryfineq/add_random_mask/" +\
              timestr + "condparam09/"
    reflog = "/home/kazu/desktop/200522/Ei42/veryfineq/condparam_09/" +\
             timestr + "calc.log"
    savefile = workdir + "difbin.pkl"
    dq = 0.025
    dE = 0.2
    eu = int(dE/dq)
    dataname = "Ei42"
    project = gather_optbinidx(num_try, workdir, timestr, reflog, savefile, dq,
                               eu)
    project.create_fig(dataname)
    project.process()


#samplerun()
