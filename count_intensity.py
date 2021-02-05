#!/usr/bin/env python
# this class helps to get the total count of INS data within a 4D space
# which is specified by the boundary in rlu and meV.
# Kazuyoshi TATSUMI 2021/02/03
import numpy as np
import h5py


class CountInt:
    def __init__(self, hdffile, headfile, boundary):
        self.hdffile = hdffile
        self.headfile = headfile
        self.boundary = boundary

    def get_initial_array_info(self):
        ln = 0
        info = []
        for line in open(self.headfile, 'r'):
            ln += 1
            if ln <= 4:
                info.append(line[:-1].split('=')[1].split(','))
            else:
                break
        self.array_info = np.array(info, dtype=np.float64)

    def calc_binnumber(self):
        self.get_initial_array_info()
        print('self.array_info[:, 0]', self.array_info[:, 0])
        binnumber_float = (self.boundary - self.array_info[:, 0]) /\
            self.array_info[:, 2]
        print('binnumber_float', binnumber_float)
        self.binnumber = np.array(np.round(binnumber_float), dtype=np.int32)
        print('binnumber', self.binnumber)

    def readhdf5(self):
        f = h5py.File(self.hdffile, 'r')
        self.data4 = np.array(f["data4"])
        self.condition = np.array(f["condition"])

    def getcount(self):
        self.calc_binnumber()
        self.readhdf5()
        xi = self.binnumber[0, 0]
        xe = self.binnumber[1, 0]+1
        yi = self.binnumber[0, 1]
        ye = self.binnumber[1, 1]+1
        zi = self.binnumber[0, 2]
        ze = self.binnumber[1, 2]+1
        ei = self.binnumber[0, 3]
        ee = self.binnumber[1, 3]+1
        print('xi,xe:', xi, xe)
        print('yi,ye:', yi, ye)
        print('zi,ze:', zi, ze)
        print('ei,ee:', ei, ee)
        print('total count:', np.sum(self.data4[xi:xe, yi:ye, zi:ze, ei:ee]))
        print('maskfrac:', np.sum(self.condition[xi:xe, yi:ye, zi:ze, ei:ee]*1.0)
              / self.condition[xi:xe, yi:ye, zi:ze, ei:ee].size*1.0)


def samplerun():
    hdffile = 'out_hw_all.hdf5'
    headfile = 'dummy_head.txt'
    boundary = np.array([[1.35, 2.125],
                         [0.75, 1.525],
                         [-0.25, 0.20],
                         [10.0, 33.0]]).T
    counter_17714 = CountInt(hdffile, headfile, boundary)
    counter_17714.getcount()


#samplerun()

