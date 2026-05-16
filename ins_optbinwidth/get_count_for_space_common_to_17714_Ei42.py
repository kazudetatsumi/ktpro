#/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import count_intensity as ci
import numpy as np

def run():
    hdffile = '/home/kazu/desktop/200204/fine/hourbyhour/10h/out_hw_all.hdf5'
    headfile = '/home/kazu/desktop/200204/fine/hourbyhour/10h/dummy_head.txt'
    boundary = np.array([[1.35, 2.125],
                         [0.75, 1.525],
                         [-0.25, 0.20],
                         [10.0, 33.0]]).T
    counter_17714 = ci.CountInt(hdffile, headfile, boundary)
    counter_17714.getcount()
    hdffile = '/home/kazu/desktop/200522/Ei42/veryfineq/14m/Output4D_00_840.hdf5'
    headfile = '/home/kazu/desktop/200522/Ei42/veryfineq/14m/Output4D_00_840.txt.head'
    boundary_tmp = np.copy(boundary)
    boundary[:, 0] = boundary_tmp[:, 1]
    boundary[:, 1] = boundary_tmp[:, 0]
    counter_Ei42 = ci.CountInt(hdffile, headfile, boundary)
    counter_Ei42.getcount()

run()
