#/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import count_intensity as ci
import numpy as np

def run():
    hdffile = '/home/kazu/desktop/200204/fine/hourbyhour/10h/out_hw_all.hdf5'
    headfile = '/home/kazu/desktop/200204/fine/hourbyhour/10h/dummy_head.txt'
    boundary = np.array([[0.475, 0.525],
                         [1.475, 1.525],
                         [-0.100, 0.100],
                         [13.0, 15.5]]).T
    print('BOUNDARY', boundary)
    counter_17714 = ci.CountInt(hdffile, headfile, boundary)
    counter_17714.getcount()

    #hdffile = '/home/kazu/desktop/200522/Ei42/dq_0.025_dw_0.5/Output4D_00_840.hdf5'
    #headfile = '/home/kazu/desktop/200522/Ei42/dq_0.025_dw_0.5/Output4D_00_840.head.txt'
    hdffile = '/home/kazu/desktop/200522/Ei42/veryfineq/14m/Output4D_00_840.hdf5'
    headfile = '/home/kazu/desktop/200522/Ei42/veryfineq/14m/Output4D_00_840.txt.head'
    boundary_tmp = np.copy(boundary)
    boundary[:, 0] = boundary_tmp[:, 1]
    boundary[:, 1] = boundary_tmp[:, 0]
    print('BOUNDARY', boundary)
    counter_Ei42 = ci.CountInt(hdffile, headfile, boundary)
    counter_Ei42.getcount()
    
    #hdffile = '/home/kazu/desktop/200522/Ei24/dq_0.025/Output4D_00_1560.hdf5'
    #headfile = '/home/kazu/desktop/200522/Ei24/dq_0.025/Output4D_00_1560.head.txt'
    hdffile = '/home/kazu/desktop/200522/Ei24/fineq/26m/Output4D_00_1560.hdf5'
    headfile = '/home/kazu/desktop/200522/Ei24/fineq/26m/Output4D_00_1560.txt.head'
    counter_Ei24 = ci.CountInt(hdffile, headfile, boundary)
    counter_Ei24.getcount()
    

run()
