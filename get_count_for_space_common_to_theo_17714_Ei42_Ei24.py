#/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import count_intensity as ci
import numpy as np

def run():
    #hdffile = '/home/kazu/desktop/200204/fine/hourbyhour/10h/out_hw_all.hdf5'
    hdffile = '/home/kazu/desktop/200312/for_cu_new/old_filled/filled_again2_ddscs/10h/eliminated_data.hdf5'
    hdffile = '/home/kazu/desktop/200312/for_cu_new/old_filled/filled_ddscs_kfki_debye/5h/eliminated_data.hdf5'
    hdffile = '/home/kazu/desktop/200312/for_cu_new/old_filled/filled_ddscs_kfki_debye/10h/eliminated_data.hdf5'
    headfile = '/home/kazu/desktop/200204/fine/hourbyhour/10h/dummy_head.txt'
    boundary = np.array([[0.475, 0.525],
                         [1.475, 1.525],
                         [-0.100, 0.100],
                         [13.0, 15.5]]).T
    counter_17714 = ci.CountInt(hdffile, headfile, boundary)
    counter_17714.getcount()

    #hdffile = '/home/kazu/desktop/200522/Ei42/dq_0.025_dw_0.5/Output4D_00_840.hdf5'
    hdffile = '/home/kazu/desktop/200701/all/filled_again_ddscs/14m/eliminated_data.hdf5'
    hdffile = '/home/kazu/desktop/200701/all/filled_again_ddscs2/4m/eliminated_data.hdf5'
    hdffile = '/home/kazu/desktop/200701/all/filled_again_ddscs2/14m/eliminated_data.hdf5'

    #headfile = '/home/kazu/desktop/200522/Ei42/dq_0.025_dw_0.5/Output4D_00_840.head.txt'
    headfile = '/home/kazu/desktop/200522/Ei42/veryfineq/14m/Output4D_00_840.txt.head'
    boundary_tmp = np.copy(boundary)
    boundary[:, 0] = boundary_tmp[:, 1]
    boundary[:, 1] = boundary_tmp[:, 0]
    counter_Ei42 = ci.CountInt(hdffile, headfile, boundary)
    counter_Ei42.getcount()
    
    #hdffile = '/home/kazu/desktop/200522/Ei24/dq_0.025/Output4D_00_1560.hdf5'
    #hdffile = '/home/kazu/desktop/200903/all/filled_ddscs_kfki_debye/12m/eliminated_data.hdf5'
    hdffile = '/home/kazu/desktop/200903/all/filled_ddscs_kfki_debye/26m/eliminated_data.hdf5'
    hdffile = '/home/kazu/desktop/200903/all/filled_gamma_corrected/54m/eliminated_data.hdf5'
    headfile = '/home/kazu/desktop/200522/Ei24/fineq/26m/Output4D_00_1560.txt.head'
    counter_Ei24 = ci.CountInt(hdffile, headfile, boundary)
    counter_Ei24.getcount()
    

run()
