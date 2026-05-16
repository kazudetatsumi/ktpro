#!/usr/bin/env python
import os
import sys
home=os.path.expanduser("~")
sys.path.append(home+"/ktpro")
import read_hwhm_class as rhc
import numpy as np


def samplerun():
        #Ms = [80, 160, 320, 640, 1280, 2560]
        #pws = [0.0625, 0.125, 0.25, 0.5, 1, 2, 5]
        #pfs = ["Boxcar", "Gauss"]
        #channels = ["", "channel"]
        #bins = ["000025io", "000010io", "0000025io", "0000003io"]
        #fracs = ["", "0875", "075", "05", "0375"]
        Ms = [160]
        pws = [1]
        pfs = ["Boxcar"]
        channels = [""]
        #bins = ["000010io"]
        #bins = ["0000003io"]
        bins = ["0000025io"]
        fracs = ["", "0875", "075", "0625",  "05", "0375", "025", "0125"]
        #fracs = ["", "0875", "075", "0625",  "05", "0375"]
        prj = rhc.read_hwhm(Ms, pws, pfs, channels, bins, fracs,
                            prefix=home +
                            "/desktop/210108/Tatsumi/winparam_exam/" +
                            "different_binwidths/")
        prj.create_array()
        prj.data = prj.hwhms.squeeze()
        prj.data = prj.data[1:, :]
        tcount = np.array([18970., 16098., 13380., 10621., 7794., 5318.,
                           2826.])
        prj.data[:, 0] = tcount[:]/tcount[0]
        prj.plotter()
        #print(prj.hwhms.squeeze())
        #print(prj.hwhms.squeeze().shape)


samplerun()
