#!/usr/bin/env python
import sys
import os
home=os.path.expanduser("~")
sys.path.append("/home/kazu/ktpro")
import read_hwhm_class as rhc
import numpy as np


def samplerun():
    if len(sys.argv) >= 2:
        if sys.argv[1] == "fixbg":
            fixbg = True
        else:
            fixbg = False

    if len(sys.argv) >= 4:
        ymin = float(sys.argv[2])
        ymax = float(sys.argv[3])
    else:
        ymin = None
        ymax = None
    if len(sys.argv) >= 5:
        figfile = sys.argv[4]
    else:
        figfile = None

    #Ms = [80, 160, 320, 640, 1280, 2560]
    #pws = [0.0625, 0.125, 0.25, 0.5, 1, 2, 5]
    #pfs = ["Boxcar", "Gauss"]
    #channels = ["", "channel"]
    #bins = ["000025io", "000010io", "0000025io", "0000003io"]
    #fracs = ["", "0875", "075", "05", "0375"]
    Ms = [160]
    pws = [0.25]
    pfs = ["Boxcar"]
    channels = [""]
    #bins = ["000010io"]
    #bins = ["0000003io"]
    bins = ["000025io"]
    fracs = ["", "0875", "075", "0625",  "05", "0375", "025", "0125"]
    #fracs = ["", "0875", "075", "0625",  "05", "025", "0125"]
    #fracs = ["", "0875", "075", "0625",  "05", "0375"]
    prj = rhc.read_hwhm(Ms, pws, pfs, channels, bins, fracs,
                        prefix=home + "/desktop/210108/Tatsumi/" +
                                      "winparam_exam/direct/",
                                      numlore=2, fixbg=fixbg)
    prj.create_array()
    prj.data = prj.hwhms.squeeze()
    #prj.data = prj.data[1:,:]  # the data "" and "0875" are performed on  the same number of events.
    tcount = np.array([18970, 18970., 16098., 13380., 10621., 7794., 5318., 2826. ])
    #tcount = np.array([18970, 18970., 16098., 13380., 10621.,  5318., 2826. ])
    prj.data[:, 0] = tcount[:]/tcount[0]
    prj.data = prj.data[1:,:]
    prj.plotter(twocomp=True, isend=False, ymin=ymin, ymax=ymax)
    prj.fracs = ["0375", "05"]
    prj.bins = ["000025lastio"]
    prj.create_array()
    prj.data = prj.hwhms.squeeze()
    tcountlast = np.array([7297., 9764.])
    prj.data[:, 0] = tcountlast[:]/tcount[0]
    prj.plotter(twocomp=True, isend=True, ymin=ymin, ymax=ymax,
                c=['salmon', 'brown', 'cornflowerblue', 'midnightblue'],
                loc='center left', bbox_to_anchor=(1, 0.5), figfile=figfile)

    #print(prj.hwhms.squeeze())
    #print(prj.hwhms.squeeze().shape)


samplerun()
