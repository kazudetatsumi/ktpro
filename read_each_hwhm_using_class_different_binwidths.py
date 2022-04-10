#!/usr/bin/env python
import os
import sys
home=os.path.expanduser("~")
sys.path.append(home+"/ktpro")
import read_each_hwhm_class as rehc


def samplerun():
    Ms = [160]
    pws = [1]
    pfs = ["Boxcar"]
    channels = [""]
    bins = ["0000003io"]
    fracs = ["", "0875", "075", "0625",  "05", "0375", "025", "0125"]
    prj = rehc.read_each_hwhm(Ms, pws, pfs, channels, bins, fracs,
                              prefix=home+"/desktop/210108/" +
                              "Tatsumi/winparam_exam_6207/" +
                              "different_binwidths",
                              resultid=0)
    prj.create_array()
    prj.data = prj.hwhms.squeeze()
    prj.plotter(dataname='kde', isend=False)
    prj.bins = ["000010io"]
    prj.pws = [0.5]
    prj.resultid = 1
    prj.prefix = home+"/desktop/210108/Tatsumi/winparam_exam_6207"
    prj.create_array()
    prj.data = prj.hwhms.squeeze()
    prj.plotter(dataname='hist', isend=True)


samplerun()
