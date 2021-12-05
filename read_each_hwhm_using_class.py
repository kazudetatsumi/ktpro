#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import read_each_hwhm_class as rehc


def samplerun():
        Ms = [160]
        pws = [0.5]
        pfs = ["Boxcar"]
        channels = [""]
        bins = ["0000025io"]
        fracs = ["", "0875", "075", "0625",  "05", "0375", "025", "0125"]
        prj = rehc.read_each_hwhm(Ms, pws, pfs, channels, bins, fracs,
                                  prefix="/home/kazu/desktop/210108/" +
                                         "Tatsumi/winparam_exam_6207",
                                  resultid=0)
        prj.create_array()
        prj.data = prj.hwhms.squeeze()
        prj.plotter(dataname='kde')
        prj.bins = ["000010io"]
        prj.resultid = 2
        prj.create_array()
        prj.data = prj.hwhms.squeeze()
        prj.plotter(dataname='hist', isend=True)


samplerun()
