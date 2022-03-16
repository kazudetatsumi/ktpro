#!/usr/bin/env python
# script to survey the 2 lorentzian components 2 lorentzian components fitting
# results.
# usage: read_hwhm_2comp_txt_using_class.py fixbg
import sys
sys.path.append("/home/kazu/ktpro")
import read_hwhm_class as rhc


def samplerun():
    if len(sys.argv) >= 2:
        if sys.argv[1] == "fixbg":
            fixbg = True
        else:
            fixbg = False

    #Ms = [80, 160, 320, 640, 1280, 2560]
    pws = [0.0625, 0.125, 0.25, 0.5, 1, 2, 5]
    #pfs = ["Boxcar", "Gauss"]
    #channels = ["", "channel"]
    #bins = ["000025io", "000010io", "0000025io", "0000003io"]
    #fracs = ["", "0875", "075", "05", "0375"]
    Ms = [80, 160, 320]
    Ms = [80, 160]
    #pws = [0.25, 0.5, 1]
    pws = [0.125, 0.25, 0.5, 1]
    pfs = ["Boxcar", "Gauss", "Cauchy"]
    channels = [""]
    #bins = ["000010io"]
    bins = ["0000003io"]
    #bins = ["0000025io"]
    #bins = ["000010io", "0000025io", "0000003io"]
    fracs = [""]
    #fracs = ["", "0875", "075", "0625",  "05", "0375"]
    prj = rhc.read_hwhm(Ms, pws, pfs, channels, bins, fracs,
                        prefix="./",
                        numlore=2, fixbg=fixbg)
    prj.create_array()
    prj.outhwhm()


samplerun()
