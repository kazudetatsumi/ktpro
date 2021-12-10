#!/usr/bin/env python
import numpy as np
import os
import sys
sys.path.append("/home/kazu/ktpro")
import qens_fit_kde_hist_class as qfkhc


class read_each_hwhm(qfkhc.sqrun_kde_hist):
    def __init__(self, Ms, pws, pfs, channels, bins, fracs, prefix="./",
                 resultid=0):
        self.Ms = Ms
        self.pws = pws
        self.pfs = pfs
        self.channels = channels
        self.bins = bins
        self.fracs = fracs
        self.prefix = prefix
        self.resultid = resultid

    def getfrac(self, frc):
        if frc == "":
            frac = 1.0
        else:
            frac = float("0."+frc[1:])
        return(frac)

    def create_array(self):
        self.hwhms = np.zeros((len(self.Ms), len(self.pws), len(self.pfs),
                               len(self.channels), len(self.bins),
                               len(self.fracs), 3))
        for im, M in enumerate(self.Ms):
            for ipw, pw in enumerate(self.pws):
                for ipf, pf in enumerate(self.pfs):
                    for ic, chn in enumerate(self.channels):
                        for ib, bn in enumerate(self.bins):
                            for ifr, frc in enumerate(self.fracs):
                                if chn == "channel":
                                    if frc == "":
                                        dirname = self.prefix + "/" + str(M) +\
                                                  "_"+chn+"_"+str(pw)+"_"+bn +\
                                                  "_"+pf
                                    else:
                                        dirname = self.prefix + "/" + str(M) +\
                                                  "_"+chn+"_"+str(pw)+"_"+bn +\
                                                  "_"+pf+"_"+frc
                                else:
                                    if frc == "":
                                        dirname = self.prefix + "/" + str(M) +\
                                                  "_"+str(pw)+"_"+bn+"_"+pf
                                    else:
                                        dirname = self.prefix + "/" + str(M) +\
                                                  "_"+str(pw)+"_"+bn+"_"+pf +\
                                                  "_"+frc
                                if frc == "":
                                    logfile = dirname + "/std-"+bn+"-" +\
                                            str(M)+"-"+str(pw)+"-"+pf+".log"
                                else:
                                    logfile = dirname + "/std-"+bn+"-" +\
                                            str(M)+"-"+str(pw)+"-"+pf +\
                                            "-"+frc+".log"
                                if os.path.exists(logfile):
                                    #print(logfile,"is found")
                                    self.hwhms[im, ipw, ipf, ic, ib, ifr, 1:] =\
                                        self.getonehwhm(logfile)
                                    self.hwhms[im, ipw, ipf, ic, ib, ifr, 0] =\
                                        self.getfrac(frc)

    def getonehwhm(self, filename):
        hwhms = []
        stdhwhms = []
        with open(filename, 'r') as f:
            lines = f.readlines()
            for il, line in enumerate(lines):
                if 'estimated constants' in line:
                    hwhms.append(float(lines[il+1].split()[1]))
                    stdhwhms.append(float(lines[il+4].split()[1]))
        #print("chk", hwhms, filename)
        # we assume the log file of the fitting contains the fitting paramters
        # in the correct order, i.e., firstly, those obtained using the kde,
        # secondly, those obtained using the histograms with the background was
        # determined by the
        # results of the fitting with the kde, and finally those obtained uisng
        # the histograms.
        return(np.array([hwhms[self.resultid], stdhwhms[self.resultid]]))


def samplerun():
        #Ms = [80, 160, 320, 640, 1280, 2560]
        #pws = [0.0625, 0.125, 0.25, 0.5, 1, 2, 5]
        #pfs = ["Boxcar", "Gauss"]
        #channels = ["", "channel"]
        #bins = ["000025io", "000010io", "0000025io", "0000003io"]
        #fracs = ["", "0875", "075", "05", "0375"]
        Ms = [160]
        pws = [0.5]
        pfs = ["Boxcar"]
        channels = [""]
        bins = ["0000025io"]
        fracs = ["", "0875", "075", "0625",  "05", "0375"]
        prj = read_hwhm(Ms, pws, pfs, channels, bins, fracs, prefix="/home/kazu/desktop/210108/Tatsumi/winparam_exam")
        prj.create_array()
        prj.data = prj.hwhms.squeeze()
        prj.plotter()
        #print(prj.hwhms.squeeze())
        #print(prj.hwhms.squeeze().shape)


#samplerun()
