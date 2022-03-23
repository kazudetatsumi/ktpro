#!/usr/bin/env python
import numpy as np
import os
import sys
sys.path.append("~/ktpro")
import qens_fit_kde_hist_class as qfkhc


class read_hwhm(qfkhc.sqrun_kde_hist):
    def __init__(self, Ms, pws, pfs, channels, bins, fracs, prefix="./", numlore=1, fixbg=True):
        self.Ms = Ms
        self.pws = pws
        self.pfs = pfs
        self.channels = channels
        self.bins = bins
        self.fracs = fracs
        self.prefix = prefix
        self.numlore = numlore
        self.fixbg = fixbg

    def getfrac(self, frc):
        if frc == "":
            frac = 1.0
        else:
            frac = float("0."+frc[1:])
        return(frac)

    def create_array(self):
        self.hwhms = np.zeros((len(self.Ms), len(self.pws), len(self.pfs),
                               len(self.channels), len(self.bins),
                               len(self.fracs), 5+4*(self.numlore-1)))
        for im, M in enumerate(self.Ms):
            for ipw, pw in enumerate(self.pws):
                for ipf, pf in enumerate(self.pfs):
                    for ic, chn in enumerate(self.channels):
                        for ib, bn in enumerate(self.bins):
                            for ifr, frc in enumerate(self.fracs):
                                if chn == "channel":
                                    if frc == "":
                                        dirname = self.prefix + "/" + str(M)+"_"+chn+"_"+str(pw) +\
                                                "_"+bn+"_"+pf
                                    else:
                                        dirname = self.prefix + "/" + str(M)+"_"+chn+"_"+str(pw) +\
                                                "_"+bn+"_"+pf+"_"+frc
                                else:
                                    if frc == "":
                                        dirname = self.prefix + "/" + str(M)+"_"+str(pw)+"_"+bn +\
                                                "_"+pf
                                    else:
                                        dirname = self.prefix + "/" + str(M)+"_"+str(pw)+"_"+bn +\
                                                "_"+pf+"_"+frc
                                if frc == "":
                                    if self.numlore == 1:
                                        logfile = dirname + "/std-"+bn+"-" +\
                                         str(M)+"-"+str(pw)+"-"+pf+".log"
                                    elif self.numlore == 2:
                                        logfile = dirname + "/std-"+bn+"-" +\
                                         str(M)+"-"+str(pw)+"-"+pf+"-2lore.log"
                                else:
                                    if self.numlore == 1:
                                        logfile = dirname + "/std-"+bn+"-" +\
                                         str(M)+"-"+str(pw)+"-"+pf+"-"+frc +\
                                         ".log"
                                    elif self.numlore == 2:
                                        logfile = dirname + "/std-"+bn+"-" +\
                                         str(M)+"-"+str(pw)+"-"+pf+"-"+frc +\
                                         "-2lore.log"
                                if os.path.exists(logfile):
                                    #print(logfile,"is found")
                                    self.hwhms[im, ipw, ipf, ic, ib, ifr, 1:] =\
                                        self.gethwhm(logfile)
                                    self.hwhms[im, ipw, ipf, ic, ib, ifr, 0] =\
                                        self.getfrac(frc)
                                #else:
                                    #print(logfile,"is NOT found")


    def subouthwhm(self, prop="HWHM", lore="1st", addidx=0):
        self.f.write("%s \n" % (prop))
        for ifr, frc in enumerate(self.fracs):
            for ib, bn in enumerate(self.bins):
                for ipf, pf in enumerate(self.pfs):
                    for ic, chn in enumerate(self.channels):
                        self.f.write("%s Lorentzian %s, %s, %s, %s \n"
                                     % (lore, pf, chn, bn, frc))
                        #self.f.write("M, 0.0625, 0.125, 0.25, 0.5, 1, 2, 5 \n")
                        line = "M"
                        for _pw in self.pws:
                            line += ", %s"%(str(_pw ))
                        line += "\n"
                        self.f.write(line)
                        for im, M in enumerate(self.Ms):
                            h = self.hwhms[im, :, ipf, ic, ib, ifr, 1+addidx]
                            line = "%d"%(M)
                            for _h in h:
                                line += ", %e"%(_h)
                            line += "\n"
                            #self.f.write("%d, %e, %e, %e, %e, %e, %e, %e \n"
                            #            % (M, h[0], h[1], h[2], h[3], h[4],
                            #               h[5], h[6]))
                            self.f.write(line)



    def outhwhm(self):
        with open('results_hwhms.txt', 'w') as self.f:
            if self.numlore == 1:
                self.subouthwhm(prop="HWHM", lore="1st", addidx=0)
                self.subouthwhm(prop="STD", lore="1st", addidx=1)
                self.subouthwhm(prop="HWHMhist", lore="1st", addidx=2)
                self.subouthwhm(prop="STDhist", lore="1st", addidx=3)
            if self.numlore == 2:
                self.subouthwhm(prop="HWHM", lore="1st", addidx=0)
                self.subouthwhm(prop="HWHM", lore="2nd", addidx=2)
                self.subouthwhm(prop="STD", lore="1st", addidx=1)
                self.subouthwhm(prop="STD", lore="2nd", addidx=3)
                self.subouthwhm(prop="HWHMhist", lore="1st", addidx=4)
                self.subouthwhm(prop="HWHMhist", lore="2nd", addidx=6)
                self.subouthwhm(prop="STDhist", lore="1st", addidx=5)
                self.subouthwhm(prop="STDhist", lore="2nd", addidx=7)


        #    f.write("HWHM  \n")
        #    for ifr, frc in enumerate(self.fracs):
        #        for ib, bn in enumerate(self.bins):
        #            for ipf, pf in enumerate(self.pfs):
        #                for ic, chn in enumerate(self.channels):
        #                    f.write("1st Lorentzian %s, %s, %s, %s \n" % (pf, chn, bn, frc))
        #                    f.write("M, 0.0625, 0.125, 0.25, 0.5, 1, 2, 5 \n")
        #                    for im, M in enumerate(self.Ms):
        #                        h = self.hwhms[im, :, ipf, ic, ib, ifr, 1]
        #                        f.write("%d, %e, %e, %e, %e, %e, %e, %e \n"
        #                                % (M, h[0], h[1], h[2], h[3], h[4],
        #                                    h[5], h[6]))
        #    if self.numlore == 2:
        #        for ifr, frc in enumerate(self.fracs):
        #            for ib, bn in enumerate(self.bins):
        #                for ipf, pf in enumerate(self.pfs):
        #                    for ic, chn in enumerate(self.channels):
        #                        f.write("2nd Lorentzian %s, %s, %s, %s \n" % (pf, chn, bn, frc))
        #                        f.write("M, 0.0625, 0.125, 0.25, 0.5, 1, 2, 5 \n")
        #                        for im, M in enumerate(self.Ms):
        #                            h = self.hwhms[im, :, ipf, ic, ib, ifr, 3]
        #                            f.write("%d, %e, %e, %e, %e, %e, %e, %e \n"
        #                                    % (M, h[0], h[1], h[2], h[3], h[4],
        #                                        h[5], h[6]))
        #    f.write("STD of HWHM  \n")
        #    for ib, bn in enumerate(self.bins):
        #        for ipf, pf in enumerate(self.pfs):
        #            for ic, chn in enumerate(self.channels):
        #                f.write("1st Lorenzian %s, %s, %s, %s \n" % (pf, chn, bn, frc))
        #                f.write("M, 0.0625, 0.125, 0.25, 0.5, 1, 2, 5 \n")
        #                for im, M in enumerate(self.Ms):
        #                    h = self.hwhms[im, :, ipf, ic, ib, ifr, 3]
        #                    f.write("%d, %e, %e, %e, %e, %e, %e, %e \n"
        #                            % (M, h[0], h[1], h[2], h[3], h[4],
        #                                h[5], h[6]))
        #    if self.numlore == 2:
        #        for ib, bn in enumerate(self.bins):
        #            for ipf, pf in enumerate(self.pfs):
        #                for ic, chn in enumerate(self.channels):
        #                    f.write("2nd Lorentzian %s, %s, %s, %s \n" % (pf, chn, bn, frc))
        #                    f.write("M, 0.0625, 0.125, 0.25, 0.5, 1, 2, 5 \n")
        #                    for im, M in enumerate(self.Ms):
        #                        h = self.hwhms[im, :, ipf, ic, ib, ifr, 4]
        #                        f.write("%d, %e, %e, %e, %e, %e, %e, %e \n"
        #                                % (M, h[0], h[1], h[2], h[3], h[4],
        #                                    h[5], h[6]))



    def gethwhm(self, filename):
        #print(filename)
        hwhms = []
        stdhwhms = []
        with open(filename, 'r') as f:
            lines = f.readlines()
            for il, line in enumerate(lines):
                if 'estimated constants' in line:
                    if self.numlore == 2:
                        if lines[il+1].split()[0] == "[":
                            hwhms.append(float(lines[il+1].split()[2]))
                            hwhms.append(float(lines[il+1].split()[4]))
                        else:
                            hwhms.append(float(lines[il+1].split()[1]))
                            hwhms.append(float(lines[il+1].split()[3]))
                        stdhwhms.append(float(lines[il+4].split()[1]))
                        stdhwhms.append(float(lines[il+6].split()[3]))
                    else:
                        #hwhms.append(float(lines[il+1].split(" ")[1]))
                        hwhms.append(float(lines[il+1].split()[1]))
                        stdhwhms.append(float(lines[il+4].split(" ")[2]))
        #print("chk", hwhms, filename)
        if len(hwhms) >= 3:
            if self.numlore == 1:
                if self.fixbg:
                    return(np.array([hwhms[0], stdhwhms[0], hwhms[1], stdhwhms[1]])
                           )
                else:
                    return(np.array([hwhms[0], stdhwhms[0], hwhms[2], stdhwhms[2]])
                           )
            elif self.numlore == 2:
                if self.fixbg:
                    #print(hwhms[0], stdhwhms[0], hwhms[1], stdhwhms[1],
                    #      hwhms[2], stdhwhms[2], hwhms[3], stdhwhms[3])
                    return(np.array([hwhms[0], stdhwhms[0], hwhms[1],
                           stdhwhms[1], hwhms[2], stdhwhms[2], hwhms[3],
                           stdhwhms[3]]))
                else:
                    #print(hwhms[0], stdhwhms[0], hwhms[1], stdhwhms[1],
                    #      hwhms[4], stdhwhms[4], hwhms[5], stdhwhms[5])
                    return(np.array([hwhms[0], stdhwhms[0], hwhms[1],
                           stdhwhms[1], hwhms[4], stdhwhms[4], hwhms[5],
                           stdhwhms[5]]))
        else:
            #print("number of hwhms is smaller than 2 check logfile")
            return(np.zeros((self.numlore*4)))


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
        prj = read_hwhm(Ms, pws, pfs, channels, bins, fracs, prefix="~/desktop/210108/Tatsumi/winparam_exam")
        prj.create_array()
        prj.data = prj.hwhms.squeeze()
        prj.plotter()
        #print(prj.hwhms.squeeze())
        #print(prj.hwhms.squeeze().shape)


#samplerun()
