#!/usr/bin/env python
import numpy as np
import os


def run():
    Ms = [80, 160, 320, 640, 1280, 2560]
    pws = [0.0625, 0.125, 0.25, 0.5, 1, 2, 5]
    pfs = ["Boxcar", "Gauss"]
    channels = ["", "channel"]
    bins = ["000025io", "000010io", "0000025io", "0000003io"]
    fracs = ["", "0875", "075", "05", "0375"]
    hwhfs = np.zeros((6, 7, 2, 2, 4, 5))
    for im, M in enumerate(Ms):
        for ipw, pw in enumerate(pws):
            for ipf, pf in enumerate(pfs):
                for ic, chn in enumerate(channels):
                    for ib, bn in enumerate(bins):
                        for ifr, frc in enumerate(fracs):
                            if chn == "channel":
                                if frc == "":
                                    dirname = str(M)+"_"+chn+"_"+str(pw)+"_"+bn+"_"+pf
                                else:
                                    dirname = str(M)+"_"+chn+"_"+str(pw)+"_"+bn+"_"+pf+"_"+frc
                            else:
                                if frc == "":
                                    dirname = str(M)+"_"+str(pw)+"_"+bn+"_"+pf
                                else:
                                    dirname = str(M)+"_"+str(pw)+"_"+bn+"_"+pf+"_"+frc
                            if frc == "":
                                logfile = dirname + "/std-"+bn+"-"+str(M)+"-"+str(pw)+"-"+pf+".log"
                            else:
                                logfile = dirname + "/std-"+bn+"-"+str(M)+"-"+str(pw)+"-"+pf+"-"+frc+".log"
                            if os.path.exists(logfile):
                                #print(logfile,"is found")
                                hwhfs[im, ipw, ipf, ic, ib, ifr] = gethwhf(logfile)
                            else:
                                #print(logfile,"is not found")
                                hwhfs[im, ipw, ipf, ic, ib, ifr] = 0.0
    with open('results_hwhfs.txt', 'w') as f:
        for ifr, frc in enumerate(fracs):
            for ib, bn in enumerate(bins):
                for ipf, pf in enumerate(pfs):
                    for ic, chn in enumerate(channels):
                        f.write("%s, %s, %s, %s\n" % (pf, chn, bn, frc))
                        f.write("M, 0.0625, 0.125, 0.25, 0.5, 1, 2, 5 \n")
                        for im, M in enumerate(Ms):
                            h = hwhfs[im, :, ipf, ic, ib, ifr]
                            f.write("%d, %e, %e, %e, %e, %e, %e, %e \n" % (M, h[0], h[1], h[2], h[3], h[4], h[5], h[6]))


def gethwhf(filename):
    hwhfs = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for il, line in enumerate(lines):
            if 'estimated constants' in line:
                if len(lines[il+1].split(" ")) == 4:
                    hwhfs.append(float(lines[il+1].split(" ")[1]))
    print("chk", hwhfs, filename)
    if len(hwhfs) >= 1:
        return(hwhfs[0])
    else:
        return(0.0000)


run()

