#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import gather_optbinindx_class as goc
import subprocess
import os
import numpy as np


class maxfracall(goc.gather_optbinidx):

    def __init__(self, tail, cond="/"):
        self.tail = tail
        self.cond = cond

    def get_all_maxfrac(self):
        dirs = subprocess.getoutput("/bin/ls -d *" + self.tail +"/").split()
        times = np.sort(np.array([int(_[:-2]) for _ in dirs]))
        maxmaskfracs = []
        dirlist = []
        for time in times:
            dirname = str(time) + self.tail
            workdir = "./" + dirname + "/" + self.cond + "/"
            if os.path.isfile(workdir + "difbin.pkl"):
                self.savefile = workdir + "difbin.pkl"
                super(maxfracall, self).get_max_maskfrac()
                maxmaskfracs.append(self.max_maskfrac)
                dirlist.append(dirname)

        if self.cond == "/":
            txtfile = "maxfrac.txt"
        else:
            txtfile = "maxfrac_" + self.cond + ".txt"
        with open(txtfile, mode='w') as f:
            for dd, mm in zip(dirlist, maxmaskfracs):
                f.write("%s %f \n" % (dd, mm))


def samplerun():
    tail = "m"
    cond = "/condparam09/"
    proj = maxfracall(tail, cond)
    proj.get_all_maxfrac()


#samplerun()
