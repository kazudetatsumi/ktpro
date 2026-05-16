#!/usr/bin/env python
import subprocess
import numpy as np


class Get_All_TCount:

    def __init__(self, workdir, tail):
        self.workdir = workdir
        self.tail = tail

    def get_dirnames(self):
        dirs = subprocess.getoutput("/bin/ls -d " + self.workdir + "/*" +
                                    self.tail + "/").split()
        times = np.sort(np.array([int(_.split("/")[-2][:-1]) for _ in dirs]))
        self.dirnames = [str(_) + self.tail for _ in times]

    def get_tcounts(self):
        with open(self.workdir + "/result.txt_vec", mode='r') as f:
            self.counts = [
                           float(_.split()[0]) for _ in
                           f.readlines()[1:len(self.dirnames)+1]
                          ]


def samplerun():
    workdir = "/home/kazu/desktop/200204/fine/hourbyhour/" +\
              "ortho_opt_without_mask"
    tail = "h"
    proj = Get_All_TCount(workdir, tail)
    proj.get_dirnames()
    proj.get_tcounts()


#samplerun()
