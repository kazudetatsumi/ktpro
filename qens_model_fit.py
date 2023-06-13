#!/usr/bin/env python
import numpy as np 
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_class import Sqens_balloon_resamples as qbr

class qens_model_fit(qbr):
    def __init__(self, runNo, qsize):
        preprefix = "/home/kazu/desktop/210108/Tatsumi/"
        self.orgprefix = preprefix + "from_pca03/wcorr/test/100/orgs/"
        self.stdprefix = preprefix + "from_pca03/wcorr/test/100/qs/"
        self.kdeprefix = preprefix + "winparam_exam_" + str(runNo) + \
            "/160_1_0000025io_Boxcar/n8000/"
        self.runNo = runNo
        self.qsize = qsize

    def getdata(self):
        self.outfile = self.orgprefix + "outhist" + str(self.runNo) + "m.pkl"
        self.loadfile()
        orghout = self.outall
        self.outfile = self.orgprefix + "outkde" + str(self.runNo) + "m.pkl"
        self.loadfile()
        orgkout = self.outall
        hout = np.zeros((2, self.qsize))
        for qidx in range(0, self.qsize):
            self.outfile = self.stdprefix + "runno" + str(self.runNo) +\
                           "/outhist.pkl." + str(qidx)
            self.loadfile()
            hout[0, qidx], hout[1, qidx] = self.stats()
        self.kdefile = self.kdeprefix + "kde2.log"
        masked, gamma = self.readlog()
        print(masked)

    def stats(self):
        orderidx1 = np.argsort(self.outall[:, 1])
        lbs = self.outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.16))], 1]
        ubs = self.outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.84))], 1]
        return (ubs - lbs) / 2., np.std(self.outall[:, 1])

    def readlog(self):
        masked = []
        gamma = []
        for line in open(self.kdefile):
            masked.append("-1" in line[:-1].split(']')[0])
            gamma.append(line[:-1].split('[')[2].split()[1])
        return np.array(masked), np.array(gamma)





def testrun():
    runNo = 6202
    qsize = 11
    prj = qens_model_fit(runNo, qsize)
    prj.getdata()


testrun()
