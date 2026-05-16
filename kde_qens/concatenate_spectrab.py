#!/usr/bin/env python
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
#from get_resampled_data_mpi_classm import Sget_qlist as sgq
from get_resampled_data_org_class import Sget_qlist as sgq


class SSconcatenate(sgq):
    def __init__(self, pklfile):
        super().__init__(pklfile=pklfile)
        self.load_pkl()
        print(self.spectrab.shape)
        _spectrab = np.array(self.spectrab)
        self.pklfile += ".bak"
        self.load_pkl()
        print(self.spectrab.shape)
        __spectrab = np.array(self.spectrab)
        self.spectrab = np.concatenate((__spectrab, _spectrab), axis=2)
        self.pklfile = self.pklfile.split(".bak")[0] + ".united"
        self.save_pkl()


def samplerun():
    runNo = 6207
    pklfile = "run" + str(runNo) + "spectraorg.pkl"
    SSconcatenate(pklfile)


samplerun()
