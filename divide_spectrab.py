#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
#from get_resampled_data_mpi_classm import Sget_qlist as sgq
from get_resampled_data_org_class import Sget_qlist as sgq


class SSdivide(sgq):
    def __init__(self, pklfile):
        super().__init__(pklfile=pklfile)
        self.load_pkl()
        print(self.spectrab.shape)
        for qidx in range(0, self.spectrab.shape[2]):
            _pklfile = self.pklfile
            self.pklfile += "." + str(qidx)
            _spectrab = self.spectrab
            self.spectrab = self.spectrab[:, :, qidx, :]
            self.save_pkl()
            self.spectrab = _spectrab
            self.pklfile = _pklfile


def samplerun():
    runNo = 3540
    #pklfile = "run" + str(runNo) + "spectraorg.pkl"
    #SSdivide(pklfile)
    pklfile = "run" + str(runNo) + "spectrab.pkl"
    SSdivide(pklfile)


samplerun()
