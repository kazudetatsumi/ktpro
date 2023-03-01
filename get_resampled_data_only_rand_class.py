#!/usr/bin/env python
# This script generates randomstates and saves them as a pickle file
# in order to draw bootstrap samples of QENS data by setting the random
# states in the file.
# It requres inputs of runNo, the number of the randomstates, and the TimeParam
# for the utsusemi data reduction, because the generated randomstates
# depends on the total neutron counts in the QENS data.
# Kazuyoshi TATSUMI 2023/02/28

import numpy as np
import pickle
import sys
sys.path.append("/home/kazu/ktpro")
from get_resampled_data_mpi_class import Sget_qlist as sq


class Get_rand(sq):
    def __init__(self, save_file=None, pklfile=None):
        self.save_file = save_file
        self.pklfile = pklfile

    def get_randomstates(self, runNo, nbs, seed=314):
        intensity1d = self.intensity.flatten().astype(int)
        nonzeroidx = np.nonzero(intensity1d)[0]
        x = np.array([idx for idx in nonzeroidx for num_repeat in
                     range(intensity1d[idx])], dtype=int)
        N = x.shape[0]
        np.random.seed(seed)
        randomstates = []
        for inb in range(nbs):
            randomstates.append(np.random.get_state())
            Nb = np.random.poisson(lam=N)
            idx = np.random.randint(0, N, Nb)
        with open('randomstates.pkl' + '.' + str(runNo) + '.' + str(nbs), 'wb') as f:
            pickle.dump(randomstates, f, -1)


def run(runNo, nbs, TimeParam):
    prj = Get_rand()
    prj.get_org_data("0.000025", runNo, TimeParam=TimeParam)
    prj.get_org_intensity_array()
    prj.get_randomstates(runNo, nbs, seed=314)


#run(6204, 30,  "-1.0/-1.0")
#run(6202, 30000, True,  "-1.0/-1.0")

