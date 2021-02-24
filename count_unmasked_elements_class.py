#!/usr/bin/env python
import numpy as np
import h5py


class Count_Unmasked_Elements:
    def __init__(self, hdffile):
        self.file = hdffile
        self.get_arrays()

    def get_arrays(self):
        f = h5py.File(self.file, 'r')
        self.data = np.array(f["data4"])
        self.condition = np.array(f["condition"])

    def get_size(self):
        return self.condition[self.condition > 0.5].size

    def get_frac(self):
        return self.condition[self.condition > 0.5].size * 1.0 \
                / self.condition.size*1.0


def samplerun():
    hdffile = "/home/kazu/desktop/200204/fine/hourbyhour/" + \
              "ortho_opt_without_mask/10h/eliminated_data.hdf5"
    pro = Count_Unmasked_Elements(hdffile)
    print("num of unmasked elements:", pro.get_size())
    print("frac of unmasked elements:", pro.get_frac())

    pro.file = "/home/kazu/desktop/200522/Ei42/veryfineq/14m/eliminated_data.hdf5"
    pro.get_arrays()
    print("num of unmasked elements:", pro.get_size())
    print("frac of unmasked elements:",
          pro.get_frac())

    pro.file = "/home/kazu/desktop/200522/Ei24/fineq/26m/eliminated_data.hdf5"
    pro.get_arrays()
    print("num of unmasked elements:", pro.get_size())
    print("frac of unmasked elements:", pro.get_frac())


#samplerun()
