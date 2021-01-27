#!/usr/bin/env python
import h5py 
import numpy as np


def save_condition(condition, outfile):
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('condition', data=condition)


def read_h5py(outfile, prop):
    f = h5py.File(outfile, 'r')
    return f[prop]

def run():
    condition = np.array(read_h5py("eliminated_data.hdf5", "condition"))
    save_condition(condition, "condition.hdf5")


run()
