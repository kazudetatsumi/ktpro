#!/usr/bin/env python
import numpy as np
import h5py
import yaml


def get_data(fn):
    f = h5py.File(fn)
    q = f["qpoint"]      # (nqps, 3)
    omega = f["frequency"] # (nqps, 3natom)


def parse_rlat(my):
    with open(my) as f:
        data = yaml.load(f)
        rlat = data["primitive_cell"]["reciprocal_lattice"]
    rlat = np.array(rlat)
    return rlat



def run():
    my = "/home/kazu/tmp13/phonopy-master/example/nacl/phonopy.yaml"
    rlat = parse_rlat(my)
    ra = (np.sum(rlat[0,:]*rlat[0,:]))**0.5
    print ra

    fn = "/home/kazu/tmp13/phonopy-master/example/nacl/mesh.hdf5"
    data = get_data(fn)



run()
