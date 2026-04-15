#!/usr/bin/env python
import numpy as np
import sys
import copy
import subprocess
from pypolymlp.core.interface_vasp import Poscar
from pypolymlp.utils.structure_utils import supercell_diagonal
from pypolymlp.utils.vasp_utils import write_poscar_file

np.random.seed(100)


def get_displacedpos(axis, atompos, atomtyp, dispmax_cart):
    atomposdisped = np.zeros_like(atompos.T)
    for aidx, (pos, typ) in enumerate(zip(atompos.T, atomtyp)):
        disp_cart = np.random.uniform(
                low=-0.5, high=0.5, size=3)*dispmax_cart[typ]
        # print(disp_cart)
        disp_frac = np.matmul(np.linalg.inv(axis), disp_cart)
        # print(disp_frac)
        # print("disp, dispmax:", np.sum(disp_cart**2)**0.5, dispmax_cart[typ])
        atomposdisped[aidx] = pos + disp_frac
    return atomposdisped.T

# file = sys.argv[1]


file = 'POSCAR'
st = Poscar(file).structure
st_sup = supercell_diagonal(st, [1, 1, 1])
# st_sup.elements[st_sup.elements == 'A'] = 'Pd'
# st_sup.elements[st_sup.elements == 'B'] = 'H'
st_sup_dispd = copy.deepcopy(st_sup)
latvec = st_sup.axis
latconst = np.sum(latvec**2, axis=0)**0.5
print(latconst)
# The length of h displ in my previous work is 1.32 angs.
dispmax_cart = np.array([1.32*0.5, 1.32])  # in Angs.
ndisps = 10
no = 0
#print(np.linspace(dispmax_cart[0], 0, 101)[:-1])
#print(dispmax_cart[0]*np.exp(-np.linspace(0, 100, 101)*np.log(1.5))[:-1])
for (dispmax1_cart, dispmax2_cart) in zip(
        np.linspace(dispmax_cart[0], 0, 101)[:-1],
        np.linspace(dispmax_cart[1], 0, 101)[:-1]):
    for didxo in range(ndisps):
        print(dispmax1_cart, dispmax2_cart)
        st_sup_dispd.positions = get_displacedpos(
                st_sup.axis,
                st_sup.positions,
                st_sup.types,
                np.array([dispmax1_cart, dispmax2_cart]))
        write_poscar_file(st_sup_dispd, filename="POSCAR_"+str(no+1).zfill(4))
        # replace "A" and "B" elements in generated POSCAR_ by "Pd" and "H",
        # respectively.
        subprocess.run('sed -e s/A/Pd/g POSCAR_'+str(no+1).zfill(4) +
                       ' > POSCAR_0'+str(no+1).zfill(4), shell=True)
        subprocess.run('sed -e s/B/H/g POSCAR_0'+str(no+1).zfill(4) +
                       ' > POSCAR_'+str(no+1).zfill(4), shell=True)
        subprocess.run('rm POSCAR_0'+str(no+1).zfill(4), shell=True)
        no += 1
