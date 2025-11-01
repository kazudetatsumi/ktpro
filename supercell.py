#!/usr/bin/env python
import numpy as np
import sys
import copy
import subprocess
from pypolymlp.core.interface_vasp import Poscar
from pypolymlp.utils.structure_utils import supercell_diagonal
from pypolymlp.utils.vasp_utils import write_poscar_file

np.random.seed(100)


def get_displacedpos(atompos, atomtyp, dispmax):
    atomposdisped = np.zeros_like(atompos.T)
    for aidx, (pos, typ) in enumerate(zip(atompos.T, atomtyp)):
        disp = np.random.uniform(low=-0.5, high=0.5, size=3)*dispmax[typ]
        print(np.sum(disp**2)**0.5, dispmax[typ])
        atomposdisped[aidx] = pos + disp
    return atomposdisped.T

#file = sys.argv[1]


file = 'POSCAR'
st = Poscar(file).structure
st_sup = supercell_diagonal(st, [2, 2, 2])
#st_sup.elements[st_sup.elements == 'A'] = 'Pd'
#st_sup.elements[st_sup.elements == 'B'] = 'H'
st_sup_dispd = copy.deepcopy(st_sup)
latvec = st_sup.axis
latconst = np.sum(latvec**2, axis=0)**0.5
print(latconst[0])
dispmax = np.array([1.32*0.1, 1.32])/latconst[0]  # The length of h displ in my previous work is 1.32 angs.
ndisps = 10
#st_sup_dispd_list = []
#for dispmax1 in np.linspace(0, dispmax[0], 3)[1:]:
#    for dispmax2 in np.linspace(0, dispmax[1], 3)[1:]:
for (dispmax1, dispmax2) in zip(np.linspace(0, dispmax[0], 2)[1:],
                                np.linspace(0, dispmax[1], 2)[1:]):
    print(dispmax1, dispmax2)
    for no in range(ndisps):
        st_sup_dispd.positions = get_displacedpos(st_sup.positions,
                                                  st_sup.types,
                                                  np.array([dispmax1, dispmax2]))
        write_poscar_file(st_sup_dispd, filename="POSCAR_disp"+str(no))
        subprocess.run('sed -e s/A/Pd/g POSCAR_disp'+str(no)+' > POSCAR_disp0'+str(no), shell=True)
        subprocess.run('sed -e s/B/H/g POSCAR_disp0'+str(no)+' > POSCAR_disp'+str(no), shell=True)
        subprocess.run('rm POSCAR_disp0'+str(no), shell=True)
        #st_sup_dispd_list.append(copy.deepcopy(st_sup_dispd))


#write_poscar_file(st_sup, filename="SPOSCAR")






