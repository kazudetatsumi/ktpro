#!/usr/bin/env python
import numpy as np
from pypolymlp.api.pypolymlp_calc import PypolymlpCalc
from phono3py.file_IO import read_fc2_from_hdf5
from scipy.optimize import minimize, Bounds
import sys
from typing import Optional

from phonopy import Phonopy, PhonopyQHA

from pypolymlp.calculator.properties import Properties
#from pypolymlp.calculator.utils.phonon_utils import is_imaginary
from pypolymlp.core.data_format import PolymlpParams, PolymlpStructure
from pypolymlp.utils.phonopy_utils import (
    phonopy_cell_to_structure,
    structure_to_phonopy_cell,
)
from pypolymlp.utils.structure_utils import isotropic_volume_change
sys.path.append("/home/kazu/ktpro")


def coth(x):
    return np.cosh(x) / np.sinh(x)



def coth_from_bose(E_meV, kBT_meV):
    # x = E/(2kBT)
    x = E_meV / (2.0 * kBT_meV)
    # n = 1/(exp(E/kBT) - 1)
    n = np.zeros_like(E_meV)
    mask = (E_meV > 1e-12)
    n[mask] = 1.0 / (np.exp(E_meV[mask] / kBT_meV) - 1.0)
    return 2*n + 1.0


#def GetTransitionMatrix(q, nebin=3000):
#    a = np.arange(nmesh)
#    pos = np.array(np.meshgrid(a, a, a)).transpose((0, 2, 1, 3)
#                                                   ).reshape(3, -1)
#    arg = (-1.0*np.matmul(q, pos)*2.*np.pi*1.j/nmesh
#            #       ).reshape(-1, nmesh, nmesh, nmesh).squeeze()
#    #mat = np.conj(self.wavefuncs)*self.wavefuncs[istate]*np.exp(arg)
#    #self.sqw = np.abs(mat.reshape(mat.shape[0], -1).sum(axis=1))**2
#    ene = np.arange(0, nebin, 1)
#    spec = np.zeros(nebin)
#    for iw, s in enumerate(sqw):
#        dE = (self.E[iw+1] - self.E[istate])
#        sigma = dE*0.02
#        spec += s*np.exp(-(ene - dE)**2/sigma**2)
#        self.dataset = {}
#        self.dataset['ene'] = ene
#        self.dataset['spec'] = spec
#        self.dataset['E'] = self.E
#        self.dataset['sqw'] = self.sqw
#        with open("./savedata_" + label + ".pkl", 'wb') as f:
#            pickle.dump(self.dataset, f, 4)


def samplerun():
    hplanc = 4.135667696  # [meV], which is used for conversion of THz -> meV
    pot = "/home/kazu/WORK/pypolymlp/pdh/222/polymlp.yaml"
    case = PypolymlpCalc(pot=pot, verbose=True)
    case.load_structures_from_files(poscars='CONTCAR')
    fc2 = read_fc2_from_hdf5('fc2.hdf5')
    unitcell = case.structures[0]
    unitcell_ph = structure_to_phonopy_cell(unitcell)
    supercell_matrix = np.array([2, 2, 2])
    mesh = (10, 10, 10)
    #temperatures = np.arange(0, 1000, 10)
    temperature = 25 # K
    kb = 8.617333262*0.01 # meV/K
    is_gamma_center = True
    with_eigenvectors = True

    ph = Phonopy(unitcell_ph, supercell_matrix)
    ph.force_constants = fc2
    ph.run_mesh(mesh, with_eigenvectors=with_eigenvectors,
                is_gamma_center=is_gamma_center)
    #th = ph.run_thermal_properties(temperatures=temperatures)
    mesh_dict = ph.get_mesh_dict()
    q_points = mesh_dict['qpoints']
    weights = mesh_dict['weights']
    eigenvalues = mesh_dict['frequencies']*hplanc
    eigenvectors = mesh_dict['eigenvectors']
    Q = np.array([1., 0., 0.])
    terms = np.sum(np.abs((Q*eigenvectors[:, :, -3:]))**2, axis=-1) \
         / eigenvalues * coth_from_bose(eigenvalues, kb*temperature)
    terms *= weights[:, np.newaxis]
    ene = np.arange(300)
    spec = np.zeros(300)
    for (dE, t) in zip(eigenvalues.flatten(), terms.flatten()):
        sigma = dE*0.03
        spec += t*np.exp(-(ene - dE)**2/sigma**2)


    import matplotlib.pyplot as plt
    plt.plot(ene, spec)
    plt.show()




samplerun()
