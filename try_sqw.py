#!/usr/bin/env python
import numpy as np
from pypolymlp.api.pypolymlp_calc import PypolymlpCalc
from phono3py.file_IO import read_fc2_from_hdf5
from scipy.optimize import minimize, Bounds
import sys
from typing import Optional

import numpy as np
from phono3py.file_IO import read_fc2_from_hdf5
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
    temperatures = np.arange(0, 1000, 10)
    is_gamma_center = True
    with_eigenvectors = True

    ph = Phonopy(unitcell_ph, supercell_matrix)
    ph.force_constants = fc2
    ph.run_mesh(mesh, with_eigenvectors=with_eigenvectors,
                is_gamma_center=is_gamma_center)
    th = ph.run_thermal_properties(temperatures=temperatures)
    mesh_dict = ph.get_mesh_dict()
    q_points = mesh_dict['qpoints']
    weights = mesh_dict['weights']
    eigenvalues = mesh_dict['frequencies']*hplanc
    eigenvectors = mesh_dict['eigenvectors']
    Q = np.array([1., 0., 0.])
    temr = np.abs(np.dot(Q, eigenvectors[:, :, -3:]))**2 / eigenvalues *  coth(



samplerun()
