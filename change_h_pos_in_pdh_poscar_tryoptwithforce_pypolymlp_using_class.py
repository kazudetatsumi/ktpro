#!/usr/bin/env python
import numpy as np
from pypolymlp.api.pypolymlp_calc import PypolymlpCalc
from scipy.optimize import minimize, Bounds
import copy
import sys
sys.path.append("/home/kazu/ktpro")
from change_h_pos_in_poscar_class import change_hpos as ch


def getene_wstr(hposs, positions, case):
    ene = []
    case.structures[0].positions[0, 15] = positions
    for hpos in hposs:
        case.structures[0].positions[:, -1] = hpos
        ene.append(case.eval()[0][0])
    return np.array(ene)


def do_opt(positions, prj, case):
    bounds = Bounds(np.zeros_like(positions), np.ones_like(positions))
    res = minimize(
        get_E0,
        positions,
        bounds=bounds,
        jac=get_grad,
        args=(prj, case,),
        method="L-BFGS-B",
        options={"maxiter": 100},
        )
    print("success:", res.success, "f*:", res.fun)


def get_E0(positions, prj, case):
    print('pos:', positions)
    if np.abs(positions-prj.positions[15, 0]) > 0.000000000000001:
        prj.positions[15, 0] = positions
        prj.potential = getene_wstr(prj.hpos, positions, case
                                    ).reshape((prj.nx, prj.nx, prj.nx))
        prj.GetVG()
        prj.GetH()
        prj.GetEigen(Issave=False, Compress=True)
    return prj.E_comp[0]


def get_grad(positions, prj, case):
    fx = []
    if np.abs(positions-prj.positions[15, 0]) > 0.000000000000001:
        prj.positions[15, 0] = positions
        prj.potential = getene_wstr(prj.hpos, positions, case
                                    ).reshape((prj.nx, prj.nx, prj.nx))
        prj.GetVG()
        prj.GetH()
        prj.GetEigen(Issave=False, Compress=True)
    prj.GetDensity()
    case.structures[0].positions[0, 15] = positions
    for hpos in prj.hpos:
        case.structures[0].positions[:, -1] = hpos
        fx.append(case.eval()[1][0][0, 15])
    f = -np.sum(np.array(fx).reshape(prj.densities[0].shape) *
                prj.densities[0])/np.sum(prj.densities[0]
                                         )*7.8968828318389797*1000.
    print('grad:', f, '(meV/frac)', f/1000./7.8968828318389797, '(eV/Angs)')
    return f


def samplerun():
    infile = 'CONTCAR'
    std = 0.5
    edgelength = 0.32
    nx = 20
    prim = False
    prj = ch(infile, std, edgelength, nx, prim=prim)
    prj.GetCrystalParamsFromPoscar()
    prj.GetAllHpos()
    prj.GetG()
    positions = copy.deepcopy(prj.positions[:-1])
    positions[15, 0] += -0.1112663224481034296
    pot = "/home/kazu/WORK/pypolymlp/pd32h1/polymlp.yaml"
    poscars = "/home/kazu/WORK/vasp/pd32h1/calc/CONTCAR"
    case = PypolymlpCalc(pot=pot, verbose=True, require_mlp=True)
    case.load_structures_from_files(poscars=poscars)
    do_opt(positions[15, 0], prj, case)


samplerun()
