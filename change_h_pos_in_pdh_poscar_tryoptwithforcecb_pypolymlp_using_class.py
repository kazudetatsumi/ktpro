#!/usr/bin/env python
import numpy as np
from pypolymlp.api.pypolymlp_calc import PypolymlpCalc
from scipy.optimize import minimize, Bounds
import copy
import sys
from dataclasses import dataclass
sys.path.append("/home/kazu/ktpro")
from change_h_pos_in_poscar_class import change_hpos as ch


@dataclass
class LogCache:
    f_meV: float | None = None
    g_meV_per_frac: float | None = None


def fun(x, prj, case, log: LogCache):
    f = float(get_E0(x, prj, case))
    log.f_meV = f
    return f


def jac(x, prj, case, log: LogCache):
    g = float(get_grad(x, prj, case))
    log.g_meV_per_frac = g
    return np.array([g], dtype=float)


def make_callback(log: LogCache):
    k = 0
    conv = 7.8968828318389797
    def cb(xk):
        nonlocal k
        k += 1
        f_meV = log.f_meV
        g_meV_per_frac = log.g_meV_per_frac
        if g_meV_per_frac is None:
            print(f"[cycle {k:3d}] f = {f_meV:.6f} meV, grad = (not computed yet)")
            return
        g_eV_per_Angs = g_meV_per_frac / 1000.0 / conv
        print(f"[cycle {k:3d}] f = {f_meV:.6f} meV, "
              f"grad = {g_meV_per_frac:.6e} meV/frac, {g_eV_per_Angs:.6e} eV/Å")
    return cb


def do_opt(positions, prj, case):
    x0 = positions
    bounds = Bounds([0.0], [1.0]) 

    log = LogCache()
    cb = make_callback(log)

    res = minimize(
        fun, x0,
        jac=jac,
        args=(prj, case, log),
        bounds=bounds,
        method="L-BFGS-B",
        callback=cb,
        options={'ftol': 1e-6, 'gtol': 1e-4, 'maxiter': 500},
    )
    print("success:", res.success, "f*:", res.fun, "nit:", res.nit)


def getene_wstr(hposs, positions, case):
    ene = []
    case.structures[0].positions[0, 15] = positions[0]
    for hpos in hposs:
        case.structures[0].positions[:, -1] = hpos
        ene.append(case.eval()[0][0]) 
    return np.array(ene)


def get_E0(positions, prj, case):
    if np.abs(positions - prj.positions[15, 0]) > 1e-15 or 'E_comp' not in dir(prj):
        print('**get_E0: processes from potential to GetEigen are to be done', positions, prj.positions[15,0])
        prj.positions[15, 0] = positions[0]
        prj.potential = getene_wstr(prj.hpos, positions, case).reshape((prj.nx, prj.nx, prj.nx))
        prj.GetVG()
        prj.GetH()
        prj.GetEigen(Issave=False, Compress=True)
    else:
        print('**get_E0: processes from potential to GetEigen are skipped')
    return prj.E_comp[0]


def get_grad(positions, prj, case):
    fx = []
    if np.abs(positions - prj.positions[15, 0]) > 1e-15 or 'E_comp' not in dir(prj):
        print('**get_grad: processeses from potential to GetEigen are to be done', positions,  prj.positions[15,0])
        prj.positions[15, 0] = positions[0]
        prj.potential = getene_wstr(prj.hpos, positions, case).reshape((prj.nx, prj.nx, prj.nx))
        prj.GetVG()
        prj.GetH()
        prj.GetEigen(Issave=False, Compress=True)
        prj.GetDensity()
    else:
        print('**get_grad: processeses from potential to GetEigen are skipped')
        prj.GetDensity()
    case.structures[0].positions[0, 15] = positions[0]
    for hpos in prj.hpos:
        case.structures[0].positions[:, -1] = hpos
        fx.append(case.eval()[1][0][0, 15])

    f_meV_per_frac = -np.sum(np.array(fx).reshape(prj.densities[0].shape) * prj.densities[0]) / np.sum(prj.densities[0]) \
                     * 7.8968828318389797 * 1000.0
    return f_meV_per_frac


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
    #prj.positions[15, 0] += -0.1112663224481034296
    positions = copy.deepcopy(prj.positions[:-1])
    #positions[15, 0] += -0.1112663224481034296
    pot = "/home/kazu/WORK/pypolymlp/pd32h1/polymlp.yaml"
    poscars = "/home/kazu/WORK/vasp/pd32h1/calc/CONTCAR"
    case = PypolymlpCalc(pot=pot, verbose=True, require_mlp=True)
    case.load_structures_from_files(poscars=poscars)
    do_opt(positions[15, 0:1], prj, case)


samplerun()
