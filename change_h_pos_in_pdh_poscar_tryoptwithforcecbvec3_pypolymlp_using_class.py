#!/usr/bin/env python
import numpy as np
from pypolymlp.api.pypolymlp_calc import PypolymlpCalc
from scipy.optimize import minimize, Bounds
import copy
import sys
import time
from dataclasses import dataclass
sys.path.append("/home/kazu/ktpro")
from change_h_pos_in_poscar_class import change_hpos as ch


@dataclass
class LogCache:
    f_meV: float | None = None
    g_meV_per_frac: float | None = None


def fun(x, prj, atomslist, case, log: LogCache):
    f = float(get_E0(x, prj, atomslist, case))
    log.f_meV = f
    return f


def jac(x, prj, atomslist, case, log: LogCache):
    timeb = time.time()
    g = get_grad(x, prj, atomslist, case)
    timea = time.time()
    print('computation time in get_grad:', timea-timeb)
    log.g_meV_per_frac = g
    return np.array(g, dtype=float)


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
        #print(f"[cycle {k:3d}] f = {f_meV:.6f} meV, "
        #      f"grad = {g_meV_per_frac:.6e} meV/frac, {g_eV_per_Angs:.6e} eV/Å")
        print(f"[cycles {k:3d}] f = {f_meV:.6f} meV",  g_meV_per_frac, "meV/frac", g_eV_per_Angs, "eV/Å")
    return cb


def do_opt(positions, atomslist, prj, case):
    x0 = positions[atomslist].flatten()
    #bounds = Bounds([0.0], [1.0]) 
    d_cart = 0.06
    bounds = [(x0[i] - d_cart, x0[i] + d_cart) for i in range(len(x0))]

    log = LogCache()
    cb = make_callback(log)

    res = minimize(
        fun, x0,
        jac=jac,  # jac is almost necessary
        args=(prj, atomslist, case, log),
        bounds=bounds,
        method="L-BFGS-B",
        callback=cb,
        options={'ftol': 1e-9, 'gtol': 1e-7, 'maxiter': 500},
    )
    print("success:", res.success, "f*:", res.fun, "nit:", res.nit)


def getene_wstr(hposs, positions, atomslist, case):
    case.structures[0].positions[:, atomslist] = positions.T
    new_structs = []
    for hpos in hposs:
        s = copy.deepcopy(case.structures[0])
        s.positions[:, -1] = hpos
        new_structs.append(s)
    case.structures = new_structs
    eval_out = case.eval()
    energies = np.asarray(eval_out[0])
    case.structures = case.structures[0]
    return energies


def _getene_wstr(hposs, positions, atomslist, case):
    ene = []
    case.structures[0].positions[:, atomslist] = positions.T
    for hpos in hposs:
        case.structures[0].positions[:, -1] = hpos
        ene.append(case.eval()[0][0])
    return np.array(ene)


def get_E0(positions, prj, atomslist, case):
    print('positions:', positions)
    positions = positions.reshape((-1, 3))
    if np.max(np.abs(positions - prj.positions[atomslist])) > 1e-15 or 'E_comp' not in dir(prj):
        print('**get_E0: processes from potential to GetEigen are to be done')
        prj.positions[atomslist, :] = positions
        prj.potential = getene_wstr(prj.hpos, positions, atomslist, case).reshape((prj.nx, prj.nx, prj.nx))
        prj.GetVG()
        prj.GetH()
        prj.GetEigen(Issave=False, Compress=True)
    else:
        print('**get_E0: processes from potential to GetEigen are skipped')
    return prj.E_comp[0]


def _get_grad(positions, prj, atomslist, case):
    fx = []
    positions = positions.reshape((-1, 3))
    if np.max(np.abs(positions - prj.positions[atomslist])) > 1e-15 or 'E_comp' not in dir(prj):
        print('**get_grad: processeses from potential to GetEigen are to be done')
        prj.positions[atomslist] = positions
        prj.potential = getene_wstr(prj.hpos, positions, atomslist, case).reshape((prj.nx, prj.nx, prj.nx))
        prj.GetVG()
        prj.GetH()
        prj.GetEigen(Issave=False, Compress=True)
        prj.GetDensity()
    else:
        print('**get_grad: processeses from potential to GetEigen are skipped')
        prj.GetDensity()
    case.structures[0].positions[:, atomslist] = positions.T
    for ih, hpos in enumerate(prj.hpos):
        case.structures[0].positions[:, -1] = hpos
        fx.append(case.eval()[1][0][:, atomslist])

    f_meV_per_frac = -np.sum(np.array(fx).T.reshape((-1,) + prj.densities[0].shape) * prj.densities[0], axis=(1,2,3)) / np.sum(prj.densities[0]) \
                     * 7.8968828318389797 * 1000.0
    print(f_meV_per_frac)
    return f_meV_per_frac


def get_grad(positions, prj, atomslist, case):
    fx = []
    positions = positions.reshape((-1, 3))
    if np.max(np.abs(positions - prj.positions[atomslist])) > 1e-15 or 'E_comp' not in dir(prj):
        print('**get_grad: processeses from potential to GetEigen are to be done')
        prj.positions[atomslist] = positions
        prj.potential = getene_wstr(prj.hpos, positions, atomslist, case).reshape((prj.nx, prj.nx, prj.nx))
        prj.GetVG()
        prj.GetH()
        prj.GetEigen(Issave=False, Compress=True)
        prj.GetDensity()
    else:
        print('**get_grad: processeses from potential to GetEigen are skipped')
        prj.GetDensity()
    case.structures[0].positions[:, atomslist] = positions.T
    new_structs = []
    for hpos in prj.hpos:
        s = copy.deepcopy(case.structures[0])  
        s.positions[:, -1] = hpos      
        new_structs.append(s)
    case.structures = new_structs  
    eval_out = case.eval() 
    for ih in range(prj.hpos.shape[0]):
        fx.append(eval_out[1][ih][:, atomslist])
    case.structures = case.structures[0]
    f_meV_per_frac = -np.sum(np.array(fx).T.reshape((-1,) + prj.densities[0].shape) * prj.densities[0], axis=(1,2,3)) / np.sum(prj.densities[0]) \
                     * 7.8968828318389797 * 1000.0
    print(f_meV_per_frac)
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
    #prj.positions[15, 1] += 0.1112663224481034296
    #prj.positions[15, 2] += 0.1112663224481034296
    positions = copy.deepcopy(prj.positions[:-1])
    #positions[15, 0] += -0.01112663224481034296
    #positions[15, 1] += 0.02112663224481034296
    #positions[15, 2] += -0.031126632244810342965
    pot = "/home/kazu/WORK/pypolymlp/pd32h1/polymlp.yaml"
    poscars = "/home/kazu/WORK/vasp/pd32h1/calc/CONTCAR"
    case = PypolymlpCalc(pot=pot, verbose=True, require_mlp=True)
    case.load_structures_from_files(poscars=poscars)
    #atomslist=[4, 8, 9, 10, 11, 20]
    atomslist=range(32)
    do_opt(positions, atomslist, prj, case)



samplerun()
