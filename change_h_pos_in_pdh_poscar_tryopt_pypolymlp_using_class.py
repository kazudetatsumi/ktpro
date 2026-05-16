#!/usr/bin/env python
import numpy as np
from pypolymlp.api.pypolymlp_calc import PypolymlpCalc
from scipy.optimize import minimize, Bounds
import copy
import sys
sys.path.append("/home/kazu/ktpro")
from change_h_pos_in_poscar_class import change_hpos as ch


def getene_wstr(hposs, positions):
    pot = "/home/kazu/WORK/pypolymlp/pd32h1/polymlp.yaml"
    poscars = "/home/kazu/WORK/vasp/pd32h1/calc/CONTCAR"
    case = PypolymlpCalc(pot=pot, verbose=True, require_mlp=True)
    case.load_structures_from_files(poscars=poscars)
    ene = []
    #case.structures[0].positions[:, :-1] = positions.T 
    case.structures[0].positions[0, 15] = positions
    #case.structures[0].axis = rcell[0]
    for hpos in hposs:
        case.structures[0].positions[:, -1] = hpos
        ene.append(case.eval()[0][0])
    return np.array(ene)


def getf(hposs, density, positions):
    pot = "/home/kazu/WORK/pypolymlp/pd32h1/polymlp.yaml"
    poscars = "/home/kazu/WORK/vasp/pd32h1/calc/CONTCAR"
    case = PypolymlpCalc(pot=pot, verbose=True, require_mlp=True)
    case.load_structures_from_files(poscars=poscars)
    fx = []
    case.structures[0].positions[0, 15] = positions
    for hpos in hposs:
        case.structures[0].positions[:, -1] = hpos
        fx.append(case.eval()[1][0][0, 15])
    return np.sum(np.array(fx).reshape(density.shape)*density)/np.sum(density)*7.8968828318389797


def minimize_positions_simple(f, positions, prj):
    #N = positions.shape[0]
    x = positions.ravel()
    bounds = Bounds(np.zeros_like(x), np.ones_like(x))

    def fun(x):
        #pos = x.reshape(N, 3)
        pos = x
        return f(pos, prj)
    #res = minimize(fun, x, method="L-BFGS-B", bounds=bounds,
    #               options={"maxiter": 200, "ftol": 1e-9, "iprint": -1})
    res = minimize(fun, x, bounds=bounds,
            #jac=getf(prj.hpos, prj.GetDensity()[0], positions),
                   options={'xatol': 1e-4, 'fatol': 1e-6, 'maxiter': 500})
    return res


def do_opti(positions, prj):
    res = minimize_positions_simple(get_E0, positions, prj)
    print("success:", res.success, "f*:", res.fun)
    positions_opt = res.x.reshape(N, 3)


def get_E0(positions, prj, jac=False):
    print(positions)
    #print("MSE:", np.mean((positions - prj.positions[:-1])**2))
    #print(np.mean((prj.cell[1][:-1]-positions)**2))
    #prj.cell[1][:-1] = positions
    #prj.GetSym(refine=False)
    #prj.GetIrreducibleShift()
    #prj.CheckShortAtomDistanceFlag()
    #prj.ene = getene_wstr(prj.irr_hpos, prj.rcell)
    #prj.GetDataOverAllHpos4()
    #prj.GetPotential()
    prj.potential = getene_wstr(prj.hpos, positions).reshape((prj.nx, prj.nx, prj.nx))
    prj.GetVG()
    #prj.GetG()
    prj.GetH()
    prj.GetEigen(Issave=False, Compress=True)
    if jac:
        prj.f = getf(prj.hpos, prj.GetDensity()[0], positions),
    #print('spacings:', np.round(prj.E[0:15] - np.min(prj.E), 6))
    #return prj.E_comp[0] 


def dongara_get_E0(positions, prj, jac=False):
    get_E0(positions, prj, jac=jac)
    return prj.E_comp[0]

#prj = ch(infile, std, edgelength, nx, prim=prim, rg=rg)
#    prj.GetSym()
#    prj.GetAllHpos()
#    prj.GetIrreducibleShift()
#    prj.CheckShortAtomDistanceFlag()
#    prj.ene = getene(prj.irr_hpos)
#    prj.GetDataOverAllHpos4()
#    prj.GetPotential()
#    prj.GetVG()
#    prj.GetG()
#    prj.GetH()
#    prj.GetEigen(Issave=False, Compress=True)

def samplerun():
    infile = 'CONTCAR'
    std = 0.5
    edgelength = 0.32
    # edgelength = 0.6
    nx = 20
    #rg = 13
    # nx = 13
    prim = False
    prj = ch(infile, std, edgelength, nx, prim=prim)
    prj.GetCrystalParamsFromPoscar()
    #prj.GetSym(refine=False)
    prj.GetAllHpos()
    prj.GetG()
    positions = copy.deepcopy(prj.positions[:-1])
    positions[15, 0] += -0.1112663224481034296
    do_opti(positions[15,0], prj)
    #print(get_E0(positions, prj))
    #print(get_E0(positions, prj))

    
    #get_E0(prj.rcell[1], prj)
    #prj.GetIrreducibleShift()
    #prj.CheckShortAtomDistanceFlag()
    #prj.ene = getene_wstr(prj.irr_hpos, prj.rcell)

    #prj.GetDataOverAllHpos4()
    #prj.GetPotential()
    #prj.GetVG()
    #prj.GetG()
    #prj.GetH()
    #prj.GetEigen(Issave=False, Compress=True)
    #prj.GetDensity()
    #prj.GenerateShiftedPoscar()
    #prj.GetDataOverAllHpos4()
    # print(prj.irr_idx)
    # print(prj.irr_idx.shape)


def samplerun2():
    infile = 'CONTCAR'
    std = 0.5
    edgelength = 0.32
    # edgelength = 0.6
    nx = 20
    #  nx = 13
    enefile = 'ENERGIES'
    prim = False
    prj = ch(infile, std, edgelength, nx, enefile=enefile, prim=prim)
    prj.GetSym()
    prj.GetAllHpos()
    prj.GetIrreducibleShift()
    prj.CheckShortAtomDistanceFlag()
    prj.GetEnergies()
    prj.GetDataOverAllHpos4()
    prj.GetPotential()
    prj.PlotPotential()
    prj.WritePotential()
    prj.GetPotFile()
    #prj.GetVG()


def samplerun3():
    infile = 'CONTCAR'
    std = 0.5
    edgelength = 0.32
    nx = 20
    enefile = 'ENERGIES'
    rg = 13
    #a = 4.071
    prim = False

    #prj = ch(infile, std, edgelength, nx, enefile=enefile, rg=rg, a=a)
    prj = ch(infile, std, edgelength, nx, enefile=enefile, rg=rg, prim=prim)
    #prj.GetEnergies()
    prj.GetSym()
    prj.GetAllHpos()
    prj.GetIrreducibleShift()
    prj.CheckShortAtomDistanceFlag()
    prj.GetEnergies()
    prj.GetDataOverAllHpos4()
    prj.GetPotential()
    prj.GetVG()
    prj.GetG()
    prj.GetH()
    prj.GetEigen(Issave=True)
    prj.GetDensity()
    #for istate in range(10):
    #    prj.GetDensityFile(istate)
    #prj.GetWavefuncs()
    #prj.GetTransitionMatrix(np.array([1.0, 1.0, 1.0]), label="111")
    #prj.GetTransitionMatrix(np.array([1.0, 1.0, 0.0]), label="110")
    #prj.GetTransitionMatrix(np.array([0.0, 0.0, 2.0]), label="001")
    #plt.show()
    #prj.PlotPotential()
    #prj.WritePotential()


def samplerun4():
    infile = 'CONTCAR'
    std = 0.5
    edgelength = 0.64
    nx = 20
    enefile = 'ENERGIES'
    rg = 13
    #a = 4.071
    prim = False

    #prj = ch(infile, std, edgelength, nx, enefile=enefile, rg=rg, a=a)
    prj = ch(infile, std, edgelength, nx, enefile=enefile, rg=rg, prim=prim)
    prj.LoadEigen()
    prj.GetG()
    prj.GetDensity()
    for istate in range(10):
        prj.GetDensityFile(istate)
    #prj.GetWavefuncs()
    #prj.GetTransitionMatrix(np.array([1.0, 1.0, 1.0]), label="111")
    #prj.GetTransitionMatrix(np.array([1.0, 1.0, 0.0]), label="110")
    #prj.GetTransitionMatrix(np.array([0.0, 0.0, 2.0]), label="001")
    #plt.show()
    #prj.PlotPotential()
    #prj.WritePotential()


samplerun()
