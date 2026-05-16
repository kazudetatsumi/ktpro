#!/usr/bin/env python
import numpy as np
from pypolymlp.api.pypolymlp_calc import PypolymlpCalc
import sys
sys.path.append("/home/kazu/ktpro")
from change_h_pos_in_poscar_class import change_hpos as ch
import time

def getene(irr_hpos):
    pot = "/home/kazu/WORK/pypolymlp/pd32h1/polymlp.yaml"
    poscars = "/home/kazu/WORK/vasp/pd32h1/calc/CONTCAR"
    case = PypolymlpCalc(pot=pot, verbose=True, require_mlp=True)
    case.load_structures_from_files(poscars=poscars)
    ene = []
    for ir in irr_hpos:
        case.structures[0].positions[:, -1] = ir
        ene.append(case.eval()[0][0])
    return np.array(ene)

def samplerun():
    infile = 'CONTCAR'
    std = 0.5
    edgelength = 0.32
    # edgelength = 0.6
    nx = 20
    rg = 13
    # nx = 13
    prim = False
    time0 = time.time()
    prj = ch(infile, std, edgelength, nx, prim=prim, rg=rg)
    time1 = time.time()
    prj.GetSym()
    time2 = time.time()
    prj.GetAllHpos()
    time3 = time.time()
    prj.GetIrreducibleShift()
    time4 = time.time()
    prj.CheckShortAtomDistanceFlag()
    time5 = time.time()
    prj.ene = getene(prj.irr_hpos)
    time6 = time.time()
    prj.GetDataOverAllHpos4()
    time7 = time.time()
    prj.GetPotential()
    time8 = time.time()
    prj.GetVG()
    time9 = time.time()
    prj.GetG()
    time10 = time.time()
    prj.GetH()
    time11 = time.time()
    prj.GetEigen(Issave=False, Compress=True)
    time12 = time.time()
    print(time1 - time0) #5.269050598144531e-05
    print(time2 - time1) #0.030260562896728516
    print(time3 - time2) #0.022807836532592773
    print(time4 - time3) #2.940394878387451  GetIrreducibleShift
    print(time5 - time4) #0.04167318344116211
    print(time6 - time5) #0.9073753356933594
    print(time7 - time6) #8.644279718399048  GetDataOverAllHpos4
    print(time8 - time7) #6.67572021484375e-05
    print(time9 - time8) #0.00047016143798828125 
    print(time10- time9) #0.005171775817871094
    print(time11 - time10) #141.2989010810852 GetH - > 4.144586563110352
    print(time12 - time11) #22.466978788375854 GetEigen
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
