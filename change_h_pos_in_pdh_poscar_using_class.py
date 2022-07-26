#!/usr/bin/env python
import numpy as np
import spglib
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
from change_h_pos_in_poscar_class import change_hpos as ch



def samplerun():
    infile = 'CONTCAR'
    std = 0.5
    edgelength = 0.4
    # edgelength = 0.6
    nx = 5
    # nx = 13
    prj = ch(infile, std, edgelength, nx)
    prj.GetSym()
    prj.GetAllHpos()
    prj.GetIrreducibleShift()
    prj.GenerateShiftedPoscar()
    prj.GetDataOverAllHpos4()
    # print(prj.irr_idx)
    # print(prj.irr_idx.shape)


def samplerun2():
    infile = 'CONTCAR'
    std = 0.5
    edgelength = 0.4
    # edgelength = 0.6
    nx = 14
    #  nx = 13
    enefile = 'ENERGIES'
    prj = ch(infile, std, edgelength, nx, enefile=enefile)
    prj.GetEnergies()
    prj.GetSym()
    prj.GetAllHpos()
    prj.GetIrreducibleShift()
    prj.GetDataOverAllHpos4()
    prj.GetPotential()
    #prj.PlotPotential()
    #prj.WritePotential()
    prj.GetVG()


def samplerun3():
    infile = 'CONTCAR'
    std = 0.5
    edgelength = 0.4
    nx = 14
    enefile = 'ENERGIES'
    rg = 5
    a = 4.07

    prj = ch(infile, std, edgelength, nx, enefile=enefile, rg=rg, a=a)
    prj.GetEnergies()
    prj.GetSym()
    prj.GetAllHpos()
    prj.GetIrreducibleShift()
    prj.GetDataOverAllHpos4()
    prj.GetPotential()
    prj.GetVG()
    prj.GetG()
    prj.GetH()
    prj.GetEigen()
    #prj.PlotPotential()
    #prj.WritePotential()


samplerun3()
