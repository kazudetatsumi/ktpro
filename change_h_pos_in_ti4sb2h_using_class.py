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
    edgelength = [0.4, 0.4, 0.133]
    # edgelength = 0.6
    nx = 10
    shift = np.array([0., 0., 0.0238267635235516])
    # nx = 13
    prj = ch(infile, std, edgelength, nx, shift=shift)
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
    edgelength = [0.4, 0.4, 0.133]
    # edgelength = 0.6
    nx = 10
    #  nx = 13
    shift = np.array([0., 0., 0.0238267635235516])
    enefile = 'ENERGIES'
    prj = ch(infile, std, edgelength, nx, enefile=enefile, shift=shift)
    prj.GetEnergies()
    prj.GetSym()
    prj.GetAllHpos()
    prj.GetIrreducibleShift()
    prj.GetDataOverAllHpos4()
    prj.GetPotential()
    prj.PlotPotential()
    #prj.WritePotential()


samplerun2()
