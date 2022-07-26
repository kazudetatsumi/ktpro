#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
from change_h_pos_in_poscar_class import change_hpos as ch


def samplerun():
    infile = 'CONTCAR'
    std = 0.5
    edgelength = [0.6, 0.6, 0.164]
    # edgelength = 0.6
    nx = 20
    shift = [0., 0., 0.0238267635235516]
    hshift = [0., 0., -3]
    # nx = 13
    print(shift)
    prj = ch(infile, std, edgelength, nx, shift=shift, hshift=hshift)
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
    edgelength = [0.6, 0.6, 0.164]
    # edgelength = 0.6
    nx = 20
    #  nx = 13
    shift = [0., 0., 0.0238267635235516]
    hshift = [0., 0., -3]
    enefile = 'ENERGIES'
    prj = ch(infile, std, edgelength, nx, enefile=enefile, shift=shift,
             hshift=hshift)
    prj.GetEnergies()
    prj.GetSym()
    prj.GetAllHpos()
    prj.GetIrreducibleShift()
    prj.GetDataOverAllHpos4()
    prj.GetPotential()
    prj.PlotPotential()
    prj.WritePotential()


samplerun2()
