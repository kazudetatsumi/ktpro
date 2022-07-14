#!/usr/bin/env python
import numpy as np
import spglib


def GetCrystalParamsFromPoscar(infile):
    with open(infile, 'r') as f:
        lines = f.readlines()
    latmag = float(lines[1].split()[0])
    lattice = np.zeros((3, 3))
    for ivec in range(0, 3):
        lattice[ivec] = np.array((lines[ivec+2].split()[0:3]), dtype=float)
    lattice = lattice * latmag
    nspc = np.asarray((lines[6].split()), dtype=int)
    numbers = [i+1 for i in range(0, nspc.shape[0]) for j in range(0, nspc[i])]
    positions = np.zeros((np.sum(nspc), 3))
    for ipos in range(0, np.sum(nspc)):
        positions[ipos] = np.asarray((lines[9+ipos].split()[0:3][0:4]),
                                     dtype=float)
    return(lines, lattice, positions, numbers)


def GetRefineCell(cell):
    lattice, positions, numbers = spglib.refine_cell(cell)
    sorted_numbers = np.sort(numbers)
    sorted_positions = positions[np.argsort(numbers)]
    cell = (lattice, sorted_positions, sorted_numbers)
    return(cell)


def GetSym(infile):
    lines, lattice, positions, numbers = GetCrystalParamsFromPoscar(infile)
    cell = (lattice, positions, numbers)
    cell = GetRefineCell(cell)
    print(spglib.get_spacegroup(cell, symprec=1e-5))
    sym = spglib.get_symmetry(cell, symprec=1e-5)
    return(sym)


def GetAllHpos(std, nx, edgelength):
    hpos = []
    dx = edgelength/nx
    for ix in range(0, nx):
        for iy in range(0, nx):
            for iz in range(0, nx):
                hpos.append([std - edgelength/2 + ix*dx,
                             std - edgelength/2 + iy*dx,
                             std - edgelength/2 + iz*dx])
    hpos = np.array(hpos)
    return(hpos)


def GetIrreducibleShift(sym, hpos):
    irr_hpos = hpos[0].reshape((1, 3))
    for _hpos in hpos:
        print(_hpos)
        for rot, trans in zip(sym['rotations'], sym['translations']):
            __hpos = np.matmul(rot, _hpos) + trans
            MatchFound = False
            for ir in irr_hpos:
                if np.sum(np.abs(__hpos % 1.0 - ir)) < 1e-5:
                    MatchFound = True
                    break
            if MatchFound:
                break
        if not MatchFound:
            irr_hpos = np.append(irr_hpos, (_hpos % 1.0).reshape((1, 3)),
                                 axis=0)
    return irr_hpos


def GetDataOverAllHpos(hpos, irr_hpos, sym):
    irr_idx = np.zeros((hpos.shape[0]), dtype=int)
    for ih, _hpos in enumerate(hpos):
        MatchFound = False
        for iridx, ir in enumerate(irr_hpos):
            for rot, trans in zip(sym['rotations'], sym['translations']):
                if np.sum(np.abs(_hpos - ((np.matmul(rot, ir) + trans) % 1.0)
                                 )) < 1e-5:
                    MatchFound = True
                    irr_idx[ih] = iridx
                    break
            if MatchFound:
                break
    return(irr_idx)


def GenerateShiftedPoscar(irr_hpos, std, infile):
    lines, lattice, positions, numbers = GetCrystalParamsFromPoscar(infile)
    cell = GetRefineCell((lattice, positions, numbers))
    lat = cell[0]
    pos = cell[1]
    for il in range(3):
        lines[il+2] = "     {:.16f}    {:.16f}    {:.16f} \n"\
            .format(lat[il, 0], lat[il, 1], lat[il, 2])
    for il in range(len(numbers)):
        lines[il+9] = "  {:.16f}  {:.16f}  {:.16f}   F   F   F\n"\
             .format(pos[il, 0], pos[il, 1], pos[il, 2])
    for iridx, ir in enumerate(irr_hpos):
        outfile = 'POSCAR_' + str(iridx+1).zfill(3)
        with open(outfile, 'w') as f:
            for il, line in enumerate(lines):
                if il <= 8 or il >= 9 + positions.shape[0]:
                    f.write(line)
                else:
                    if np.sum(np.abs(np.asarray(line.split()[0:3],
                                                dtype=float) - std)) < 1e-5:
                        f.write("  {:.16f}  {:.16f}  {:.16f}   F   F   F\n"
                                .format(ir[0], ir[1], ir[2]))
                    else:
                        f.write(line)


def run():
    infile = 'CONTCAR'
    sym = GetSym(infile)
    std = 0.5
    edgelength = 0.4
    nx = 16
    hpos = GetAllHpos(std, nx, edgelength)
    irr_hpos = GetIrreducibleShift(sym, hpos)
    GenerateShiftedPoscar(irr_hpos, std, infile)
    #irr_idx = GetDataOverAllHpos(hpos, irr_hpos, sym)


run()
