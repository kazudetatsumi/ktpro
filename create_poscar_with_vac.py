#i!/usr/bin/env python
import numpy as np
import spglib
import sys
sys.path.append("/home/kazu/ktpro")
#from change_h_pos_in_poscar_class import change_hpos as ch


class create_vacancy():
    def __init__(self, infile, target_atmidx):
        self.infile = infile
        self.target_atmidx = target_atmidx

    def GetSym(self):
        self.GetCrystalParamsFromPoscar()
        self.cell = (self.lattice, self.positions, self.numbers)
        self.sym = spglib.get_symmetry(self.cell, symprec=1e-5)

    def GetCrystalParamsFromPoscar(self):
        with open(self.infile, 'r') as f:
            self.lines = f.readlines()
        latmag = float(self.lines[1].split()[0])
        lattice = np.zeros((3, 3))
        for ivec in range(0, 3):
            lattice[ivec] = np.array((self.lines[ivec+2].split()[0:3]),
                                     dtype=float)
        self.lattice = lattice * latmag
        self.nspc = np.asarray((self.lines[6].split()), dtype=int)
        self.numbers = [i+1 for i in range(0, self.nspc.shape[0])
                        for j in range(0, self.nspc[i])]
        self.positions = np.zeros((np.sum(self.nspc), 3))
        for ipos in range(0, np.sum(self.nspc)):
            self.positions[ipos] = np.asarray((self.lines[8+ipos]
                                               .split()[0:3][0:4]
                                               ), dtype=float)

    def GetAllNonequivAtoms(self):
        vecids = []
        for atidx in range(np.sum(self.nspc[0:self.target_atmidx]),
                           np.sum(self.nspc[0:self.target_atmidx+1])):
            if self.sym['equivalent_atoms'][atidx] == atidx:
                vecids.append(atidx)
        self.vecids = np.array(vecids)
        print(self.vecids)

    def GeneratePoscarwithVac(self):
        lines = self.lines
        nspc = self.nspc
        nspc[self.target_atmidx] -= 1
        lines[6] = ""
        for sp in nspc:
            lines[6] += "  "
            lines[6] += str(sp)
        lines[6] += "\n"
        for idx, vecid in enumerate(self.vecids):
            with open('POSCAR_'+str(idx), 'w') as f:
                for il, line in enumerate(lines):
                    if il - 8 != vecid:
                        f.write(line)

