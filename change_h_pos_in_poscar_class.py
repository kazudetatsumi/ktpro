#!/usr/bin/env python
import numpy as np
import spglib
import matplotlib.pyplot as plt


class change_hpos():
    def __init__(self, infile, std, edgelength, nx, enefile=None, shift=None):
        self.infile = infile
        self.std = std
        self.edgelength = edgelength
        self.nx = nx
        self.enefile = enefile
        self.shift = shift

    def GetCrystalParamsFromPoscar(self):
        with open(self.infile, 'r') as f:
            self.lines = f.readlines()
        latmag = float(self.lines[1].split()[0])
        lattice = np.zeros((3, 3))
        for ivec in range(0, 3):
            lattice[ivec] = np.array((self.lines[ivec+2].split()[0:3]),
                                     dtype=float)
        self.lattice = lattice * latmag
        nspc = np.asarray((self.lines[6].split()), dtype=int)
        self.numbers = [i+1 for i in range(0, nspc.shape[0])
                        for j in range(0, nspc[i])]
        self.positions = np.zeros((np.sum(nspc), 3))
        for ipos in range(0, np.sum(nspc)):
            self.positions[ipos] = np.asarray((self.lines[9+ipos]
                                               .split()[0:3][0:4]
                                               ), dtype=float)

    def ShiftAllAtompos(self):
        if type(self.shift) != 'NoneType':
            print('All atom pos are shifted by', self.shift)
            pos_temp = self.cell[1]
            pos_shift = pos_temp - self.shift
            self.cell = (self.cell[0], pos_shift, self.cell[2])

    def GetRefineCell(self):
        lattice, positions, numbers = spglib.refine_cell(self.cell)
        sorted_numbers = np.sort(numbers)
        sorted_positions = positions[np.argsort(numbers)]
        sorted_positions[np.abs(sorted_positions) < 1e-10] = 0.
        self.cell = (lattice, sorted_positions, sorted_numbers)

    def GetSym(self):
        self.GetCrystalParamsFromPoscar()
        self.cell = (self.lattice, self.positions, self.numbers)
        self.ShiftAllAtompos()
        self.GetRefineCell()
        print(spglib.get_spacegroup(self.cell, symprec=1e-5))
        self.sym = spglib.get_symmetry(self.cell, symprec=1e-5)

    def GetAllHpos(self):
        hpos = []
        if type(self.edgelength) is list:
            dx = np.array(self.edgelength)/self.nx
            for ix in range(0, self.nx):
                for iy in range(0, self.nx):
                    for iz in range(0, self.nx):
                        hpos.append([(self.std - self.edgelength[0]/2
                                      + ix*dx[0]) % 1.0,
                                     (self.std - self.edgelength[1]/2
                                      + iy*dx[1]) % 1.0,
                                     (self.std - self.edgelength[2]/2
                                      + iz*dx[2]) % 1.0])
            self.hpos = np.array(hpos)
            print(self.hpos)
        else:
            print(self.edgelength)
            print(type(self.edgelength))
            # dx = self.edgelength/(self.nx-1)
            dx = self.edgelength/self.nx
            for ix in range(0, self.nx):
                for iy in range(0, self.nx):
                    for iz in range(0, self.nx):
                        hpos.append([(self.std - self.edgelength/2
                                      + ix*dx) % 1.0,
                                     (self.std - self.edgelength/2
                                      + iy*dx) % 1.0,
                                     (self.std - self.edgelength/2
                                      + iz*dx) % 1.0])
            self.hpos = np.array(hpos)
            print(self.hpos)

    def GetIrreducibleShift_old(self):
        irr_hpos = self.hpos[0].reshape((1, 3))
        for ih, _hpos in enumerate(self.hpos):
            for rot, trans in zip(self.sym['rotations'],
                                  self.sym['translations']):
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
        print(irr_hpos)
        print(irr_hpos.shape)
        self.irr_hpos = irr_hpos

    def GetIrreducibleShift(self):
        # print(sym['rotations'].shape)
        # print(hpos.shape)
        # print(np.matmul(hpos, sym['rotations']).transpose(1, 0, 2).shape)
        # print(sym['translations'].shape)
        # print((np.matmul(hpos, sym['rotations']).transpose(1, 0, 2)
        # + sym['translations']).shape)
        irr_hpos = np.array(self.hpos)
        isym = 0
        for rot, trans in zip(self.sym['rotations'][1:],
                              self.sym['translations'][1:]):
            isym += 1
            sym_hpos = ((np.matmul(rot, irr_hpos.T)).T + trans) % 1.0
            cond = np.ones((irr_hpos.shape[0]), dtype=bool)
            for i, _irr_hpos in enumerate(irr_hpos):
                if cond[i]:
                    cond *= np.sum(np.abs(_irr_hpos - sym_hpos) % 1.0,
                                   axis=1) >= 1e-5
                    cond[i] = True
            irr_hpos = irr_hpos[cond]
        print(irr_hpos.shape)
        print(irr_hpos)
        self.irr_hpos = irr_hpos

    def GetIrreducibleShift2(self):
        # this algorithm is too slow and buggy.
        test = (np.matmul(self.hpos,
                          self.sym['rotations'][1:]).transpose(1, 0, 2)
                + self.sym['translations'][1:]) % 1.0
        test2 = np.repeat(np.expand_dims(test, 2), self.hpos.shape[0], axis=2)
        cond = np.prod(np.sum(np.abs((self.hpos - test2) % 1.0),
                              axis=3) >= 1e-5, axis=1, dtype=bool)
        #irr_idx = [0]
        #for icol in range(1, cond.shape[0]):
        #    if np.prod(cond[0:icol, icol], dtype=bool):
        #        irr_idx.append(icol)
        #self.irr_hpos = self.hpos[irr_idx]

    def GetDataOverAllHpos(self):
        self.irr_idx = np.zeros((self.hpos.shape[0]), dtype=int)
        for ih, _hpos in enumerate(self.hpos):
            MatchFound = False
            for iridx, ir in enumerate(self.irr_hpos):
                for rot, trans in zip(self.sym['rotations'],
                                      self.sym['translations']):
                    if np.sum(np.abs(_hpos - ((np.matmul(rot, ir) + trans)
                                              % 1.0))) < 1e-5:
                        MatchFound = True
                        self.irr_idx[ih] = iridx
                        break
                if MatchFound:
                    break
            if not MatchFound:
                self.irr_idx[ih] = 999
        print(self.irr_idx)
        print(self.irr_idx.shape)

    def GetDataOverAllHpos2(self):
        self.irr_idx = np.zeros((self.hpos.shape[0]), dtype=int)
        for ih, _hpos in enumerate(self.hpos):
            cond = np.ones((self.irr_hpos.shape[0]), dtype=bool)
            for rot, trans in zip(self.sym['rotations'],
                                  self.sym['translations']):
                sym_irr_hpos = (np.matmul(self.irr_hpos, rot) + trans) % 1.0
                cond *= np.sum(np.abs((_hpos - sym_irr_hpos) % 1.0),
                               axis=1) >= 1e-5
            self.irr_idx[ih] = np.arange(0, self.irr_hpos.shape[0]
                                         )[np.invert(cond)]
        print(self.irr_idx.shape)

    def GetDataOverAllHpos3(self):
        self.irr_idx = np.zeros((self.hpos.shape[0]), dtype=int)
        sym_irr_hpos = (np.matmul(self.irr_hpos,
                                  self.sym['rotations']).transpose(1, 0, 2)
                        + self.sym['translations']) % 1.0
        for ih, _hpos in enumerate(self.hpos):
            cond = np.prod(np.sum(np.abs((_hpos - sym_irr_hpos) % 1.0),
                           axis=2) >= 1e-5, axis=1, dtype=bool)
            self.irr_idx[ih] = np.arange(0, self.irr_hpos.shape[0]
                                         )[np.invert(cond)]
        print(self.irr_idx)

    def GetDataOverAllHpos4(self):
        # print(self.irr_hpos.shape)
        # print(self.sym['rotations'].shape)
        # print(np.matmul(self.irr_hpos,
        #       self.sym['rotations']).transpose(1, 0, 2).shape)
        # print(np.matmul(self.sym['rotations'],
        #                 self.irr_hpos.T).transpose(2, 0, 1).shape)
        # print(self.sym['translations'].shape)
        # test = (np.matmul(self.irr_hpos,
        #                   self.sym['rotations']).transpose(1, 0, 2)
        #         + self.sym['translations']) % 1.0
        test = (np.matmul(self.sym['rotations'],
                          self.irr_hpos.T).transpose(2, 0, 1)
                + self.sym['translations']) % 1.0
        test2 = np.repeat(np.expand_dims(test, 2), self.hpos.shape[0], axis=2)
        cond = np.prod(np.sum(np.abs(self.hpos - test2) % 1.0,
                              axis=3) >= 1e-5, axis=1, dtype=bool)
        pairs = np.where(np.invert(cond))
        self.irr_idx = pairs[0][np.argsort(pairs[1])]
        print(self.irr_idx.shape)
        print(self.irr_idx)

    def GenerateShiftedPoscar(self):
        lat = self.cell[0]
        pos = self.cell[1]
        lines = self.lines
        for il in range(3):
            lines[il+2] = "     {:.16f}    {:.16f}    {:.16f} \n"\
                .format(lat[il, 0], lat[il, 1], lat[il, 2])
        for il in range(len(self.numbers)):
            if il <= 3:
                lines[il+9] = "  {:.16f}  {:.16f}  {:.16f}   F   F   F\n"\
                     .format(pos[il, 0], pos[il, 1], pos[il, 2])
            else:
                lines[il+9] = "  {:.16f}  {:.16f}  {:.16f}   T   T   T\n"\
                     .format(pos[il, 0], pos[il, 1], pos[il, 2])

        for iridx, ir in enumerate(self.irr_hpos):
            outfile = 'POSCAR_' + str(iridx+1).zfill(3)
            with open(outfile, 'w') as f:
                for il, line in enumerate(lines):
                    if il <= 8 or il >= 9 + self.positions.shape[0]:
                        f.write(line)
                    else:
                        if np.sum(np.abs(np.asarray(line.split()[0:3],
                                         dtype=float) - self.std)) < 1e-5:
                            f.write("  {:.16f}  {:.16f}  {:.16f}   F   F   F\n"
                                    .format(ir[0], ir[1], ir[2]))
                        else:
                            f.write(line)

    def GetEnergies(self):
        ene = []
        with open(self.enefile, 'r') as f:
            lines = f.readlines()
        for line in lines:
            ene.append(line.split()[3])
        self.ene = np.array(ene, dtype=float)
        print(np.min(self.ene))

    def GetPotential(self):
        self.potential = self.ene[self.irr_idx].reshape((self.nx, self.nx,
                                                         self.nx))

    def PlotPotential(self):
        #plt.pcolor(self.potential[:, :, 5] - np.min(self.potential))
        #plt.colorbar()
        plt.plot(self.potential[5,5,:]-np.min(self.potential))
        # plt.plot(self.hpos[:,0].reshape((self.nx, self.nx, self.nx))[:,7,7],
        # self.potential[:, 7, 7] - np.min(self.potential), marker='o')
        # y = np.zeros((self.nx))
        # x = np.zeros((self.nx))
        print(np.min(self.potential))
        print(self.potential[5, 5, 5])
        # for i in range(0, self.nx):
        #     y[i] = self.potential[i, i, 7] - np.min(self.potential)
        #     x[i] = ((self.hpos[:,0].reshape((self.nx, self.nx,
        #                                      self.nx))[i, i, 7])**2 +
        #             (self.hpos[:,1].reshape((self.nx, self.nx,
        #                                      self.nx))[i, i, 7])**2 +
        #             (self.hpos[:,2].reshape((self.nx, self.nx,
        #                                      self.nx))[i, i, 7])**2 )**0.5
        # plt.plot(x, y, marker='o')
        dist = ((self.hpos[:, 0].reshape((self.nx, self.nx, self.nx)) -
                 self.hpos[:, 0].reshape((self.nx, self.nx,
                                          self.nx))[7, 7, 7])**2 +
                (self.hpos[:, 1].reshape((self.nx, self.nx, self.nx)) -
                 self.hpos[:, 1].reshape((self.nx, self.nx,
                                          self.nx))[7, 7, 7])**2 +
                (self.hpos[:, 2].reshape((self.nx, self.nx, self.nx)) -
                 self.hpos[:, 2].reshape((self.nx, self.nx,
                                          self.nx))[7, 7, 7])**2)**0.5
        #plt.plot([dist[i, i, 7] for i in range(self.nx-1, 6, -1)],
        #         [self.potential[i, i, 7] - np.min(self.potential)
        #         for i in range(self.nx-1, 6, -1)], marker='o', label='110')
        #plt.plot([dist[i, 7, 7] for i in range(self.nx-1, 6, -1)],
        #         [self.potential[i, 7, 7] - np.min(self.potential)
        #         for i in range(self.nx-1, 6, -1)], marker='x', label='100')
        #plt.plot([dist[i, i, i] for i in range(self.nx-1, 6, -1)],
        #         [self.potential[i, i, i] - np.min(self.potential)
        #         for i in range(self.nx-1, 6, -1)], marker='.', label='111')
        #plt.legend()

        plt.show()

    def WritePotential(self):
        np.savetxt('hpot.txt', self.potential.flatten())


def samplerun():
    infile = 'CONTCAR'
    std = 0.5
    edgelength = 0.4
    # edgelength = 0.6
    nx = 14
    # nx = 13
    prj = change_hpos(infile, std, edgelength, nx)
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
    prj = change_hpos(infile, std, edgelength, nx, enefile=enefile)
    prj.GetEnergies()
    prj.GetSym()
    prj.GetAllHpos()
    prj.GetIrreducibleShift()
    prj.GetDataOverAllHpos4()
    prj.GetPotential()
    #prj.PlotPotential()
    prj.WritePotential()


#samplerun()
