#!/usr/bin/env python
import numpy as np
import spglib
import matplotlib.pyplot as plt


class change_hpos():
    def __init__(self, infile, std, edgelength, nx, enefile=None,
                 shift=[0., 0., 0.], hshift=[0., 0., 0.], rg=None):
        self.infile = infile
        self.std = std
        self.edgelength = edgelength
        self.nx = nx
        self.enefile = enefile
        self.shift = np.array(shift)
        # print(shift)
        self.hshift = np.array(hshift)
        self.rg = rg
        self.a0 = 0.5291772
        self.mh = 1836
        self.Eh = 27.211

    def GetCrystalParamsFromPoscar(self):
        with open(self.infile, 'r') as f:
            self.lines = f.readlines()
        latmag = float(self.lines[1].split()[0])
        lattice = np.zeros((3, 3))
        for ivec in range(0, 3):
            lattice[ivec] = np.array((self.lines[ivec+2].split()[0:3]),
                                     dtype=float)
        self.lattice = lattice * latmag
        #if type(self.edgelength) is list:
        #    self.a = np.sum(self.lattice[0, :]**2)**0.5/self.a0*self.edgelength[0]
        #    self.b = np.sum(self.lattice[1, :]**2)**0.5/self.a0*self.edgelength[1]
        #    self.c = np.sum(self.lattice[2, :]**2)**0.5/self.a0*self.edgelength[2]
        #    self.edgelengthina0 = np.sum(self.lattice**2, axis=1)**0.5 /\
        #            self.a0*np.array(self.edgelength)
        #else:
            #self.a = np.sum(self.lattice[0, :]**2)**0.5/self.a0*self.edgelength
            #self.b = np.sum(self.lattice[1, :]**2)**0.5/self.a0*self.edgelength
            #self.c = np.sum(self.lattice[2, :]**2)**0.5/self.a0*self.edgelength
        self.edgelengthina0 = np.sum(self.lattice**2, axis=1)**0.5/self.a0 *\
            np.array(self.edgelength)
        #print("chk, a=", self.a)
        #print("chk, b=", self.b)
        #print("chk, c=", self.c)
        print("chk, test", self.edgelengthina0)
        nspc = np.asarray((self.lines[6].split()), dtype=int)
        self.numbers = [i+1 for i in range(0, nspc.shape[0])
                        for j in range(0, nspc[i])]
        self.positions = np.zeros((np.sum(nspc), 3))
        for ipos in range(0, np.sum(nspc)):
            self.positions[ipos] = np.asarray((self.lines[9+ipos]
                                               .split()[0:3][0:4]
                                               ), dtype=float)

    def ShiftAllAtompos(self):
        if np.sum(np.abs(self.shift)) > 1e-5:
            print('All atom pos are shifted by', self.shift)
            pos_temp = self.cell[1]
            pos_shift = pos_temp - np.array(self.shift)
            self.cell = (self.cell[0], pos_shift, self.cell[2])
            self.positions = pos_shift
        # else:
        #     print('Not shifted!')

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
        # print(spglib.get_spacegroup(self.cell, symprec=1e-5))
        tmpsym = spglib.get_symmetry(self.cell, symprec=1e-5)
        # print(tmpsym['rotations'].shape)
        # for rot, trans in zip(tmpsym['rotations'], tmpsym['translations']):
        #         print(rot, trans)
        rpositions = np.zeros((self.positions.shape[0]-1, 3))
        rnumbers = []
        irp = 0
        for ipos, pos in enumerate(self.positions):
            if np.sum(np.abs(pos - self.std) % 1.0) >= 1e-5:
                rpositions[irp] = pos
                rnumbers.append(self.numbers[ipos])
                irp += 1
        # symmetry on the atoms without the target h atom
        rcell = (self.lattice, rpositions, rnumbers)
        # print(spglib.get_spacegroup(rcell, symprec=1e-5))
        self.sym = spglib.get_symmetry(rcell, symprec=1e-5)
        # print(self.sym['rotations'].shape)
        # for rot, trans in zip(self.sym['rotations'], self.sym['translations']):
        #         print(rot, trans)
        #self.sym = spglib.get_symmetry(self.cell, symprec=1e-5)

    def GetAllHpos(self):
        # Obtain all H positions in a 2D np array [atomidx, x/y/z/ axis],
        # whose 1st axis is in the C order, i.e., the z coordination is most
        # rapidly changed.
        hpos = []
        if type(self.edgelength) is list:
            dx = np.array(self.edgelength)/self.nx
            for ix in range(0, self.nx):
                for iy in range(0, self.nx):
                    for iz in range(0, self.nx):
                        hpos.append([(self.std - self.edgelength[0]/2
                                      + (ix+self.hshift[0])*dx[0]) % 1.0,
                                     (self.std - self.edgelength[1]/2
                                      + (iy+self.hshift[1])*dx[1]) % 1.0,
                                     (self.std - self.edgelength[2]/2
                                      + (iz+self.hshift[2])*dx[2]) % 1.0])
            self.hpos = np.array(hpos)
        else:
            # print(self.edgelength)
            # print(type(self.edgelength))
            # dx = self.edgelength/(self.nx-1)
            dx = self.edgelength/self.nx
            for ix in range(0, self.nx):
                for iy in range(0, self.nx):
                    for iz in range(0, self.nx):
                        hpos.append([(self.std - self.edgelength/2
                                      + (ix+self.hshift[0])*dx) % 1.0,
                                     (self.std - self.edgelength/2
                                      + (iy+self.hshift[1])*dx) % 1.0,
                                     (self.std - self.edgelength/2
                                      + (iz+self.hshift[2])*dx) % 1.0])
            self.hpos = np.array(hpos)
        #self.hpos += self.hshift
        # print(self.hpos)

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
        # print(irr_hpos)
        # print(irr_hpos.shape)
        self.irr_hpos = irr_hpos

    def GetIrreducibleShift(self):
        # Obtain irreducible H positions w.r.t. symmetry operations.
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
        # print(irr_hpos.shape)
        # print(irr_hpos)
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
        # Inverse process to the process 'IrreducibleShift'.
        # Here the irreducible H index for each of all the H positions is
        # clarified.
        # GetDataOverAllHpos4 is fastest.
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
        # print(self.irr_idx)
        # print(self.irr_idx.shape)

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
        # print(self.irr_idx.shape)

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
        # print(self.irr_idx)

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
        # print(self.irr_idx.shape)
        # print(self.irr_idx)

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
            if self.irr_hpos.shape[0] > 1000:
                outfile = 'POSCAR_' + str(iridx+1).zfill(4)
            else:
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
        # print(np.min(self.ene))

    def GetPotential(self):
        self.potential = self.ene[self.irr_idx].reshape((self.nx, self.nx,
                                                         self.nx))

    def PlotPotential(self):
        #plt.pcolor(self.potential[10, :, :] - np.min(self.potential))
        #plt.colorbar()
        plt.plot(self.potential[10,10,:]-np.min(self.potential))
        # plt.plot(self.hpos[:,0].reshape((self.nx, self.nx, self.nx))[:,7,7],
        # self.potential[:, 7, 7] - np.min(self.potential), marker='o')
        # y = np.zeros((self.nx))
        # x = np.zeros((self.nx))
        print(np.min(self.potential))
        print(np.unravel_index(np.argmin(self.potential), self.potential.shape))
        print(self.potential[10,10,12])
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
        # Write the potential in a column, which is in the Fortran order,
        # because Maple uses the Fortran order as its default.
        np.savetxt('hpot.txt', self.potential.flatten(order='F'))

    def GetVG(self):
        self.z = np.fft.fftn((self.potential-np.min(self.potential))
                             / self.Eh, norm='forward')

    def GetG(self):
        Gs = []
        if type(self.edgelength) is list:
            for i in range(-self.nx//2, self.nx//2+1):
                for j in range(-self.nx//2, self.nx//2+1):
                    for k in range(-self.nx//2, self.nx//2+1):
                        #if (i/self.edgelengthina0[0])**2 + (j/self.edgelengthina0[1])**2 + (k/self.edgelengthina0[2])**2 < self.rg**2:
                        if np.sum((np.array([i, j, k])/self.edgelengthina0)**2) < self.rg**2:
                            Gs.append([i, j, k])
        else:
            for i in range(-self.nx//2, self.nx//2+1):
                for j in range(-self.nx//2, self.nx//2+1):
                    for k in range(-self.nx//2, self.nx//2+1):
                        if i**2 + j**2 + k**2 < self.rg**2:
                            Gs.append([i, j, k])
        self.Gs = np.array(Gs).reshape((-1, 3))
        print('chk Gs:', self.Gs.shape)

    def GetH(self):
        self.H = np.zeros((self.Gs.shape[0], self.Gs.shape[0]), dtype='cdouble'
                          )
        for i, gi in enumerate(self.Gs):
            for j, gj in enumerate(self.Gs):
                if i == j:
                    _gi = gi/self.edgelengthina0
                    K = (2*np.pi)**2*np.dot(_gi, _gi)/(2.*self.mh)
                else:
                    K = 0.
                dg = (gi - gj) % self.nx
                V = self.z[dg[0], dg[1], dg[2]]
                self.H[i, j] = K + V

    def GetEigen(self):
        self.E, self.U = np.linalg.eigh(self.H)
        self.E *= self.Eh * 1000.
        #print(self.E[0:10] - np.min(self.E))
        #plt.plot(self.E[0:13] - np.min(self.E), marker='o')
        #plt.show()

    def GetWavefuncs(self):
        nmesh = self.nx*2
        a = np.arange(nmesh)
        pos = np.array(np.meshgrid(a, a, a)).transpose((0, 2, 1, 3)).reshape((3, -1))
        arg = 1.0*np.matmul(self.Gs, pos)*2.*np.pi*1.j/nmesh
        self.wavefuncs = np.matmul(self.U.T, np.exp(arg)).reshape(
                (-1, nmesh, nmesh, nmesh))

    def GetDensity(self):
        self.GetWavefuncs()
        self.densities = np.imag(self.wavefuncs)**2+np.real(self.wavefuncs)**2
        #plt.pcolor(self.densities[0,:,:,14])
        #plt.show()
        #phi = np.zeros((self.nx, self.nx), dtype='cdouble')
        #for ix in range(self.nx):
        #    for iy in range(self.nx):
        #        for ig in range(self.Gs.shape[0]):
        #            g = self.Gs[ig]
        #            phi[ix, iy] += self.U[ig, 0]*np.exp(-2.0*np.pi*1.j*(1.0*g[0]*ix+1.0*g[1]*iy+g[2]*7.0)/self.nx)
        #den = np.imag(phi)**2+np.real(phi)**2
        #plt.pcolor(den)
        #plt.show()

    def GetDensityFile(self, no):
        outfile = 'density' + str(no) + '.out'
        den = self.densities[no].flatten(order='F')
        out = ""
        for irow in range(den.shape[0] // 5):
            out += "{:.11e} {:.11e} {:.11e} {:.11e} {:.11e} \n".format(
                    den[irow*5], den[irow*5+1], den[irow*5+2], den[irow*5+3],
                    den[irow*5+4])
        last = ""
        for il in range(den.shape[0] % 5):
            last += "{:.11e} ".format(den[il+(den.shape[0] // 5)*5])
        out += last + "\n"
        with open(outfile, 'w') as f:
            f.write(out)

    #def GetTransitionMatrix(self, q, domegas, Isplot=True, label=None):
    def GetTransitionMatrix(self, q, Isplot=True, label=None):
        nmesh = self.nx*2
        a = np.arange(nmesh)
        pos = np.array(np.meshgrid(a, a, a)).transpose((0, 2, 1, 3)).reshape(3, -1)
        arg = (-1.0*np.matmul(q, pos)*2.*np.pi*1.j/nmesh).reshape(-1, nmesh, nmesh, nmesh).squeeze()
        #mat = np.conj(self.wavefuncs)*self.wavefuncs[0]*np.exp(np.repeat(np.expand_dims(arg, 1), self.wavefuncs.shape[0], axis=1))
        mat = np.conj(self.wavefuncs)*self.wavefuncs[0]*np.exp(arg)
        #mat2 = (mat.transpose((1, 2, 3, 4, 0))*domegas).transpose((0, 4, 1, 2, 3))
        self.sqw = np.abs(mat.reshape(mat.shape[0], -1).sum(axis=1))**2
        #self.sqw = (np.abs(mat.reshape(mat2.shape[0], mat2.shape[1], -1).sum(axis=1))**2).sum(axis=1)
        if Isplot:
            ene = np.arange(0, 3000, 1)
            spec = np.zeros(3000)
            for iw, s in enumerate(self.sqw[1:]):
                dE = (self.E[iw+1] - self.E[0])
                sigma = dE*0.02
                spec += s*np.exp(-(ene - dE)**2/sigma**2)
            plt.plot(ene, spec, label=label)
            plt.bar(self.E[1:] - self.E[0], self.sqw[1:])
            plt.xlim((0, 400))

    def Integrateoverallsolidangle(self, qlength):
        nmesh = 10
        thetas = np.pi*np.arange(0, nmesh)/nmesh
        dtheta = thetas[1] - thetas[0]
        phis = 2.0*np.pi*np.arange(0, nmesh*2)/nmesh*2
        dphi = phis[1] - phis[0]
        print(dtheta, dphi)
        print(thetas)
        print(phis)
        #nmesh = 10
        #a = np.arange(nmesh)
        #pos = np.array(np.meshgrid(a, a)).reshape(2, -1)
        #thetas = np.pi*pos[0]/nmesh
        #phis = 2.0*np.pi*pos[1]/nmesh
        #dtheta = np.pi/nmesh
        #dphi = 2.0*np.pi/nmesh
        #domegas = np.sin(thetas)*dtheta*dphi
        #qs = np.zeros((thetas.shape[0], 3))
        #qs[:,0] = np.sin(thetas)*np.cos(phis)
        #qs[:,1] = np.sin(thetas)*np.sin(phis)
        #qs[:,2] = np.cos(thetas)
        #self.GetTransitionMatrix(qs, domegas,  Isplot=True)


        for it, theta in enumerate(thetas):
            for ip, phi in enumerate(phis):
                print(it, ip)
                domega = np.sin(theta)*dtheta*dphi
                #q = qlength * np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
                q = qlength * np.array([np.cos(theta), np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi)])
                self.GetTransitionMatrix(q, Isplot=False)
                if it == 0 and ip == 0:
                    IntegratedSqw = np.zeros_like(self.sqw)
                IntegratedSqw += self.sqw*domega
        ene = np.arange(0, 3000, 1)
        spec = np.zeros(3000)
        for iw, s in enumerate(IntegratedSqw[1:]):
            dE = (self.E[iw+1] - self.E[0])
            sigma = dE*0.02
            spec += s*np.exp(-(ene - dE)**2/sigma**2)
        plt.plot(ene, spec)
        plt.bar(self.E[1:] - self.E[0], IntegratedSqw[1:])
        plt.xlim((0, 400))

    def Integrateoverallsolidangle2(self, qlength):
        nmesh = 10
        coss, weights = np.polynomial.legendre.leggauss(nmesh)
        thetas = np.arccos(coss)
        phis = 2.0*np.pi*np.arange(0, nmesh*2)/nmesh*2
        #dphi = phis[1] - phis[0]
        #print(dtheta, dphi)
        #print(thetas)
        #print(phis)
        #nmesh = 10
        #a = np.arange(nmesh)
        #pos = np.array(np.meshgrid(a, a)).reshape(2, -1)
        #thetas = np.pi*pos[0]/nmesh
        #phis = 2.0*np.pi*pos[1]/nmesh
        #dtheta = np.pi/nmesh
        #dphi = 2.0*np.pi/nmesh
        #domegas = np.sin(thetas)*dtheta*dphi
        #qs = np.zeros((thetas.shape[0], 3))
        #qs[:,0] = np.sin(thetas)*np.cos(phis)
        #qs[:,1] = np.sin(thetas)*np.sin(phis)
        #qs[:,2] = np.cos(thetas)
        #self.GetTransitionMatrix(qs, domegas,  Isplot=True)


        for it, theta in enumerate(thetas):
            for ip, phi in enumerate(phis):
                print(it, ip)
                #domega = np.sin(theta)*dtheta*dphi
                #q = qlength * np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
                q = qlength * np.array([np.cos(theta), np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi)])
                self.GetTransitionMatrix(q, Isplot=False)
                if it == 0 and ip == 0:
                    IntegratedSqw = np.zeros_like(self.sqw)
                IntegratedSqw += self.sqw*weights[it]
        ene = np.arange(0, 3000, 1)
        spec = np.zeros(3000)
        for iw, s in enumerate(IntegratedSqw[1:]):
            dE = (self.E[iw+1] - self.E[0])
            sigma = dE*0.02
            spec += s*np.exp(-(ene - dE)**2/sigma**2)
        plt.plot(ene, spec)
        plt.bar(self.E[1:] - self.E[0], IntegratedSqw[1:])
        plt.xlim((0, 400))


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
