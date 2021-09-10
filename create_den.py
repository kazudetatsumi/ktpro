#!/usr/bin/env python
import numpy as np


class create_den:
    def __init__(self, outfile, protfile, mesh, latconst, atomd, atomw):
        self.outfile = outfile
        self.protfile = protfile
        self.mesh = mesh
        self.latconst = latconst
        self.atomd = atomd
        self.atomw = atomw

    def get_all_coord(self):
        tmpx = np.arange(0, self.mesh/2+1)
        Y, Z, X = np.meshgrid(tmpx, tmpx, tmpx) # I dont know why this order is.
        _X = X[X>=Y]
        _Y = Y[X>=Y]
        _Z = Z[X>=Y]
        self.X = np.ravel(_X)
        self.Y = np.ravel(_Y)
        self.Z = np.ravel(_Z)
        print(self.X)
        print(self.Y)
        print(self.Z)

    def get_density(self):
        x = self.X*1.0/self.mesh*self.latconst
        y = self.Y*1.0/self.mesh*self.latconst
        z = self.Z*1.0/self.mesh*self.latconst
        d = self.atomd
        w = self.atomw
        self.den = (np.pi*w**2)**(-1.5)*(np.exp(-(z-d/2.0)**2/w**2 - y**2/w**2
                                         - x**2/w**2) +
                                         np.exp(-(z+d/2.0)**2/w**2 - y**2/w**2
                                         - x**2/w**2))

    def output(self):
        dummy = 0.0
        with open(self.protfile, mode='r') as g:
            lines = g.readlines()
        with open(self.outfile, mode='w') as f:
            for line in lines[0:2]:
                f.write("%s" % (line))
            f.write("%8d \n" % (self.den.shape[0]))
            for X, Y, Z, d in zip(self.X, self.Y, self.Z, self.den):
                f.write("%5d%5d%5d%20.10e%20.10f \n" % (X, Y, Z, d, dummy))
            for line in lines[-11:]:
                f.write("%s" % (line))

    def checkdensity(self):
        with open('dnsity.den', mode='r') as g:
            lines = g.readlines()
        denlines = lines[3:-11]
        density = []
        for line in denlines:
            density.append(float(line.split()[3])*float(line.split()[5]))

        density = np.array(density)
        print("sum", np.sum(density))

    def checkdensity2(self):
        with open('outsave.grd', mode='r') as g:
        #with open('dnsitysave.grd', mode='r') as g:
            lines = g.readlines()
        density = []
        denlines = lines[3:]
        print(denlines[0])
        for line in denlines:
            for val in line[:-1].split():
                density.append(float(val))

        density = np.array(density)
        print(np.sum(density))
        print(np.sum(density)/128.0**3)
        
def samplerun():
    mesh = 128
    latconst = 10.0
    atomd = 1.1
    atomw = 0.2
    outfile = 'out.den'
    protfile = 'prot.den'
    proj = create_den(outfile, protfile, mesh, latconst, atomd, atomw)
    proj.get_all_coord()
    proj.get_density()
    proj.output()
    proj.checkdensity2()
    #proj.read_prot()


samplerun()



