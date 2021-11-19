# Kazuyoshi TATSUMI 2021/09/04
import numpy as np
import re


class hkl:
    def __init__(self, hkl_file, prot_file, out_file, latconst, atomd, atomw,
                 imag, nmag, atomw2=None, single=False):
        self.hkl_file = hkl_file
        self.prot_file = prot_file
        self.out_file = out_file
        self.latconst = latconst
        self.atomd = atomd
        self.atomw = atomw
        self.atomw2 = atomw2
        self.imag = imag
        self.nmag = nmag
        self.single = single

    def read_hklfile(self):
        f = open(self.hkl_file)
        _h = []
        _k = []
        _l = []
        _m = []
        for line in f:
            if not re.compile(r'[A-z]').search(line):
                values = line.split()
                _h.append(values[0])
                _k.append(values[1])
                _l.append(values[2])
                _m.append(values[9])
        f.close()
        self.hklm = np.zeros((5, len(_m)))
        self.hklm[0] = _h
        self.hklm[1] = _k
        self.hklm[2] = _l
        self.hklm[3] = _m
        self.hklm = np.unique(self.hklm, axis=1)

    def calc_strfac(self):
        h = self.hklm[0]
        k = self.hklm[1]
        l = self.hklm[2]
        multi = self.hklm[3]
        a = self.latconst
        w = self.atomw
        d = self.atomd
        if self.single:
            print("single")
            w2 = self.atomw2
            self.hklm[4] = self.imag*np.exp(-np.pi**2*(w**2*(h**2+k**2) +
                                                       w2**2*l**2)/(a**2))
        else:
            print("twoatoms")
            self.hklm[4] = 2.0*self.imag*np.cos(np.pi*d*l/a) \
                * np.exp(-w**2*np.pi**2*(h**2+k**2+l**2)/(a**2))

        self.hklm = self.hklm[:, np.lexsort((l, k, h, h**2+k**2+l**2))]

    def output(self):
        dummy = 0.00
        dummy2 = 0.003
        with open(self.out_file, mode='w') as f:
            g = open(self.prot_file, mode='r')
            head = g.readlines()[0:4]
            for item in head:
                f.write(item)
            if self.single:
                f.write("%15.7f%15.7f%15.7f \n" % (self.imag, dummy, dummy2))
            else:
                f.write("%15.7f%15.7f%15.7f \n" % (2.0*self.imag, dummy, dummy2))
            f.write("%5d\n" % self.hklm.shape[1])
            for hklidx in range(0, self.hklm.shape[1]):
                f.write("%5d%5d%5d%15.7f%15.7f%15.7f \n" %
                        (self.hklm[0][hklidx], self.hklm[1][hklidx],
                            self.hklm[2][hklidx], self.hklm[4][hklidx] +
                            np.random.randn()*np.abs(self.hklm[4][hklidx]
                                                     * self.nmag),
                            dummy, np.abs(self.hklm[4][hklidx]*self.nmag))
                        )


def samplerun():
    datadir = "/home/kazu/desktop/210904/"
    hkl_file = datadir + "H2.txt"
    prot_file = datadir + "prot.mem"
    out_file = datadir + "out2.mem"
    imag = 100.0
    nmag = 0.01
    proj = hkl(hkl_file, prot_file, out_file, 10.0, 1.1, 0.2, imag, nmag)
    proj.read_hklfile()
    proj.calc_strfac()
    proj.output()


samplerun()
