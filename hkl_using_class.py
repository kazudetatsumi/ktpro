#!/usr/bin/env python
# Kazuyoshi TATSUMI 2021/09/04
import sys
print("CHK")
sys.path.append("/home/kazu/ktpro")
import hkl_class as hc


def run(atomd):
    datadir = "/home/kazu/desktop/210904/"
    hkl_file = datadir + "H2.txt"
    prot_file = datadir + "prot.mem"
    out_file = datadir + "out_" + str(atomd) + ".mem"
    imag = 100.0
    nmag = 0.01
    proj = hc.hkl(hkl_file, prot_file, out_file, 5.0, atomd, 0.2, imag, nmag)
    proj.read_hklfile()
    proj.calc_strfac()
    proj.output()


def runsingle(atomw2):
    datadir = "/home/kazu/desktop/210904/"
    hkl_file = datadir + "H2.txt"
    prot_file = datadir + "prot.mem"
    out_file = datadir + "out_single_" + str(atomw2) + ".mem"
    imag = 100.0
    nmag = 0.01
    print("CHK1")
    proj = hc.hkl(hkl_file, prot_file, out_file, 5.0, 1.1, 0.2, imag, nmag,
                  atomw2=atomw2)
    proj.read_hklfile()
    proj.calc_strfac()
    proj.output()


runsingle(0.2)
#run(1.1)
