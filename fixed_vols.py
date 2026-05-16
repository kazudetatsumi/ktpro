#!/usr/bin/env python

import numpy as np
import subprocess

fs = '/home/kazu/gamma-si3n4-unit/fixed_vols/MASTER_POSCAR'

dvol = 0.4*4
nmax = 5

def get_lat(filename):
    data = []
    j = 0
    with open(filename) as f:
        for line in f:
            if j >= 2 and j <= 4:
                data.append([float(x) for x in line.split()])
            j += 1 
    return np.array(data)

def creat_poscar(inum):
    latvec=get_lat(fs)
    volume = np.dot(np.cross(latvec[0,:],latvec[1,:]),latvec[2,:])
    latfrac=((volume+dvol*inum)/volume)**(1.0/3.0)
    cmd = "sed -e s/NUM/%s/ %s > POSCAR%s" % (str(latfrac), fs, str(inum))
    subprocess.call( cmd, shell=True )

def run():
    for i in range(nmax*(-1),nmax+1):
        creat_poscar(i)
    creat_poscar(-1*nmax-2)
    creat_poscar(-1*nmax-4)
    creat_poscar(nmax+9)
    creat_poscar(nmax+7)
    creat_poscar(nmax+5)
    creat_poscar(nmax+4)
    creat_poscar(nmax+3)
    creat_poscar(nmax+1)
    creat_poscar(nmax+2)

run()
