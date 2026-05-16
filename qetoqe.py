#!/usr/bin/env python
#import xml.etree.ElementTree as ET
from lxml import etree 
import numpy as np
import sys

atoangs = 0.529177210903
tree = etree.parse('./pwscf.save/data-file-schema.xml')
root = tree.getroot()
outfile = "outfile.in"

species2 = tree.xpath('//output/atomic_species')
NumEle = int(species2[0].attrib["ntyp"])
AM = tree.xpath('//output/atomic_species/species/mass')
pseu = tree.xpath('//output/atomic_species/species/pseudo_file')
species = tree.xpath('//output/atomic_species/species')
Elename = []
for i in species:
    Elename.append(i.attrib["name"]) 

cell = tree.findall('//output/atomic_structure/cell/')
values0 = cell[0].text.split()
values1 = cell[1].text.split()
values2 = cell[2].text.split()
latmat = np.array([values0, values1, values2]).astype(np.float)
inv_latmat = np.linalg.inv(latmat)
latmat_angs = latmat *  atoangs

atoms = tree.xpath('//output/atomic_structure/atomic_positions/atom')
Ele = []
vecs = []
for i in atoms:
    Ele.append(i.attrib["name"])
    vecs.append(i.text.split())
vecs = np.dot(np.array(vecs).astype(np.float),inv_latmat)

NumAtom = len(Ele)


#OUTPUT
if len(sys.argv) >= 2:
    if sys.argv[1] == "-vc":
        with open(outfile, mode='w') as f:
            f.write("&control\n")
            f.write("    calculation = 'vc-relax'\n")
            f.write("    tprnfor = .true.\n")
            f.write("    tstress = .true.\n")
            f.write("    pseudo_dir = '/home/kazu/code/q-e-qe-6.4.1/pseudo/'\n")
            f.write("/\n")
            f.write("&system\n")
            f.write("    ibrav = 0\n")
            f.write("    nat = %d\n" % (NumAtom))
            f.write("    ntyp = %d \n" % (NumEle))
            f.write("    ecutwfc = 70.0\n")
            f.write("    occupations='smearing',\n")
            f.write("    smearing='mv',\n")
            f.write("    degauss=0.01,\n")
            f.write("/\n")
            f.write("&electrons\n")
            f.write("    diagonalization = 'david'\n")
            f.write("    conv_thr = 1.0d-9\n")
            f.write("    mixing_beta = 0.1\n")
            f.write("    max_electronstep = 400\n")
            f.write("    mixing_ndim = 45\n")
            f.write("/\n")
            f.write("&ions\n")
            f.write("/\n")
            f.write("&cell\n")
            f.write("/\n")
            f.write("ATOMIC_SPECIES\n")
            for e, m, p in zip(Elename, AM. pseu):
                f.write(" %s %12.8f  %s\n" % (e, m, p))
            f.write("CELL_PARAMETERS angstrom \n")
            f.write("%22.16f%22.16f%22.16f \n" % (latmat_angs[0,0], latmat_angs[1,0], latmat_angs[2,0]))
            f.write("%22.16f%22.16f%22.16f \n" % (latmat_angs[0,1], latmat_angs[1,1], latmat_angs[2,1]))
            f.write("%22.16f%22.16f%22.16f \n" % (latmat_angs[0,2], latmat_angs[1,2], latmat_angs[2,2]))
            f.write("K_POINTS automatic\n")
            f.write(" 5 5 3 0 0 0\n")
else:
    print "FUCK"
