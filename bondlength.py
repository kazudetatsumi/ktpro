#!/usr/bin/env python
import numpy as np

max_length = 20
tolerance = 0.0001
nbins = 50
sigma = 0.10
asigma = 0.005

car = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/POSCAR"


def get_pos(filename):
    with open(filename) as f:
        latvec = []
        ele = []
        num = []
        pos = []
        elepos = []
        atomid = []
        j = 0
        for line in f:
            if j >= 2 and j <= 4:
                latvec.append([float(x) for x in line.split()])
            if j == 5:
                ele = [str(x) for x in line.split()]
            if j == 6:
                num = [int(x) for x in line.split()]
            if j >= 8 and j <= 7 + sum(num):
                pos.append([float(x) for x in line.split()])
                atomid.append(j)
            j += 1
        limatom = 0
        for numa, elea in zip(num, ele):
            for k in range(limatom, limatom+numa):
                elepos.append(elea)
            limatom += numa
        return np.array(latvec), np.array(ele), np.array(num), np.array(pos),\
            np.array(elepos), np.array(atomid)


def expand_lattice(posa, elea, atomid, nlat):
    expos = []
    exele = []
    exid = []
    for ix in range(-1*nlat[0], nlat[0]):
        for iy in range(-1*nlat[1], nlat[1]):
            for iz in range(-1*nlat[2], nlat[2]):
                for pos, ele, aid in zip(posa, elea, atomid):
                    expos.append([ix + pos[0], iy + pos[1], iz + pos[2]])
                    exele.append(ele)
                    exid.append(aid)
    return np.array(expos), np.array(exele), np.array(exid)


def count_bondlength(expos, exele, exid, pos, ele, atomid, lat, el):
    data = []
    cartdata = []
    for e1 in el:
        for e2 in el:
            print(e1, e2)
            tmpdata = []
            tmpcartdispos = []
            for pos1, ele1, id1 in zip(pos, ele, atomid):
                for pos2, ele2, id2 in zip(expos, exele, exid):
                    if ele1 == e1 and ele2 == e2:
                        dispos = pos2 - pos1
                        cartdispos = np.dot(dispos, lat)
                        length = np.linalg.norm(cartdispos)
                        if length <= max_length and length > 0:
                            if ele1 == ele2 and pos2[0] >= 0 and pos2[0] <= 1\
                               and pos2[1] >= 0 and pos2[1] <= 1 and\
                               pos2[2] >= 0 and pos2[2] <= 1:
                                if id2 > id1:
                                    tmpdata.append(length)
                                    tmpcartdispos.append(cartdispos)
                            else:
                                tmpdata.append(length)
                                tmpcartdispos.append(cartdispos)
            data.append(tmpdata)
            cartdata.append(tmpcartdispos)
    return data, cartdata


def run(filename, numlat, rmin, rmax):
    lattice, element, numatom, posatom, eleatom, atomid = get_pos(filename)
    exposatom, exeleatom, exatomid = expand_lattice(posatom, eleatom, atomid,
                                                    numlat)
    distdata, cartposdata = count_bondlength(exposatom, exeleatom, exatomid,
                                             posatom, eleatom, atomid, lattice,
                                             element)


run(car, [2, 2, 3], 0.0, 4.)
