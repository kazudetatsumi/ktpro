#!/usr/bin/env python
# This script calculates an integrated squared error between pdf and histogram. 
# The most simple version.
# You should change deltax, deltay, deltaz, deltaw by hands and prepare pdf.hdf5 and hist.hdf5.
# Kazuyoshi TATSUMI 2020.
import numpy as np
import h5py


deltax=4
deltay=3
deltaz=4
deltaw=3
xi = 120
xe = 172
yi = 61
ye = 145
zi = 16
ze = 53
ei = 20
ee = 70
#tcount = 6./3.0*53*5.6e+04
#tcount = 2.*tcount


def read_h5py(outfile, term):
    f=h5py.File(outfile, 'r')
    return np.array(f[term])

def run():
    proddelta = deltax*deltay*deltaz*deltaw
    pdf = read_h5py("pdf.hdf5", "pdf")
    pdfs = pdf[xi:xe, yi:ye, zi:ze, ei:ee]
    datahist = read_h5py("hist.hdf5", "data4")
    condhist = read_h5py("hist.hdf5","condition")
    #pdfsize = pdf.shape
    #print(pdfs.shape)
    #print(datahist.shape)
    histsize  = datahist.shape
    #print(histsize)
    ise = 0.
    sumhist = 0.
    sumpdfs = 0.
    sumpdf = np.sum(pdf)
    print(np.sum(datahist))
    sumhist = np.sum(condhist[condhist == np.max(condhist)]/np.max(condhist)*1.0)
    sumdata = np.sum(datahist[condhist == np.max(condhist)])
    for hx in range(0, histsize[0]):
        for hy in range(0, histsize[1]):
            for hz in range(0, histsize[2]):
                for hw in range(0, histsize[3]):
                    if condhist[hx, hy, hz, hw] == np.max(condhist):
                        sumpdfs += np.sum(pdfs[hx*deltax:(hx+1)*deltax, hy*deltay:(hy+1)*deltay, hz*deltaz:(hz+1)*deltaz,hw*deltaw:(hw+1)*deltaw])
    print("test",sumdata/sumpdfs*sumpdf*proddelta)
    ##check the numbers
    #if abs(sumdata/sumpdfs*sumpdf*1.0-tcount)/tcount > 0.01:
    #    print("total count estimated from pdf is more differed by 1 % from the hard-coded tcount value")
    for hx in range(0, histsize[0]):
        for hy in range(0, histsize[1]):
            for hz in range(0, histsize[2]):
                for hw in range(0, histsize[3]):
                    if condhist[hx, hy, hz, hw] == np.max(condhist):
                        for px in range(hx*deltax, (hx+1)*deltax):
                            for py in range(hy*deltay, (hy+1)*deltay):
                                for pz in range(hz*deltaz, (hz+1)*deltaz):
                                    for pw in range(hw*deltaw, (hw+1)*deltaw):
                                        ise += (pdfs[px, py, pz, pw] - datahist[hx, hy, hz, hw]/(1.0*sumdata/sumpdfs*sumpdf*proddelta))**2
    ise = ise / (proddelta*sumhist*1.0)
    print(sumdata/sumpdfs*sumpdf*1.0)
    print("average_pdfs:",sumpdfs/((proddelta*sumhist*1.0)))
    print("ise",ise)
    print("sq_ise/average_pdf (%)", ise**0.5/(sumpdfs/((proddelta*sumhist*1.0)))*100)


run()
