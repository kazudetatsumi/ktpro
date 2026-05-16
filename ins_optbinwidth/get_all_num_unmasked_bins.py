#!/usr/bin/env python
# This script generates a histogram file with optimized binwidths,
# if it doesnot exsists.
# Then, this script gets the number of unmasked bins in the histogram.
# For analysis on the effect of an additional mask on the optimized bin-widths.
# Kazuyoshi TATSUMI 2021/02/24
import sys
sys.path.append("/home/kazu/ktpro")
import gather_optbinindx_class as goc
import gethist4d_class as ghc
import count_unmasked_elements_class as cuec


class SGenhist(ghc.gethist4d_class, goc.gather_optbinidx,
               cuec.Count_Unmasked_Elements):
    def __init__(self, head, logfile):
        self.nw = super(SGenhist, self).getbinidx(logfile)
        super(SGenhist, self).histprocess(head)
        self.file = head + "hist_eliminated.hdf5"
        super(SGenhist, self).get_arrays()

    def nbin(self):
        return super(SGenhist, self).get_size()


def samplerun():
    parelent = "/home/kazu/desktop/200522/Ei24/fineq/"
    heads = [parelent + str(i*2) + "m/" for i in range(1, 14)]
    logfiles = [heads[i-1] + "std-" + str(i*2) + ".log" for i in range(1, 14)]
    nbins = [SGenhist(_h, _f).nbin() for _h, _f in zip(heads, logfiles)]
    print(nbins)


#samplerun()
