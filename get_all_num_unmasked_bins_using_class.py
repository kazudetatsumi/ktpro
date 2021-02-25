#!/usr/bin/env python
# This script generates a histogram file with optimized binwidths,
# if it doesnot exsists.
# Then, this script gets the number of unmasked bins in the histogram.
# For analysis on the effect of an additional mask on the optimized bin-widths.
# Kazuyoshi TATSUMI 2021/02/24
import sys
sys.path.append("/home/kazu/ktpro")
import get_all_num_unmasked_bins_class as ganubc


def run():
    parelent = "/home/kazu/desktop/200522/Ei24/fineq/"
    heads = [parelent + str(i*2) + "m/" for i in range(1, 14)]
    logfiles = [heads[i-1] + "std-" + str(i*2) + ".log" for i in range(1, 14)]
    nbins = [ganubc.SGenhist(_h, _f).nbin() for _h, _f in zip(heads, logfiles)]
    print(nbins)


run()
