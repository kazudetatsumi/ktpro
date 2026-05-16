#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import get_all_maxfrac_class as gamc


def run():
    tail = "h"
    cond = "condparam09"
    proj = gamc.maxfracall(tail)
    proj.get_all_maxfrac()
    proj = gamc.maxfracall(tail, cond)
    proj.get_all_maxfrac()


run()
