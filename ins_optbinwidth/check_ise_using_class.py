#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
import check_ise_class as cic
import os



def run():
    timen = os.getcwd().split("/")[-1]
    tryn = os.getcwd().split("/")[-2].split("try")[-1]
    log = "./std-ise-" + tryn +"_" + timen + ".log"
    reflog = "./std-" + timen + ".log"
    prj = cic.SCheckiselog(log, reflog)
    if os.getcwd()[19:25] == "200903":
        if int(timen[:-1]) <= 2:
            shift = np.array([3, 2, 3, 0])
        else:
            shift = np.array([0, 0, 1, 0])
    else:
        shift = np.zeros(4, dtype=int)
    print("shift:", shift)
    prj.getiselogdata(shift)
    prj.plot_ise_dcount()

    plt.show()

run()
