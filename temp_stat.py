#!/usr/bin/env python
# This script read an OSZICAR file of a vasp MD run and print the average and std of the temperatures.
import os
import numpy as np


def get_tempstat():
    os.system('grep T= OSZICAR > temperature')
    data = np.genfromtxt('temperature', dtype=float)
    temp = data[:, 2]        #  [K]
    print("temperature average [K]:",np.average(temp), "temperature standard deviation [K]:",np.std(temp))


get_tempstat()
