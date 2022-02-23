#!/usr/bin/env python
import numpy as np
from ctypes import *
import os
lib = CDLL(os.environ["HOME"]+"/ktpro/sskernel.so")


def calc_sskernel_f90(x, tinsize):
    lib.ssk.restype = c_void_p
    lib.ssk.argtypes = [
                        POINTER(c_double),
                        POINTER(c_int),
                        POINTER(c_int),
                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)
                        ]
    xsize = x.shape[0]
    yh = np.zeros((tinsize))
    optw = c_double()
    lib.ssk(
            byref(optw),
            c_int(xsize),
            c_int(tinsize),
            x,
            yh
            )
    print('yh=', yh)
    print('optw=',optw.value)

def gettinsize(x):
    T = np.max(x) - np.min(x)
    dx = np.sort(np.diff(np.sort(x)))
    dt_samp = dx[np.nonzero(dx)][0]
    return int(min(np.ceil(T / dt_samp), 1e3))

def testrun():
    xdat = np.array([4.37, 3.87, 4.00, 4.03, 3.50, 4.08, 2.25, 4.70, 1.73, 4.93, 1.73, 4.62, 3.43, 4.25, 1.68, 3.92, 3.68, 3.10, 4.03, 1.77,
    4.08, 1.75, 3.20, 1.85, 4.62, 1.97, 4.50, 3.92, 4.35, 2.33, 3.83, 1.88, 4.60, 1.80, 4.73, 1.77, 4.57, 1.85, 3.52, 4.00, 3.70,
    3.72, 4.25, 3.58, 3.80, 3.77, 3.75, 2.50, 4.50, 4.10, 3.70, 3.80, 3.43, 4.00, 2.27, 4.40, 4.05, 4.25, 3.33, 2.00, 4.33, 2.93,
    4.58, 1.90, 3.58, 3.73, 3.73, 1.82, 4.63, 3.50, 4.00, 3.67, 1.67, 4.60, 1.67, 4.00, 1.80, 4.42, 1.90, 4.63, 2.93, 3.50, 1.97,
    4.28, 1.83, 4.13, 1.83, 4.65, 4.20, 3.93, 4.33, 1.83, 4.53, 2.03, 4.18, 4.43, 4.07, 4.13, 3.95, 4.10, 2.72, 4.58, 1.90, 4.50,
    1.95, 4.83, 4.12])
    tinsize=gettinsize(xdat)
    print('tinsize=',tinsize)
    calc_sskernel_f90(xdat, tinsize)

testrun()
