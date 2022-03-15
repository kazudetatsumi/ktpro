#!/usr/bin/env python
import numpy as np
from ctypes import *
import os, sys
#import matplotlib.pyplot as plt
lib = CDLL(os.environ["HOME"]+"/ktpro/ssvkernel_f90.so")

def calc_ssvkernel_f90(x, tin, M0, winparam):
    lib.ssvk.restype = c_void_p
    lib.ssvk.argtypes = [
                        POINTER(c_int),
                        POINTER(c_double),
                        POINTER(c_int),
                        POINTER(c_int),
                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)
                        ]
    xsize = x.shape[0]
    tinsize = tin.shape[0]
    yopt = np.zeros((tinsize))
    optw = np.zeros((tinsize))
    lib.ssvk(
            c_int(M0),
            c_double(winparam),
            c_int(xsize),
            c_int(tinsize),
            x,
            tin,
            optw,
            yopt
            )
    #print('yopt=', yopt)
    #print('optw=',optw)
    return yopt, tin, optw

def calc_ssvkernel_f90_tst(x, tinsize):
    lib.ssvk.restype = c_void_p
    lib.ssvk.argtypes = [
                        POINTER(c_int),
                        POINTER(c_double),
                        POINTER(c_int),
                        POINTER(c_int),
                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)
                        ]
    M0 = 80
    winparam = 5.
    xsize = x.shape[0]
    yopt = np.zeros((tinsize))
    optw = np.zeros((tinsize))
    lib.ssvk(
            c_int(M0),
            c_double(winparam),
            c_int(xsize),
            c_int(tinsize),
            x,
            optw,
            yopt
            )
    #print('yopt=', yopt)
    #print('optw=',optw)
    return yopt, optw


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
    #tinsize = gettinsize(xdat)
    #print('tinsize=',tinsize)
    #yopt, optw = calc_ssvkernel_f90_tst(xdat, tinsize)
    #plt.plot(yopt)
    #plt.show()
    ssvkernel_fort(xdat)


#def ssvkernel_fort(x, tin=None, M=80, winparam=5, WinFunc='Gauss'):
def ssvkernel(x, tin=None, M=80, winparam=5, WinFunc='Gauss'):
    if tin is None:
        T = np.max(x) - np.min(x)
        dx = np.sort(np.diff(np.sort(x)))
        dt_samp = dx[np.nonzero(dx)][0]
        tin = np.linspace(np.min(x), np.max(x), int(min(np.ceil(T / dt_samp), 1e3)))
        t = tin
        x_ab = x[(x >= min(tin)) & (x <= max(tin))]
    else:
        T = np.max(x) - np.min(x)
        x_ab = x[(x >= min(tin)) & (x <= max(tin))]
        dx = np.sort(np.diff(np.sort(x)))
        dt_samp = dx[np.nonzero(dx)][0]
        if dt_samp > min(np.diff(tin)):
            t = np.linspace(min(tin), max(tin), min(np.ceil(T / dt_samp), 1e3))
            sys.exit("t is NOT tin, which is not implemented yet in the fortran code!")
        else:
            t = tin
    return calc_ssvkernel_f90(x, t, M, winparam)

#testrun()
