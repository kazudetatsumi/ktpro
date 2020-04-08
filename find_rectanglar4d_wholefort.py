#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
lib = ctypes.CDLL("/home/kazu/ktpro/rectanglar.so")
def calc_rectanglar_f90(condition):
    class result(ctypes.Structure):
        _fields_ = [("lb_0", ctypes.c_int), ("ub_0", ctypes.c_int), ("lb_1",
                    ctypes.c_int), ("ub_1", ctypes.c_int), ("lb_2",
                    ctypes.c_int), ("ub_2", ctypes.c_int)]
    lib.rectanglar.restype = result
    lib.rectanglar.argtypes = [ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               #np.ctypeslib.ndpointer(dtype=np.intat64, ndim=3)]
                               #np.ctypeslib.ndpointer(dtype=np.int32, ndim=3)]
                               np.ctypeslib.ndpointer(dtype=np.int32, ndim=4)]
    Nmax0 = condition.shape[0]
    Nmax1 = condition.shape[1]
    Nmax2 = condition.shape[2]
    Nmax3 = condition.shape[3]
    result = lib.rectanglar(
                        ctypes.byref(ctypes.c_int(Nmax0)),
                        ctypes.byref(ctypes.c_int(Nmax1)),
                        ctypes.byref(ctypes.c_int(Nmax2)),
                        ctypes.byref(ctypes.c_int(Nmax3)),
                        condition)
    result_lb_0 = result.lb_0
    result_ub_0 = result.ub_0
    result_lb_1 = result.lb_1
    result_ub_1 = result.ub_1
    result_lb_2 = result.lb_2
    result_lb_2 = result.ub_2
    lb = np.array([lb_0, lb_1, lb_2])
    ub = np.array([ub_0, ub_1, ub_2])
    return lb, ub


def run():
    maskfile = "/home/kazu/desktop/200204/coarse/hourbyhour/1h/out_hw_all.hdf5"
    f = h5py.File(maskfile)
    #condition = np.array(f["condition"][:,:,:,35], dtype=np.int32)
    condition = np.array(f["condition"][:,:,:,:], dtype=np.int32)
    lb, ub = calc_rectanglar_f90(np.squeeze(condition))

    #print lb
    #print ub 

    #print "---delete_array_pointer---"
    #lib.delete_array_pointer.restype = None
    #lib.delete_array_pointer.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
    #                             np.ctypeslib.ndpointer(dtype=np.float64, ndim=4), np.ctypeslib.ndpointer(dtype=np.float64, ndim=4), np.ctypeslib.ndpointer(dtype=np.float64, ndim=4)]
    #lib.delete_array_pointer(ctypes.byref(ctypes.c_int(lb_0)), ctypes.byref(ctypes.c_int(ub_0)), ctypes.byref(ctypes.c_int(lb_1)), ctypes.byref(ctypes.c_int(ub_1))
    #                         ctypes.byref(ctypes.c_int(lb_2)), ctypes.byref(ctypes.c_int(ub_2)))




run()
#lib.delete_array.restype = None
#lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]
#lib.delete_array(ctypes.byref(ctypes.c_int(klen)), k)
#plt.show()



