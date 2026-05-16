#!/usr/bin/env python
from ctypes import * 
import numpy as np
from mpi4py import MPI
fmath = np.ctypeslib.load_library("fmath_mpif90.so",".")
fmath.scal.argtypes = [
                      np.ctypeslib.ndpointer(dtype=np.float64), 
                      POINTER(c_double), 
                      POINTER(c_int64), 
                      POINTER(c_int32)
                      ]
fmath.scal.restype = c_void_p
comm = MPI.COMM_WORLD 
comm = comm.py2f()
n = 10000000
x = np.ones(n, dtype=np.float64) 
rank = MPI.COMM_WORLD.Get_rank()
if rank==0:
   print(np.sum(x))
alpha = c_double(0.00001)
len = c_int64(x.size)
fmath.scal(x, alpha, len, c_int32(comm))
if rank==0:
   print(np.sum(x))
