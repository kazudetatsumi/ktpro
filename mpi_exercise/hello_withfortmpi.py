#!/usr/bin/env python
import ctypes
import numpy as np
from mpi4py import MPI


lib = ctypes.CDLL("./hello_by_py.so")
lib.hello.argtypes = [ ctypes.POINTER(ctypes.c_int32) ]
lib.hello.restype = ctypes.c_void_p


comm = MPI.COMM_WORLD
comm = comm.py2f()
result = lib.hello( ctypes.byref(ctypes.c_int32(comm)) )
rank = MPI.COMM_WORLD.Get_rank()

if rank == 0:
    print(result)

