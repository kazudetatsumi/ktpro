#!/usr/bin/env python
import ctypes
from mpi4py import MPI

nmax = 32000000
pi = ctypes.c_double(0.0)

lib = ctypes.CDLL("./pifort.so")
lib.pif.argtypes = [ ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_int),
                    ctypes.POINTER(ctypes.c_double)
                   ]
lib.pif.restype = ctypes.c_void_p

comm = MPI.COMM_WORLD
comm = comm.py2f()
lib.pif( 
        ctypes.c_int(nmax),
        ctypes.c_int(comm),
        pi
       )
rank = MPI.COMM_WORLD.Get_rank()
if rank == 0:
    print(pi.value, "in python")
