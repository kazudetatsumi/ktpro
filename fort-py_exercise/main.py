#!/usr/bin/env python
import ctypes
import numpy as np

def call_fortran(n, A):
    f = np.ctypeslib.load_library("libfort.so", ".")
    f.test_.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        np.ctypeslib.ndpointer(dtype=np.float64)
        ]
    f.test_.restype = ctypes.c_void_p

    fn = ctypes.byref(ctypes.c_int32(n))

    f.test_(fn, A)

print("** calling fortran from python")

n = 3
A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=np.float64)

print("** write from python")
print(n)
print(A[0,0], A[0,1], A[0,2])
print(A[1,0], A[1,1], A[1,2])
print(A[2,0], A[2,1], A[2,2])


call_fortran(n, A)
