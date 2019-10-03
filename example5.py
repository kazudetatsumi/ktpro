import numpy as np
import ctypes

lib = ctypes.CDLL("./example5.so")

#class result(ctypes.Structure):
#    _fields_ =[("len", ctypes.c_int), ("arr", ctypes.POINTER(ctypes.c_int))]
class result(ctypes.Structure):
    _fields_ =[("len", ctypes.c_int), ("arr", ctypes.POINTER(ctypes.c_double))]

# declaration of prototypes using ndpointer
#lib.extract_plus.restype = result
#lib.extract_plus.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.int32, ndim=1)]
lib.plus1.restype = result
lib.plus1.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]

lib.delete_array.restype = None
#lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.int32, ndim=1)]
lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]

#if __name__ == "__main__":
#array_length = 10
#arr_np = np.random.randint(low=-100, high=100, size=array_length, dtype=np.int32)
#result = lib.extract_plus(ctypes.byref(ctypes.c_int(array_length)), arr_np)
#result = lib.extract_plus(ctypes.byref(ctypes.c_int(array_length)), arr_np)
#result_len = result.len
#result_vec = np.ctypeslib.as_array(result.arr, shape=(result_len,))

#for i in result_vec:
#    print i

#lib.delete_array(ctypes.byref(ctypes.c_int(result_len)), result_vec)
A_length = 10
A = np.ones((A_length), dtype=np.float64)*3**(0.5)
print A.dtype
result = lib.plus1(ctypes.byref(ctypes.c_int(A_length)), A)
B_len =  result.len
B = np.ctypeslib.as_array(result.arr, shape=(B_len,))
print B.dtype

for i in B:
    print i

lib.delete_array(ctypes.byref(ctypes.c_int(B_len)), B)
