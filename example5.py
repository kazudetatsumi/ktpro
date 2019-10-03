import numpy as np
import ctypes

lib = ctypes.CDLL("./example5.so")

#class result(ctypes.Structure):
#    _fields_ =[("len", ctypes.c_int), ("arr", ctypes.POINTER(ctypes.c_int))]
class result(ctypes.Structure):
    _fields_ =[("len0", ctypes.c_int), ("len1", ctypes.c_int), ("arr", ctypes.POINTER(ctypes.c_double))]

# declaration of prototypes using ndpointer
#lib.extract_plus.restype = result
#lib.extract_plus.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.int32, ndim=1)]
lib.plus1.restype = result
lib.plus1.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=2)]

lib.delete_array.restype = None
#lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.int32, ndim=1)]
lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=2)]

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
#A_length0 = 10
#A_length1 = 20
#A = np.ones((A_length0, A_length1), dtype=np.float64)*3**(0.5)
A = np.array([[1.0, 2.0, 3.0, 4.0],[10.0, 20.0, 30.0, 40.0], [100.0, 200.0, 300.0, 400.0]])
print "A.shape[0]",A.shape[0]
print "A.shape[1]",A.shape[1]
print A
result = lib.plus1(ctypes.byref(ctypes.c_int(A.shape[0])), ctypes.byref(ctypes.c_int(A.shape[1])), A)
B_len0 =  result.len0
B_len1 =  result.len1
print B_len0
print B_len1
B = np.ctypeslib.as_array(result.arr, shape=(B_len0, B_len1))

print B.shape[0]
print B.shape[1]
print B

lib.delete_array(ctypes.byref(ctypes.c_int(B_len0)), ctypes.byref(ctypes.c_int(B_len1)), B)
