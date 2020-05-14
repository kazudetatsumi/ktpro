import ctypes
import os

lib = ctypes.CDLL("./example4.so")


lib.sum_all_func.restype = ctypes.c_int
lib.sum_all_func.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]

if __name__=="__main__":
    offset = 3
    array_length = 10
    input_arr_shape = ctypes.c_int * array_length
    input_arr = input_arr_shape(*range(array_length))
    input_arr_pointer = ctypes.cast(input_arr, ctypes.POINTER(ctypes.c_int))

    result_arr = input_arr_shape()
    result_arr_pointer = ctypes.cast(result_arr, ctypes.POINTER(ctypes.c_int))


    print ("--- sum_all_func ---")
    result_val = lib.sum_all_func(ctypes.byref(ctypes.c_int(array_length)), \
                                  input_arr_pointer)
    print result_val


