#!/usr/bin/env python
import ctypes


lib = ctypes.CDLL("./example4.so")

#example of fortran subroutine
#lib.sum_all_sub.restype = None
#lib.sum_all_sub.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), \
#                               ctypes.POINTER(ctypes.c_int)]


#example of fortran function
lib.sum_all_func.restype = ctypes.c_int
lib.sum_all_func.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]


if __name__=="__main__":
    array_length = 5
    input_arr_shape = ctypes.c_int * array_length
    input_arr = input_arr_shape(*range(array_length))
    input_arr_pointer = ctypes.cast(input_arr, ctypes.POINTER(ctypes.c_int))


    #example of fortran subroutine
    #resultd = ctypes.c_int()
    #lib.sum_all_sub(ctypes.byref(ctypes.c_int(array_length)), \
    #                input_arr_pointer, ctypes.byref(resultd))
    #print input_arr
    #print resultd.value

    #example of fortran function
    result_val = lib.sum_all_func(ctypes.byref(ctypes.c_int(array_length)), \
                                  input_arr_pointer)
    print result_val 
