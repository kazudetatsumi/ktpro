#!/usr/bin/env python
import ctypes

lib = ctypes.CDLL("./example3.so")

class day_type(ctypes.Structure):
    _fields_ = [("month", ctypes.c_int), ("date", ctypes.c_int)]
input_day = day_type(month=6, date=20)

input_val = ctypes.c_int(3)

array_length = 10
input_arr_shape = ctypes.c_int * array_length
input_arr = input_arr_shape(*range(array_length))

lib.input.restype = None
lib.input.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), \
                    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(day_type)]

if __name__=="__main__":

    print ("--- input ---")
    input_arr_pointer = ctypes.cast(input_arr, ctypes.POINTER(ctypes.c_int))
    lib.input(ctypes.byref(input_val), ctypes.byref(ctypes.c_int(array_length)), \
            input_arr_pointer, ctypes.byref(input_day))
