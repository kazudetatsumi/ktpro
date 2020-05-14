#!/usr/bin/env python

import ctypes

lib = ctypes.CDLL("./example2.so")

# scalar
int_val = ctypes.c_int(3)

# array
array_length = 10
int_arr_shape = ctypes.c_int * array_length
int_arr = int_arr_shape(*range(array_length))

# prottype
lib.add.restype = ctypes.c_int
lib.add.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]


# cast
#print ("array object type: ", type(int_arr))
print "array object type: ", type(int_arr)
for item in int_arr:
    print item
#for item in int_arr:
#    print (item, end="  ")
#else:
#    print ()
