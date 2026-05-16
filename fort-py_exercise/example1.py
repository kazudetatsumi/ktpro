#!/usr/bin/env python
import ctypes

lib = ctypes.CDLL("./example1.so")

lib.print_hello()

