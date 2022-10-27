## This is a very basic demo of using SZx in Python
from ctypes import *
import os
import ctypes
import numpy as np
import random

sotest = cdll.LoadLibrary("../build/lib/libSZx.so")

# sotest.SZ_fast_compress_args.argtypes = [c_float, c_float]
sotest.SZ_fast_compress_args.restype = c_char_p

outSize = c_size_t()

INPUT = c_float * 1000
data = INPUT()
for i in range(1000):
    data[i] = random.random()


addr = sotest.SZ_fast_compress_args(1, 0, data, byref(outSize),
                                    0, c_float(0.8), c_float(0),
                                    c_size_t(0), c_size_t(0), c_size_t(0), c_size_t(0), c_size_t(1000))

print(type(addr))
print(addr)

print(outSize)
