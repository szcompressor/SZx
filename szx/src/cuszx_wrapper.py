import numpy as np
import ctypes
from ctypes import *

LIB_PATH = './cuszx_wrapper.so'

def get_cuszx_compress():
    dll = ctypes.CDLL(LIB_PATH, mode=ctypes.RTLD_GLOBAL)
    func = dll.cuSZx_integrated_compress
    # unsigned char *bytes, float *data, float r2r_threshold, float r2r_err, size_t nbEle, int blockSize, size_t *outSize
    func.argtypes = [POINTER(c_char), POINTER(c_float), c_float, c_float, c_size_t, c_int, POINTER(c_size_t)]
    return func

def get_cuszx_decompress():
    dll = ctypes.CDLL(LIB_PATH, mode=ctypes.RTLD_GLOBAL)
    func = dll.cuSZx_integrated_decompress
    # unsigned char *bytes, float *data, float r2r_threshold, float r2r_err, size_t nbEle, int blockSize, size_t *outSize
    func.argtypes = [POINTER(c_float), POINTER(c_char), c_size_t]
    return func

__cuszx_compress = get_cuszx_compress()
__cuszx_decompress=get_cuszx_decompress()

def cuszx_integrated_compress(bytes, data, r2r_threshold, r2r_error, nbEle, blockSize, outSize):
    bytes_p = bytes.ctypes.data_as(POINTER(c_char))
    data_p = data.ctypes.data_as(POINTER(c_float))
    r2r_threshold_p = r2r_threshold.ctypes.data_as(c_float)
    r2r_error_p = r2r_error.ctypes.data_as(c_float)
    nbEle_p = nbEle.ctypes.data_as(c_size_t)
    blockSize_p = blockSize.ctypes.data_as(c_int)
    outSize_p = outSize.ctypes.data_as(POINTER(c_size_t))

    __cuszx_compress(bytes_p, data_p, r2r_threshold_p, r2r_error_p, nbEle_p, blockSize_p, outSize_p)

    return bytes, outSize

def cuszx_integrated_decompress(data,bytes,nbEle):
    data_p = data.ctypes.data_as(POINTER(c_float))
    bytes_p = bytes.ctypes.data_as(POINTER(c_char))
    nbEle_p = nbEle.ctypes.data_as(c_size_t)
    __cuszx_decompress(data_p, bytes_p, nbEle_p)

    return data