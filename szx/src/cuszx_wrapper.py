import numpy as np
import ctypes
from ctypes import *
import random
import cupy as cp

LIB_PATH = './cuszx_wrapper.so'

def get_cuszx_compress():
    dll = ctypes.CDLL(LIB_PATH, mode=ctypes.RTLD_GLOBAL)
    func = dll.cuSZx_integrated_compress
    # unsigned char *bytes, float *data, float r2r_threshold, float r2r_err, size_t nbEle, int blockSize, size_t *outSize
    func.argtypes = [POINTER(c_ubyte), POINTER(c_float), c_float, c_float, c_size_t, c_int, POINTER(c_size_t)]
    return func

def get_cuszx_decompress():
    dll = ctypes.CDLL(LIB_PATH, mode=ctypes.RTLD_GLOBAL)
    func = dll.cuSZx_integrated_decompress
    # unsigned char *bytes, float *data, float r2r_threshold, float r2r_err, size_t nbEle, int blockSize, size_t *outSize
    func.argtypes = [POINTER(c_float), POINTER(c_char), c_size_t]
    return func

def get_device_compress():
    dll = ctypes.CDLL(LIB_PATH, mode=ctypes.RTLD_GLOBAL)
    func = dll.cuSZx_device_compress
    # Returns: unsigned char *bytes
    # Needs: float *oriData, size_t *outSize, float absErrBound, size_t nbEle, int blockSize, float threshold
    func.argtypes = [POINTER(c_float), POINTER(c_size_t), c_float, c_size_t, c_int, c_float]
    func.restype = POINTER(c_ubyte)
    return func

def get_device_decompress():
    dll = ctypes.CDLL(LIB_PATH, mode=ctypes.RTLD_GLOBAL)
    func = dll.cuSZx_device_decompress
    # Returns: unsigned char *bytes
    # Needs: float *oriData, size_t *outSize, float absErrBound, size_t nbEle, int blockSize, float threshold
    func.argtypes = [POINTER(POINTER(c_float)), c_size_t, POINTER(c_ubyte)]
    return func

__cuszx_compress = get_cuszx_compress()
__cuszx_device_compress = get_device_compress()
__cuszx_device_decompress=get_device_decompress()
#__cuszx_decompress=get_cuszx_decompress()

def cuszx_integrated_compress(bytes_a, data, r2r_threshold, r2r_error, nbEle, blockSize, outSize):
    #bytes_p = bytes_a.ctypes.data_as(POINTER(c_ubyte))
    bytes_p=bytes_a
    data_p = data.ctypes.data_as(POINTER(c_float))
    #r2r_threshold_p = r2r_threshold.ctypes.data_as(c_float)
    #r2r_error_p = r2r_error.ctypes.data_as(c_float)
    #nbEle_p = nbEle.ctypes.data_as(c_size_t)
    #blockSize_p = blockSize.ctypes.data_as(c_int)
    #outSize_p = outSize.ctypes.data_as(POINTER(c_size_t))
    #print(type(bytes_p))
    print(type(outSize))
    __cuszx_compress(bytes_p, data_p, r2r_threshold, r2r_error, nbEle, blockSize, outSize)

    return bytes_p, outSize

def cuszx_device_compress(oriData, outSize, absErrBound, nbEle, blockSize,threshold):
    oriData_p = ctypes.cast(oriData.data.ptr, ctypes.POINTER(c_float))
    
    o_bytes = __cuszx_device_compress(oriData_p, outSize, absErrBound, nbEle, blockSize, threshold)
    print(o_bytes)
    return o_bytes, outSize

def cuszx_integrated_decompress(data,bytes,nbEle):
    data_p = data.ctypes.data_as(POINTER(c_float))
    bytes_p = bytes.ctypes.data_as(POINTER(c_char))
    nbEle_p = nbEle.ctypes.data_as(c_size_t)
    __cuszx_decompress(data_p, bytes_p, nbEle_p)

    return data

def cuszx_device_decompress(nbEle, cmpBytes):
    
    newData_p = ctypes.c_float()
    newData = ctypes.pointer(ctypes.pointer(newData_p))
    nbEle_p = ctypes.c_size_t(nbEle)
    __cuszx_device_decompress(newData,nbEle_p,cmpBytes)
    return newData

if __name__ == "__main__":
    
    DATA_SIZE = 1024
    MAX_D = 10.0
    MIN_D = -10.0
    RANGE = MAX_D - MIN_D
    r2r_threshold = 0.1
    r2r_error = 0.1

    in_vector = np.zeros((DATA_SIZE,))
    for i in range(0,int(DATA_SIZE/4)):
        in_vector[i] = 0.0
    for i in range(int(DATA_SIZE/4), int(2*DATA_SIZE/4)):
        in_vector[i] = 5.0
    for i in range(int(2*DATA_SIZE/4), int(3*DATA_SIZE/4)):
        in_vector[i] = random.uniform(MIN_D, MAX_D)
    for i in range(int(3*DATA_SIZE/4), int(3*DATA_SIZE/4)+6):
        in_vector[i] = -7.0
    for i in range(int(3*DATA_SIZE/4)+6, DATA_SIZE):
        in_vector[i] = 0.001

    #a = b[0,:,:]

    in_vector = in_vector.astype('float32')
    in_vector_gpu = cp.asarray(in_vector)
    outbytes = POINTER(c_ubyte)()
    variable = ctypes.c_size_t(10)
    outSize = ctypes.pointer(variable)
    #print(outSize.contents.value)
    # print(outbytes.contents)
    # outbytes, f_size = cuszx_integrated_compress(outbytes, in_vector, np.float32(r2r_threshold), np.float32(r2r_error), np.ulonglong(DATA_SIZE), np.int32(256), outSize)
    o_bytes, outSize = cuszx_device_compress(in_vector_gpu, outSize, np.float32(r2r_error), np.ulonglong(DATA_SIZE), np.int32(256),np.float32(r2r_threshold))
    #print(outSize.contents.value)
    #print(o_bytes.contents.value)
    print("Compress Success...starting decompress ")
    d_bytes = cuszx_device_decompress(DATA_SIZE, o_bytes)
    print("Decompress Success")
