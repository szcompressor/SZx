import numpy as np
import ctypes
from ctypes import *
import random
import cupy as cp

LIB_PATH = './cuszx_wrapper.so'


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
    # Returns: float *newData
    # Needs: size_t nbEle, unsigned char *cmpBytes
    func.argtypes = [c_size_t, POINTER(c_ubyte)]
    func.restype = POINTER(c_float)
    return func

__cuszx_device_compress = get_device_compress()
__cuszx_device_decompress=get_device_decompress()


def cuszx_device_compress(oriData, outSize, absErrBound, nbEle, blockSize,threshold):
    oriData_p = ctypes.cast(oriData.data.ptr, ctypes.POINTER(c_float))
    
    o_bytes = __cuszx_device_compress(oriData_p, outSize,np.float32(absErrBound), np.ulonglong(DATA_SIZE), np.int32(blockSize),np.float32(threshold))
    print(o_bytes)
    return o_bytes, outSize


def cuszx_device_decompress(nbEle, cmpBytes):
    
    nbEle_p = ctypes.c_size_t(nbEle)
    newData = __cuszx_device_decompress(nbEle_p,cmpBytes)
    return newData

### Example of device compress/decompress wrapper usage

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


    in_vector = in_vector.astype('float32')
    in_vector_gpu = cp.asarray(in_vector)

    variable = ctypes.c_size_t(0)
    outSize = ctypes.pointer(variable)

    o_bytes, outSize = cuszx_device_compress(in_vector_gpu, outSize, r2r_error, DATA_SIZE, 256, r2r_threshold)
    print("Compress Success...starting decompress ")
    d_bytes = cuszx_device_decompress(DATA_SIZE, o_bytes)
    print("Decompress Success")
