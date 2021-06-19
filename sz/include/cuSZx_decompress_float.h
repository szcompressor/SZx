#ifndef CUSZX_DECOMPRESS_FLOAT_H
#define CUSZX_DECOMPRESS_FLOAT_H

#include <cuda_runtime.h>

// Utilities and system includes
#include <helper_cuda.h>  // helper function CUDA error checking and initialization
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples

__global__ void decompress_float(unsigned char *data, int bs, size_t nc, size_t mSize); 

#endif /* ----- #ifndef CUSZX_DECOMPRESS_FLOAT_H  ----- */
