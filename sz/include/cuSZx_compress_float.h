#ifndef CUSZX_COMPRESS_FLOAT_H
#define CUSZX_COMPRESS_FLOAT_H

#include <cuda_runtime.h>

// Utilities and system includes
#include <helper_cuda.h>  // helper function CUDA error checking and initialization
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples

#define FULL_MASK 0xffffffff

__device__
void reduction(double sum1, double sum2,
        double minDiff, double maxDiff, double sumDiff, double sumOfDiffSquare, 
        double minErr, double maxErr, double sumErr, double sumErrSqr);

__global__ void compress_float(float *oriData, unsigned char *meta, unsigned char *midBytes, float absErrBound, int bs, size_t nb, size_t mSize, int *test); 

#endif /* ----- #ifndef CUSZX_COMPRESS_FLOAT_H  ----- */
