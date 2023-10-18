#ifndef CUSZX_FLOAT_H
#define CUSZX_FLOAT_H

namespace szx{

__device__
void reduction(double sum1, double sum2,
        double minDiff, double maxDiff, double sumDiff, double sumOfDiffSquare, 
        double minErr, double maxErr, double sumErr, double sumErrSqr);

__global__ void compress_float(float *oriData, unsigned char *meta, short *offsets, unsigned char *midBytes, float absErrBound, int bs, size_t nb, size_t mSize); 

}
#endif /* ----- #ifndef CUSZX_COMPRESS_FLOAT_H  ----- */
