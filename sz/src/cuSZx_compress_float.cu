#include <stdio.h>
#include <math.h>
#include "cuSZx_compress_float.h"

#include <cooperative_groups.h>

namespace cg = cooperative_groups;

__device__
void reduction(double sum1, double sum2,
        double minDiff, double maxDiff, double sumDiff, double sumOfDiffSquare, 
        double minErr, double maxErr, double sumErr, double sumErrSqr, double *results){

    //static __shared__ double shared[10*10];
    //dynamic shared mem
    extern __shared__ float shared[];

    int lane = threadIdx.x;
    int wid = threadIdx.y;


    for (int offset = warpSize/2; offset > 0; offset /= 2) 
    {
        minDiff = min(minDiff, __shfl_xor_sync(FULL_MASK, minDiff, offset));
        maxDiff = max(maxDiff, __shfl_xor_sync(FULL_MASK, maxDiff, offset));
        minErr = min(minErr, __shfl_xor_sync(FULL_MASK, minErr, offset));
        maxErr = max(maxErr, __shfl_xor_sync(FULL_MASK, maxErr, offset));
        sum1 += __shfl_down_sync(FULL_MASK, sum1, offset);
        sum2 += __shfl_down_sync(FULL_MASK, sum2, offset);
        sumDiff += __shfl_down_sync(FULL_MASK, sumDiff, offset);
        sumOfDiffSquare += __shfl_down_sync(FULL_MASK, sumOfDiffSquare, offset);
        sumErr += __shfl_down_sync(FULL_MASK, sumErr, offset);
        sumErrSqr += __shfl_down_sync(FULL_MASK, sumErrSqr, offset);
    }

    if (lane==0){
        shared[wid] = minDiff;
        shared[blockDim.y+wid] = maxDiff;
        shared[blockDim.y*2+wid] = minErr;
        shared[blockDim.y*3+wid] = maxErr;
        shared[blockDim.y*4+wid] = sum1;
        shared[blockDim.y*5+wid] = sum2;
        shared[blockDim.y*6+wid] = sumDiff;
        shared[blockDim.y*7+wid] = sumOfDiffSquare;
        shared[blockDim.y*8+wid] = sumErr;
        shared[blockDim.y*9+wid] = sumErrSqr;
    }

    __syncthreads();                  

    //if (wid==0)printf("ddata%i=%e:%e\n", 32*6+lane, shared[32*6+lane], ySum);

    if (wid==0){
        if (threadIdx.x < blockDim.y){
            minDiff = shared[lane];
            maxDiff = shared[blockDim.y+lane];
            minErr = shared[blockDim.y*2+lane];
            maxErr = shared[blockDim.y*3+lane];
            sum1 = shared[blockDim.y*4+lane];
            sum2 = shared[blockDim.y*5+lane];
            sumDiff = shared[blockDim.y*6+lane];
            sumOfDiffSquare = shared[blockDim.y*7+lane];
            sumErr = shared[blockDim.y*8+lane];
            sumErrSqr = shared[blockDim.y*9+lane];
        }else{
            minDiff = shared[0];  
            maxDiff = shared[blockDim.y]; 
            minErr = shared[blockDim.y*2]; 
            maxErr = shared[blockDim.y*3]; 
            sum1 = 0; 
            sum2 = 0;
            sumDiff = 0; 
            sumOfDiffSquare = 0;
            sumErr = 0;
            sumErrSqr = 0;
        }

        for (int offset = warpSize/2; offset > 0; offset /= 2) 
        {
            minDiff = min(minDiff, __shfl_xor_sync(FULL_MASK, minDiff, offset));
            maxDiff = max(maxDiff, __shfl_xor_sync(FULL_MASK, maxDiff, offset));
            minErr = min(minErr, __shfl_xor_sync(FULL_MASK, minErr, offset));
            maxErr = max(maxErr, __shfl_xor_sync(FULL_MASK, maxErr, offset));
            sum1 += __shfl_down_sync(FULL_MASK, sum1, offset);
            sum2 += __shfl_down_sync(FULL_MASK, sum2, offset);
            sumDiff += __shfl_down_sync(FULL_MASK, sumDiff, offset);
            sumOfDiffSquare += __shfl_down_sync(FULL_MASK, sumOfDiffSquare, offset);
            sumErr += __shfl_down_sync(FULL_MASK, sumErr, offset);
            sumErrSqr += __shfl_down_sync(FULL_MASK, sumErrSqr, offset);
        }
        
        if (lane==0){
            results[blockIdx.x] = minDiff;
            results[gridDim.x+blockIdx.x] = minErr;
            results[gridDim.x*2+blockIdx.x] = maxDiff;
            results[gridDim.x*3+blockIdx.x] = maxErr;
            results[gridDim.x*4+blockIdx.x] = sum1;
            results[gridDim.x*5+blockIdx.x] = sum2;
            results[gridDim.x*6+blockIdx.x] = sumDiff;
            results[gridDim.x*7+blockIdx.x] = sumOfDiffSquare;
            results[gridDim.x*8+blockIdx.x] = sumErr;
            results[gridDim.x*9+blockIdx.x] = sumErrSqr;
        }
    }
    //if(lane==0){
    //    if (sum1>0.0)printf("test%i,%i,%i:%e\n",lane,wid,blockIdx.x, sum1);

    //}

}

__device__
void gridReduction_cg(double *results) 
{
    int tidx = threadIdx.x;
    int tidy = threadIdx.y;
    int bid = blockIdx.x;

    if (bid==0){
        double data = results[tidy*gridDim.x+tidx];

        for (int i=(tidx+blockDim.x); i<gridDim.x; i+=blockDim.x){
            if (tidy<2) data = min(data, results[tidy*gridDim.x+i]);
            else if (tidy<4) data = max(data, results[tidy*gridDim.x+i]);
            else data += results[tidy*gridDim.x+i];
        }
        __syncthreads();                  

        for (int offset = warpSize/2; offset > 0; offset /= 2) 
        {
            if (tidy<2) data = min(data, __shfl_xor_sync(FULL_MASK, data, offset));
            else if (tidy<4) data = max(data, __shfl_xor_sync(FULL_MASK, data, offset));
            else data += __shfl_down_sync(FULL_MASK, data, offset);
        }

        if (tidx==0) results[tidy*gridDim.x] = data;
    }
}

__global__ void compress_float(float *oriData, unsigned char *meta, unsigned char *midBytes, float absErrBound, int bs, size_t nb, size_t mSize) 
{
    int tidx = threadIdx.x;
    int tidy = threadIdx.y;
    int bid = blockIdx.x;

    float data, radius, medianValue;
    unsigned mask;
    unsigned char state;
    extern __shared__ float shared[];

    for (int b=bid; b<nb; b+=gridDim.x){
        data = oriData[b*bs+tidy*warpSize+tidx];
        float Min = data;
        float Max = data;

        int lane = threadIdx.x;
        int wid = threadIdx.y;

        for (int offset = warpSize/2; offset > 0; offset /= 2) 
        {
            Min = min(Min, __shfl_xor_sync(FULL_MASK, Min, offset));
            Max = max(Max, __shfl_xor_sync(FULL_MASK, Max, offset));
        }
        if (tidx==0){
            shared[tidy] = Min;
            shared[blockDim.y+tidy] = Max;
        }
        __syncthreads();                  

        if (tidy==0){
            if (tidx < blockDim.y){
                Min = shared[tidx];
                Max = shared[blockDim.y+tidx];
            }
           // else{
           //     minDiff = shared[0];  
           //     maxDiff = shared[blockDim.y]; 
           //     minErr = shared[blockDim.y*2]; 
           //     maxErr = shared[blockDim.y*3]; 
           //     sum1 = 0; 
           //     sum2 = 0;
           //     sumDiff = 0; 
           //     sumOfDiffSquare = 0;
           //     sumErr = 0;
           //     sumErrSqr = 0;
           // }

            mask = __ballot_sync(FULL_MASK, tidx < blockDim.y);
            for (int offset = blockDim.y/2; offset > 0; offset /= 2) 
            {
                Min = min(Min, __shfl_xor_sync(mask, Min, offset));
                Max = max(Max, __shfl_xor_sync(mask, Max, offset));
            }
            
            if (tidx==0){
                shared[0] = (Max - Min)/2;
                shared[1] = Min + radius;
            }
        }
        __syncthreads();                  

        radius = shared[0];
        medianValue = shared[1];
        state = radius <= absErrBound ? 0 : 1;
        if (tidx==0) meta[b*mSize] = state;

    }

}
