#include <stdio.h>
#include <math.h>
#include "cuSZx_compress_float.h"

#include <cooperative_groups.h>

namespace cg = cooperative_groups;

//__device__
//void reduction(double sum1, double sum2,
//        double minDiff, double maxDiff, double sumDiff, double sumOfDiffSquare, 
//        double minErr, double maxErr, double sumErr, double sumErrSqr, double *results){
//
//    //static __shared__ double shared[10*10];
//    //dynamic shared mem
//    extern __shared__ float shared[];
//
//    int lane = threadIdx.x;
//    int wid = threadIdx.y;
//
//
//    for (int offset = warpSize/2; offset > 0; offset /= 2) 
//    {
//        minDiff = min(minDiff, __shfl_xor_sync(FULL_MASK, minDiff, offset));
//        maxDiff = max(maxDiff, __shfl_xor_sync(FULL_MASK, maxDiff, offset));
//        minErr = min(minErr, __shfl_xor_sync(FULL_MASK, minErr, offset));
//        maxErr = max(maxErr, __shfl_xor_sync(FULL_MASK, maxErr, offset));
//        sum1 += __shfl_down_sync(FULL_MASK, sum1, offset);
//        sum2 += __shfl_down_sync(FULL_MASK, sum2, offset);
//        sumDiff += __shfl_down_sync(FULL_MASK, sumDiff, offset);
//        sumOfDiffSquare += __shfl_down_sync(FULL_MASK, sumOfDiffSquare, offset);
//        sumErr += __shfl_down_sync(FULL_MASK, sumErr, offset);
//        sumErrSqr += __shfl_down_sync(FULL_MASK, sumErrSqr, offset);
//    }
//
//    if (lane==0){
//        shared[wid] = minDiff;
//        shared[blockDim.y+wid] = maxDiff;
//        shared[blockDim.y*2+wid] = minErr;
//        shared[blockDim.y*3+wid] = maxErr;
//        shared[blockDim.y*4+wid] = sum1;
//        shared[blockDim.y*5+wid] = sum2;
//        shared[blockDim.y*6+wid] = sumDiff;
//        shared[blockDim.y*7+wid] = sumOfDiffSquare;
//        shared[blockDim.y*8+wid] = sumErr;
//        shared[blockDim.y*9+wid] = sumErrSqr;
//    }
//
//    __syncthreads();                  
//
//    //if (wid==0)printf("ddata%i=%e:%e\n", 32*6+lane, shared[32*6+lane], ySum);
//
//    if (wid==0){
//        if (threadIdx.x < blockDim.y){
//            minDiff = shared[lane];
//            maxDiff = shared[blockDim.y+lane];
//            minErr = shared[blockDim.y*2+lane];
//            maxErr = shared[blockDim.y*3+lane];
//            sum1 = shared[blockDim.y*4+lane];
//            sum2 = shared[blockDim.y*5+lane];
//            sumDiff = shared[blockDim.y*6+lane];
//            sumOfDiffSquare = shared[blockDim.y*7+lane];
//            sumErr = shared[blockDim.y*8+lane];
//            sumErrSqr = shared[blockDim.y*9+lane];
//        }else{
//            minDiff = shared[0];  
//            maxDiff = shared[blockDim.y]; 
//            minErr = shared[blockDim.y*2]; 
//            maxErr = shared[blockDim.y*3]; 
//            sum1 = 0; 
//            sum2 = 0;
//            sumDiff = 0; 
//            sumOfDiffSquare = 0;
//            sumErr = 0;
//            sumErrSqr = 0;
//        }
//
//        for (int offset = warpSize/2; offset > 0; offset /= 2) 
//        {
//            minDiff = min(minDiff, __shfl_xor_sync(FULL_MASK, minDiff, offset));
//            maxDiff = max(maxDiff, __shfl_xor_sync(FULL_MASK, maxDiff, offset));
//            minErr = min(minErr, __shfl_xor_sync(FULL_MASK, minErr, offset));
//            maxErr = max(maxErr, __shfl_xor_sync(FULL_MASK, maxErr, offset));
//            sum1 += __shfl_down_sync(FULL_MASK, sum1, offset);
//            sum2 += __shfl_down_sync(FULL_MASK, sum2, offset);
//            sumDiff += __shfl_down_sync(FULL_MASK, sumDiff, offset);
//            sumOfDiffSquare += __shfl_down_sync(FULL_MASK, sumOfDiffSquare, offset);
//            sumErr += __shfl_down_sync(FULL_MASK, sumErr, offset);
//            sumErrSqr += __shfl_down_sync(FULL_MASK, sumErrSqr, offset);
//        }
//        
//        if (lane==0){
//            results[blockIdx.x] = minDiff;
//            results[gridDim.x+blockIdx.x] = minErr;
//            results[gridDim.x*2+blockIdx.x] = maxDiff;
//            results[gridDim.x*3+blockIdx.x] = maxErr;
//            results[gridDim.x*4+blockIdx.x] = sum1;
//            results[gridDim.x*5+blockIdx.x] = sum2;
//            results[gridDim.x*6+blockIdx.x] = sumDiff;
//            results[gridDim.x*7+blockIdx.x] = sumOfDiffSquare;
//            results[gridDim.x*8+blockIdx.x] = sumErr;
//            results[gridDim.x*9+blockIdx.x] = sumErrSqr;
//        }
//    }
//    //if(lane==0){
//    //    if (sum1>0.0)printf("test%i,%i,%i:%e\n",lane,wid,blockIdx.x, sum1);
//
//    //}
//
//}

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



//inline short getExponent_float(float value)
//{
//	//int ivalue = floatToBigEndianInt(value);
//
//	lfloat lbuf;
//	lbuf.value = value;
//	int ivalue = lbuf.ivalue;
//	
//	int expValue = (ivalue & 0x7F800000) >> 23;
//	expValue -= 127;
//	return (short)expValue;
//}
//
//inline void computeReqLength_float(double realPrecision, short radExpo, int* reqLength, float* medianValue)
//{
//	short reqExpo = getPrecisionReqLength_double(realPrecision);
//	*reqLength = 9+radExpo - reqExpo+1; //radExpo-reqExpo == reqMantiLength
//	if(*reqLength<9)
//		*reqLength = 9;
//	if(*reqLength>32)
//	{	
//		*reqLength = 32;
//		*medianValue = 0;
//	}			
//}

__device__ int _compute_reqLength(int redius, int absErrBound)
{
    int radExpo = (redius & 0x7F800000) >> 23;
    radExpo -= 127;
    int reqExpo = (absErrBound & 0x7F800000) >> 23;
    reqExpo -= 127;
    return 9+radExpo-reqExpo+1;
}

__device__ int _shfl_scan(int lznum, int *sums)
{
    // Below is the basic structure of using a shfl instruction
    // for a scan.
    // Record "value" as a variable - we accumulate it along the way
    int value = lznum;

    // Now accumulate in log steps up the chain
    // compute sums, with another thread's value who is
    // distance delta away (i).  Note
    // those threads where the thread 'i' away would have
    // been out of bounds of the warp are unaffected.  This
    // creates the scan sum.

#pragma unroll
    for (int i = 1; i <= warpSize; i *= 2) {
        unsigned int mask = 0xffffffff;
        int n = __shfl_up_sync(mask, value, i);

        if (threadIdx.x >= i) value += n;
                      
    }

    // value now holds the scan value for the individual thread
    // next sum the largest values for each warp

    // write the sum of the warp to smem
    if (threadIdx.x == warpSize - 1) {
        sums[threadIdx.y] = value;
    }
    __syncthreads();

    //
    // scan sum the warp sums
    // the same shfl scan operation, but performed on warp sums
    //
    if (threadIdx.y == 0 && threadIdx.x < blockDim.y) {
        int warp_sum = sums[threadIdx.x];

        int mask = (1 << blockDim.y) - 1;
        for (int i = 1; i <= blockDim.y; i *= 2) {
            //int n = __shfl_up_sync(mask, warp_sum, i, blockDim.y);
            int n = __shfl_up_sync(mask, warp_sum, i);
            if (threadIdx.x >= i) warp_sum += n;
        }

        sums[threadIdx.x] = warp_sum;
    }
    __syncthreads();

    // perform a uniform add across warps in the block
    // read neighbouring warp's sum and add it to threads value
    int blockSum = 0;
    if (threadIdx.y > 0) {
        blockSum = sums[threadIdx.y - 1];
    }
    value += blockSum;

    return value;
}

__device__ int _compute_oneBlock(int base, int reqLength, float *value, int *ivalue, uchar4 *cvalue, int *sums, unsigned char *meta, unsigned char *midBytes, bool bi)
{
	int reqBytesLength;
	int rightShiftBits;


	if (reqLength%8 != 0)
	{
		reqBytesLength = reqLength/8+1;		
		rightShiftBits = 8 - reqLength%8;
    }else{
		reqBytesLength = reqLength/8;		
		rightShiftBits = 0;
    }

    //if (bi==true) printf("%i:%i:%i:%u\n", blockIdx.x, threadIdx.x, threadIdx.y, cvalue[threadIdx.y*blockDim.x+threadIdx.x].w);
    int cur_ivalue = (ivalue[threadIdx.y*blockDim.x+threadIdx.x] >> rightShiftBits) & ((1<<(32-rightShiftBits))-1);
    ivalue[threadIdx.y*blockDim.x+threadIdx.x] = cur_ivalue;
    __syncthreads();                  

    int pre_ivalue = 0;
    if (threadIdx.x!=0 || threadIdx.y!=0) pre_ivalue = ivalue[threadIdx.y*blockDim.x+threadIdx.x-1];
    pre_ivalue = cur_ivalue ^ pre_ivalue;
    __syncthreads();                  

    unsigned char leadingNum = 0;
    if (reqBytesLength == 2)
    {
        if (pre_ivalue >> 16 == 0) leadingNum = 2;
        else if (pre_ivalue >> 24 == 0) leadingNum = 1;
    }else if (reqBytesLength == 3)
    {
        if (pre_ivalue >> 8 == 0) leadingNum = 3;
        else if (pre_ivalue >> 16 == 0) leadingNum = 2;
        else if (pre_ivalue >> 24 == 0) leadingNum = 1;
    }else if (reqBytesLength == 1)
    {
        if (pre_ivalue >> 24 == 0) leadingNum = 1;

    }else if (reqBytesLength == 4)
    {
        if (pre_ivalue == 0) leadingNum = 4;
        else if (pre_ivalue >> 8 == 0) leadingNum = 3;
        else if (pre_ivalue >> 16 == 0) leadingNum = 2;
        else if (pre_ivalue >> 24 == 0) leadingNum = 1;
    }
    //if (bi==true) printf("%i:%i:%i:%u\n", blockIdx.x, threadIdx.x, threadIdx.y, leadingNum);
    //midBytes[base+threadIdx.y*blockDim.x+threadIdx.x] = leadingNum; 

    int midByte_size = reqBytesLength - (int)leadingNum;
    int midByte_sum = _shfl_scan(midByte_size, sums);
    uchar4 cur_cvalue = cvalue[threadIdx.y*blockDim.x+threadIdx.x];
    if (reqBytesLength == 2)
    {
        if (midByte_size == 1){
            midBytes[base+midByte_sum-1] = cur_cvalue.z; 
            if (bi==true) printf("%i:%i:%i:%u\n", blockIdx.x, threadIdx.x, threadIdx.y, cur_cvalue.z);
        }else if (midByte_size == 2){
            midBytes[base+midByte_sum-1] = cur_cvalue.w; 
            if (bi==true) printf("%i:%i:%i:%u\n", blockIdx.x, threadIdx.x, threadIdx.y, cur_cvalue.z);
            midBytes[base+midByte_sum-2] = cur_cvalue.z;
            if (bi==true) printf("%i:%i:%i:%u\n", blockIdx.x, threadIdx.x, threadIdx.y, cur_cvalue.w);
        }
    }else if (reqBytesLength == 3)
    {
        if (midByte_size == 1){
            midBytes[base+midByte_sum-1] = cur_cvalue.y; 
        }else if (midByte_size == 2){
            midBytes[base+midByte_sum-1] = cur_cvalue.z; 
            midBytes[base+midByte_sum-2] = cur_cvalue.y; 
        }else if (midByte_size == 3){
            midBytes[base+midByte_sum-1] = cur_cvalue.w; 
            midBytes[base+midByte_sum-2] = cur_cvalue.z; 
            midBytes[base+midByte_sum-3] = cur_cvalue.y; 
        }
    }else if (reqBytesLength == 1)
    {
        if (midByte_size == 1)
            midBytes[base+midByte_sum-1] = cur_cvalue.w; 
    }else if (reqBytesLength == 4)
    {
        if (midByte_size == 1){
            midBytes[base+midByte_sum-1] = cur_cvalue.x; 
        }else if (midByte_size == 2){
            midBytes[base+midByte_sum-1] = cur_cvalue.y; 
            midBytes[base+midByte_sum-2] = cur_cvalue.x; 
        }else if (midByte_size == 3){
            midBytes[base+midByte_sum-1] = cur_cvalue.z; 
            midBytes[base+midByte_sum-2] = cur_cvalue.y; 
            midBytes[base+midByte_sum-3] = cur_cvalue.x; 
        }else if (midByte_size == 4){
            midBytes[base+midByte_sum-1] = cur_cvalue.w; 
            midBytes[base+midByte_sum-2] = cur_cvalue.z; 
            midBytes[base+midByte_sum-3] = cur_cvalue.y; 
            midBytes[base+midByte_sum-4] = cur_cvalue.x; 
        }
    }
}

__global__ void compress_float(float *oriData, unsigned char *meta, unsigned char *midBytes, float absErrBound, int bs, size_t nb, size_t mSize, int *test) 
{
    int tidx = threadIdx.x;
    int tidy = threadIdx.y;
    int bid = blockIdx.x;

    float data, radius, medianValue;
    unsigned mask;
    unsigned char state;
    extern __shared__ float shared[];
    float* value = shared;
    int* ivalue = (int*)shared;
    uchar4* cvalue = (uchar4*)shared;
    int* sums = &ivalue[bs];


    bool bi = false;
    if (bid==73) bi=true;
    for (int b=bid; b<nb; b+=gridDim.x){
        data = oriData[b*bs+tidy*warpSize+tidx];
        float Min = data;
        float Max = data;

        for (int offset = warpSize/2; offset > 0; offset /= 2) 
        {
            Min = min(Min, __shfl_xor_sync(FULL_MASK, Min, offset));
            Max = max(Max, __shfl_xor_sync(FULL_MASK, Max, offset));
        }
        if (tidx==0){
            value[tidy] = Min;
            value[blockDim.y+tidy] = Max;
        }
        __syncthreads();                  

        if (tidy==0){
            if (tidx < blockDim.y){
                Min = value[tidx];
                Max = value[blockDim.y+tidx];
            }

            mask = __ballot_sync(FULL_MASK, tidx < blockDim.y);
            for (int offset = blockDim.y/2; offset > 0; offset /= 2) 
            {
                Min = min(Min, __shfl_xor_sync(mask, Min, offset));
                Max = max(Max, __shfl_xor_sync(mask, Max, offset));
            }
            
            if (tidx==0){
                radius = (Max - Min)/2;
                value[0] = radius;
                value[1] = Min + radius;
                value[2] = absErrBound;
            }
        }
        __syncthreads();                  

        radius = value[0];
        medianValue = value[1];
        state = radius <= absErrBound ? 0 : 1;
        if (tidx==0) meta[b*mSize] = cvalue[0].x;
        //if (tidx==0) test[b] = ivalue[0];
        __syncthreads();                  

        if (state==1){
            int reqLength = _compute_reqLength(ivalue[0], ivalue[2]);
            if (tidx==0) test[b] = reqLength;
            __syncthreads();                  
            value[tidy*blockDim.x+tidx] = data - medianValue;
            __syncthreads();                  
            _compute_oneBlock(b*bs*sizeof(float), reqLength, value, ivalue, cvalue, sums, meta, midBytes, bi);
            bi = false;
        }

    }

}
