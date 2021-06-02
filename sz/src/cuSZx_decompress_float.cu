#include <stdio.h>
#include <math.h>
#include "cuSZx_decompress_float.h"

#include <cooperative_groups.h>

namespace cg = cooperative_groups;

__device__ int _decompress_oneBlock(int bbase, int mbase, int obase, int reqLength, float *value, int *ivalue, uchar4 *cvalue, int *sums, unsigned char *meta, short *offsets, unsigned char *midBytes, bool bi)
{
	//int reqBytesLength;
	//int rightShiftBits;


	//if (reqLength%8 != 0)
	//{
	//	reqBytesLength = reqLength/8+1;		
	//	rightShiftBits = 8 - reqLength%8;
    //}else{
	//	reqBytesLength = reqLength/8;		
	//	rightShiftBits = 0;
    //}

    ////if (bi==true) printf("%i:%i:%i:%u\n", blockIdx.x, threadIdx.x, threadIdx.y, cvalue[threadIdx.y*blockDim.x+threadIdx.x].w);
    //int cur_ivalue = (ivalue[threadIdx.y*blockDim.x+threadIdx.x] >> rightShiftBits) & ((1<<(32-rightShiftBits))-1);
    //ivalue[threadIdx.y*blockDim.x+threadIdx.x] = cur_ivalue;
    //__syncthreads();                  

    //int pre_ivalue = 0;
    //if (threadIdx.x!=0 || threadIdx.y!=0) pre_ivalue = ivalue[threadIdx.y*blockDim.x+threadIdx.x-1];
    //pre_ivalue = cur_ivalue ^ pre_ivalue;
    //__syncthreads();                  

    //int leadingNum = 0;
    //if (reqBytesLength == 2)
    //{
    //    if (pre_ivalue >> 16 == 0) leadingNum = 2;
    //    else if (pre_ivalue >> 24 == 0) leadingNum = 1;
    //}else if (reqBytesLength == 3)
    //{
    //    if (pre_ivalue >> 8 == 0) leadingNum = 3;
    //    else if (pre_ivalue >> 16 == 0) leadingNum = 2;
    //    else if (pre_ivalue >> 24 == 0) leadingNum = 1;
    //}else if (reqBytesLength == 1)
    //{
    //    if (pre_ivalue >> 24 == 0) leadingNum = 1;

    //}else if (reqBytesLength == 4)
    //{
    //    if (pre_ivalue == 0) leadingNum = 4;
    //    else if (pre_ivalue >> 8 == 0) leadingNum = 3;
    //    else if (pre_ivalue >> 16 == 0) leadingNum = 2;
    //    else if (pre_ivalue >> 24 == 0) leadingNum = 1;
    //}
    ////if (bi==true) printf("%i:%i:%i:%u\n", blockIdx.x, threadIdx.x, threadIdx.y, leadingNum);
    ////midBytes[bbase+threadIdx.y*blockDim.x+threadIdx.x] = leadingNum; 

    //int midByte_size = reqBytesLength - leadingNum;
    //int midByte_sum = _shfl_scan(midByte_size, sums);
    //uchar4 cur_cvalue = cvalue[threadIdx.y*blockDim.x+threadIdx.x];
    //if (reqBytesLength == 2)
    //{
    //    if (midByte_size == 1){
    //        midBytes[bbase+midByte_sum-1] = cur_cvalue.z; 
    //        //if (bi==true) printf("%i:%i:%i:%u\n", blockIdx.x, threadIdx.x, threadIdx.y, cur_cvalue.z);
    //    }else if (midByte_size == 2){
    //        midBytes[bbase+midByte_sum-1] = cur_cvalue.w; 
    //        //if (bi==true) printf("%i:%i:%i:%u\n", blockIdx.x, threadIdx.x, threadIdx.y, cur_cvalue.z);
    //        midBytes[bbase+midByte_sum-2] = cur_cvalue.z;
    //        //if (bi==true) printf("%i:%i:%i:%u\n", blockIdx.x, threadIdx.x, threadIdx.y, cur_cvalue.w);
    //    }
    //}else if (reqBytesLength == 3)
    //{
    //    if (midByte_size == 1){
    //        midBytes[bbase+midByte_sum-1] = cur_cvalue.y; 
    //    }else if (midByte_size == 2){
    //        midBytes[bbase+midByte_sum-1] = cur_cvalue.z; 
    //        midBytes[bbase+midByte_sum-2] = cur_cvalue.y; 
    //    }else if (midByte_size == 3){
    //        midBytes[bbase+midByte_sum-1] = cur_cvalue.w; 
    //        midBytes[bbase+midByte_sum-2] = cur_cvalue.z; 
    //        midBytes[bbase+midByte_sum-3] = cur_cvalue.y; 
    //    }
    //}else if (reqBytesLength == 1)
    //{
    //    if (midByte_size == 1)
    //        midBytes[bbase+midByte_sum-1] = cur_cvalue.w; 
    //}else if (reqBytesLength == 4)
    //{
    //    if (midByte_size == 1){
    //        midBytes[bbase+midByte_sum-1] = cur_cvalue.x; 
    //    }else if (midByte_size == 2){
    //        midBytes[bbase+midByte_sum-1] = cur_cvalue.y; 
    //        midBytes[bbase+midByte_sum-2] = cur_cvalue.x; 
    //    }else if (midByte_size == 3){
    //        midBytes[bbase+midByte_sum-1] = cur_cvalue.z; 
    //        midBytes[bbase+midByte_sum-2] = cur_cvalue.y; 
    //        midBytes[bbase+midByte_sum-3] = cur_cvalue.x; 
    //    }else if (midByte_size == 4){
    //        midBytes[bbase+midByte_sum-1] = cur_cvalue.w; 
    //        midBytes[bbase+midByte_sum-2] = cur_cvalue.z; 
    //        midBytes[bbase+midByte_sum-3] = cur_cvalue.y; 
    //        midBytes[bbase+midByte_sum-4] = cur_cvalue.x; 
    //    }
    //}

    //if (threadIdx.x==0 && threadIdx.y==0) meta[mbase] = (unsigned char)reqLength;
    //if (threadIdx.x==blockDim.x-1 && threadIdx.y==blockDim.y-1) offsets[obase] = (short)midByte_sum;
    //_IntArray2ByteArray(leadingNum, mbase+1, meta, bi);

}

__global__ void decompress_float(int *offsets, unsigned char *ncBytes, float *data, int bs, size_t nc, size_t mSize, int *test) 
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
    for (int b=bid; b<nc; b+=gridDim.x){
        
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
        if (tidx==0){
            meta[b] = state;
            meta[nb+b*mSize] = cvalue[1].x;
            meta[nb+b*mSize+1] = cvalue[1].y;
            meta[nb+b*mSize+2] = cvalue[1].z;
            meta[nb+b*mSize+3] = cvalue[1].w;
        } 
        //if (tidx==0) test[b] = ivalue[0];
        __syncthreads();                  

        if (state==1){
            //int reqLength = _compute_reqLength(ivalue[0], ivalue[2]);
            //if (tidx==0) test[b] = reqLength;
            //__syncthreads();                  
            //value[tidy*blockDim.x+tidx] = data - medianValue;
            //__syncthreads();                  
            //_decompress_oneBlock(b*bs*sizeof(float), nb+b*mSize+4, b, reqLength, value, ivalue, cvalue, sums, meta, offsets, midBytes, bi);
            bi = false;
        }

    }

}
