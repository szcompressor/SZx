#include <stdio.h>
#include <math.h>
#include "cuSZx_decompress_float.h"

#include <cooperative_groups.h>

namespace cg = cooperative_groups;

__device__ int _deshfl_scan(int lznum, int *sums)
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

__global__ void decompress_float(unsigned char *data, int bs, size_t nc, size_t mSize, int *test) 
{
    int tidx = threadIdx.x;
    int tidy = threadIdx.y;
    int tid = tidy*warpSize+tidx;
    int bid = blockIdx.x;

    float radius, medianValue;
    unsigned mask;
    unsigned char leadingNum;
    extern __shared__ float shared[];
    float* value = shared;
    int* ivalue = (int*)shared;
    uchar4* c4value = (uchar4*)shared;
    unsigned char* cvalue = (unsigned char*)shared;
    int* sums = &ivalue[bs];
    int reqLength;
    uchar4* uc4bytes = (uchar4*)data;
    float* fbytes = (float*)data;
	int reqBytesLength;
	int rightShiftBits;


    bool bi = false;
    if (bid==73) bi=true;
    for (int b=bid; b<nc; b+=gridDim.x){
        c4value[tid] = uc4bytes[b*bs+tid];
        __syncthreads();                  
        medianValue = value[0];
        reqLength = (int)cvalue[4];
        //if (b<2&&tidx==0&&tidy==0) printf("sss%d:%d\n",b, reqLength);
        if (reqLength%8 != 0)
        {
            reqBytesLength = reqLength/8+1;		
            rightShiftBits = 8 - reqLength%8;
        }else{
            reqBytesLength = reqLength/8;		
            rightShiftBits = 0;
        }
        leadingNum = cvalue[5+(tid>>2)];
        leadingNum = (leadingNum >> (6-((tid&0x03)<<1))) & 0x03;
        if (b<1) printf("sss%d:%d,%d,%u\n",reqBytesLength,tidx,tidy,leadingNum);
        int midByte_size = reqBytesLength - leadingNum;
        int midByte_sum = _deshfl_scan(midByte_size, sums);

        uchar4 tmp;
        tmp.x = 0;
        tmp.y = 0;
        tmp.z = 0;
        tmp.w = 0;
        if (reqBytesLength == 2)
        {
            if (midByte_size == 1){
                tmp.z = cvalue[mSize+midByte_sum-1]; 
                //if (bi==true) printf("%i:%i:%i:%u\n", blockIdx.x, threadIdx.x, threadIdx.y, cur_cvalue.z);
            }else if (midByte_size == 2){
                tmp.w = cvalue[mSize+midByte_sum-1]; 
                //if (bi==true) printf("%i:%i:%i:%u\n", blockIdx.x, threadIdx.x, threadIdx.y, cur_cvalue.z);
                tmp.z = cvalue[mSize+midByte_sum-2];
                //if (bi==true) printf("%i:%i:%i:%u\n", blockIdx.x, threadIdx.x, threadIdx.y, cur_cvalue.w);
            }
        }else if (reqBytesLength == 3)
        {
            if (midByte_size == 1){
                tmp.y = cvalue[mSize+midByte_sum-1]; 
            }else if (midByte_size == 2){
                tmp.z = cvalue[mSize+midByte_sum-1]; 
                tmp.y = cvalue[mSize+midByte_sum-2]; 
            }else if (midByte_size == 3){
                tmp.w = cvalue[mSize+midByte_sum-1]; 
                tmp.z = cvalue[mSize+midByte_sum-2]; 
                tmp.y = cvalue[mSize+midByte_sum-3]; 
            }
        }else if (reqBytesLength == 1)
        {
            if (midByte_size == 1)
                tmp.w = cvalue[mSize+midByte_sum-1]; 
        }else if (reqBytesLength == 4)
        {
            if (midByte_size == 1){
                tmp.x = cvalue[mSize+midByte_sum-1]; 
            }else if (midByte_size == 2){
                tmp.y = cvalue[mSize+midByte_sum-1]; 
                tmp.x = cvalue[mSize+midByte_sum-2]; 
            }else if (midByte_size == 3){
                tmp.z = cvalue[mSize+midByte_sum-1]; 
                tmp.y = cvalue[mSize+midByte_sum-2]; 
                tmp.x = cvalue[mSize+midByte_sum-3]; 
            }else if (midByte_size == 4){
                tmp.w = cvalue[mSize+midByte_sum-1]; 
                tmp.z = cvalue[mSize+midByte_sum-2]; 
                tmp.y = cvalue[mSize+midByte_sum-3]; 
                tmp.x = cvalue[mSize+midByte_sum-4]; 
            }
        }
        //if (tidx<2 && tidy==0)
        //    ivalue[tidx] = offsets[b+tidx];
        //__syncthreads();                  
        //obase = ivalue[0];
        //osize = (ivalue[1]-ivalue[0])%4==0 ? (ivalue[1]-ivalue[0])/4 : (ivalue[1]-ivalue[0])/4+1;
        //int* uc4bytes = (int*)(ncBytes+obase); 
        //__syncthreads();                  
        //if (b==0&&tidx==0&&tidy==0) printf("test:%d\n", osize);
        //for (int t=tid; t<osize; t+=blockDim.y*blockDim.x){
        //    int tmp = uc4bytes[t];
        //    ivalue[t] = tmp;
        //    if (b==0) printf("sss:%u\n", t);
        //}
        //__syncthreads();                  
        //medianValue = value[0];
        //if (b==0&&tidx==0&&tidy==0)
        //    printf("median:%f\n", medianValue);


        
        //data = oriData[b*bs+tidy*warpSize+tidx];
        //float Min = data;
        //float Max = data;

        //for (int offset = warpSize/2; offset > 0; offset /= 2) 
        //{
        //    Min = min(Min, __shfl_xor_sync(FULL_MASK, Min, offset));
        //    Max = max(Max, __shfl_xor_sync(FULL_MASK, Max, offset));
        //}
        //if (tidx==0){
        //    value[tidy] = Min;
        //    value[blockDim.y+tidy] = Max;
        //}
        //__syncthreads();                  

        //if (tidy==0){
        //    if (tidx < blockDim.y){
        //        Min = value[tidx];
        //        Max = value[blockDim.y+tidx];
        //    }

        //    mask = __ballot_sync(FULL_MASK, tidx < blockDim.y);
        //    for (int offset = blockDim.y/2; offset > 0; offset /= 2) 
        //    {
        //        Min = min(Min, __shfl_xor_sync(mask, Min, offset));
        //        Max = max(Max, __shfl_xor_sync(mask, Max, offset));
        //    }
        //    
        //    if (tidx==0){
        //        radius = (Max - Min)/2;
        //        value[0] = radius;
        //        value[1] = Min + radius;
        //        value[2] = absErrBound;
        //    }
        //}
        //__syncthreads();                  

        //radius = value[0];
        //medianValue = value[1];
        //state = radius <= absErrBound ? 0 : 1;
        //if (tidx==0){
        //    meta[b] = state;
        //    meta[nb+b*mSize] = cvalue[1].x;
        //    meta[nb+b*mSize+1] = cvalue[1].y;
        //    meta[nb+b*mSize+2] = cvalue[1].z;
        //    meta[nb+b*mSize+3] = cvalue[1].w;
        //} 
        ////if (tidx==0) test[b] = ivalue[0];
        //__syncthreads();                  

        //if (state==1){
        //    //int reqLength = _compute_reqLength(ivalue[0], ivalue[2]);
        //    //if (tidx==0) test[b] = reqLength;
        //    //__syncthreads();                  
        //    //value[tidy*blockDim.x+tidx] = data - medianValue;
        //    //__syncthreads();                  
        //    //_decompress_oneBlock(b*bs*sizeof(float), nb+b*mSize+4, b, reqLength, value, ivalue, cvalue, sums, meta, offsets, midBytes, bi);
        //    bi = false;
        //}

    }

}
