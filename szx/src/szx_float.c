/**
 *  @file szx_float.c
 *  @author Sheng Di, Kai Zhao
 *  @date Aug, 2022
 *  @brief SZ_Init, Compression and Decompression functions
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "szx.h"
#include "szx_float.h"
#include "szx_BytesToolkit.h"
#include "szx_TypeManager.h"
#include <assert.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#if defined(__AVX__) || defined(__AVX2__)  || defined(__AVX512F__)
#include <immintrin.h>
#endif

unsigned char *
SZ_fast_compress_args_with_prediction_float(float *pred, float *data, size_t *outSize, float absErrBound, size_t r5,
                                            size_t r4, size_t r3, size_t r2, size_t r1, float medianValue,
                                            float radius) {
    size_t dataLength = computeDataLength(r5, r4, r3, r2, r1);
    float *delta = (float *) malloc(sizeof(float) * dataLength);
    size_t i = 0;
    for (i = 0; i < dataLength; i++)
        delta[i] = data[i] - pred[i];
    unsigned char *output = SZ_fast_compress_args_unpredictable_float(delta, outSize, absErrBound, r5, r4, r3, r2, r1,
                                                                      medianValue, radius);
    return output;
}

inline void SZ_fast_compress_args_unpredictable_one_block_float(float *oriData, size_t nbEle, float absErrBound,
                                                                unsigned char *outputBytes, int *outSize,
                                                                unsigned char *leadNumberArray_int, float medianValue,
                                                                float radius) {
    size_t totalSize = 0, i = 0;

    int reqLength;

    //compute median, value range, and radius

    short radExpo = getExponent_float(radius);
    computeReqLength_float(absErrBound, radExpo, &reqLength, &medianValue);

    int reqBytesLength = reqLength / 8;
    int resiBitsLength = reqLength % 8;
    int rightShiftBits = 0;

    size_t leadNumberArray_size = nbEle % 4 == 0 ? nbEle / 4 : nbEle / 4 + 1;

    register lfloat lfBuf_pre;
    register lfloat lfBuf_cur;
    lfBuf_pre.ivalue = 0;

    unsigned char *leadNumberArray = outputBytes + 1 + sizeof(float);

    unsigned char *exactMidbyteArray = leadNumberArray + leadNumberArray_size;

    if (resiBitsLength != 0) {
        rightShiftBits = 8 - resiBitsLength;
        reqBytesLength++;
    }
    
    sumReqNbBits += reqLength;
    sumReqNbBytes += reqBytesLength;    

    register unsigned char leadingNum = 0;
    size_t residualMidBytes_size = 0;
    if (sysEndianType == LITTLE_ENDIAN_SYSTEM) {
        if (reqBytesLength == 2) {
            for (i = 0; i < nbEle; i++) {
                leadingNum = 0;
                lfBuf_cur.value = oriData[i] - medianValue;

                lfBuf_cur.ivalue = lfBuf_cur.ivalue >> rightShiftBits;

                lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;

                if (lfBuf_pre.ivalue >> 8 == 0)
                    leadingNum = 3;
                else if (lfBuf_pre.ivalue >> 16 == 0)
                    leadingNum = 2;
                else if (lfBuf_pre.ivalue >> 24 == 0)
                    leadingNum = 1;

                leadNumberArray_int[i] = leadingNum;

                if (leadingNum == 0) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[2];
                    exactMidbyteArray[residualMidBytes_size + 1] = lfBuf_cur.byte[3];
                    residualMidBytes_size += 2;
                } else if (leadingNum == 1) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[2];
                    residualMidBytes_size++;
                }

                lfBuf_pre = lfBuf_cur;
            }
        } else if (reqBytesLength == 3) {
            for (i = 0; i < nbEle; i++) {
                leadingNum = 0;
                lfBuf_cur.value = oriData[i] - medianValue;

                lfBuf_cur.ivalue = lfBuf_cur.ivalue >> rightShiftBits;

                lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;

                if (lfBuf_pre.ivalue >> 8 == 0)
                    leadingNum = 3;
                else if (lfBuf_pre.ivalue >> 16 == 0)
                    leadingNum = 2;
                else if (lfBuf_pre.ivalue >> 24 == 0)
                    leadingNum = 1;

                leadNumberArray_int[i] = leadingNum;

                if (leadingNum == 0) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[1];
                    exactMidbyteArray[residualMidBytes_size + 1] = lfBuf_cur.byte[2];
                    exactMidbyteArray[residualMidBytes_size + 2] = lfBuf_cur.byte[3];
                    residualMidBytes_size += 3;
                } else if (leadingNum == 1) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[1];
                    exactMidbyteArray[residualMidBytes_size + 1] = lfBuf_cur.byte[2];
                    residualMidBytes_size += 2;
                } else if (leadingNum == 2) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[1];
                    residualMidBytes_size++;
                }

                lfBuf_pre = lfBuf_cur;
            }
        } else if (reqBytesLength == 1) {
            for (i = 0; i < nbEle; i++) {
                leadingNum = 0;
                lfBuf_cur.value = oriData[i] - medianValue;

                lfBuf_cur.ivalue = lfBuf_cur.ivalue >> rightShiftBits;

                lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;

                if (lfBuf_pre.ivalue >> 8 == 0)
                    leadingNum = 3;
                else if (lfBuf_pre.ivalue >> 16 == 0)
                    leadingNum = 2;
                else if (lfBuf_pre.ivalue >> 24 == 0)
                    leadingNum = 1;

                leadNumberArray_int[i] = leadingNum;

                if (leadingNum == 0) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[3];
                    residualMidBytes_size++;
                }

                lfBuf_pre = lfBuf_cur;
            }
        } else //reqBytesLength == 4
        {
            for (i = 0; i < nbEle; i++) {
                leadingNum = 0;
                lfBuf_cur.value = oriData[i] - medianValue;

                lfBuf_cur.ivalue = lfBuf_cur.ivalue >> rightShiftBits;

                lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;

                if (lfBuf_pre.ivalue >> 8 == 0)
                    leadingNum = 3;
                else if (lfBuf_pre.ivalue >> 16 == 0)
                    leadingNum = 2;
                else if (lfBuf_pre.ivalue >> 24 == 0)
                    leadingNum = 1;

                leadNumberArray_int[i] = leadingNum;

                if (leadingNum == 0) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[0];
                    exactMidbyteArray[residualMidBytes_size + 1] = lfBuf_cur.byte[1];
                    exactMidbyteArray[residualMidBytes_size + 2] = lfBuf_cur.byte[2];
                    exactMidbyteArray[residualMidBytes_size + 3] = lfBuf_cur.byte[3];
                    residualMidBytes_size += 4;
                } else if (leadingNum == 1) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[0];
                    exactMidbyteArray[residualMidBytes_size + 1] = lfBuf_cur.byte[1];
                    exactMidbyteArray[residualMidBytes_size + 2] = lfBuf_cur.byte[2];
                    residualMidBytes_size += 3;
                } else if (leadingNum == 2) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[0];
                    exactMidbyteArray[residualMidBytes_size + 1] = lfBuf_cur.byte[1];
                    residualMidBytes_size += 2;
                } else //leadingNum == 3
                {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[0];
                    residualMidBytes_size++;
                }

                lfBuf_pre = lfBuf_cur;
            }
        }

        convertIntArray2ByteArray_fast_2b_args(leadNumberArray_int, nbEle, leadNumberArray);
        int k = 0;

        unsigned char reqLengthB = (unsigned char) reqLength;
        outputBytes[k] = reqLengthB;
        k++;
        floatToBytes(&(outputBytes[k]), medianValue);
        k += sizeof(float);
        //sizeToBytes(&(outputBytes[k]), leadNumberArray_size);
        //outputBytes[k] = leadNumberArray_size;  //leadNumberArray_size can be calculated based on block size (=blockSize/4)

        totalSize = 1 + sizeof(float) + leadNumberArray_size + residualMidBytes_size;
        sum_leadNumberArray_size += leadNumberArray_size;
        sum_residualMidBytes_size += residualMidBytes_size;
        
        int s_actual_leadNumbers = 0;
        for(int j = 0;j<nbEle;j++)
		{
			if(leadNumberArray_int[j] >= reqBytesLength)
				s_actual_leadNumbers += reqBytesLength;
			else
				s_actual_leadNumbers += leadNumberArray_int[j];
		}	
        
        sum_actual_leadNumbers += s_actual_leadNumbers;
        //printf("sum_resiMidBytes=%f\n", (reqBytesLength*128-sum_actual_leadNumbers));
        //printf("test");
    } else {

    }

    *outSize = totalSize;

}

float computeRadiusBuffer_float(float *oriData, size_t nbEle, int samplingRate, int blockSize, float** radiusArray, float** mediusArray, float** buffer)
{
	size_t i = 0, j = 0;
	size_t offset = 0;
	size_t stride = samplingRate*blockSize;
	size_t nbBlocks = nbEle/stride;
	*radiusArray = (float*)malloc(sizeof(float)*(nbBlocks+1));
	*mediusArray = (float*)malloc(sizeof(float)*(nbBlocks+1));
	*buffer = (float*)malloc(sizeof(float)*(nbBlocks+1)*blockSize);
	float* p = *buffer;
	float g_min = oriData[offset], g_max = oriData[offset];
	for(i=0;i<nbBlocks;i++)
	{
		memcpy(p, &oriData[offset], blockSize*sizeof(float));
		float min = oriData[offset];
		float max = oriData[offset];
        for (j = 1; j < blockSize; j++) {
		float v = oriData[offset + j];
		if (min > v)
			min = v;
		else if (max < v)
			max = v;
		}
		if(g_max < max)
			g_max = max;
		if(g_min > min)
			g_min = min;
        float valueRange = max - min;
        float radius = valueRange / 2;		
		(*radiusArray)[i] = radius;
		(*mediusArray)[i] = min + radius;
		offset += stride;
		p += blockSize;
	}
/*	size_t remainder = nbEle%stride;
	if(remainder)
	{
		memcpy(*buffer, &oriData[offset], remainder);
		float min = oriData[offset];
		float max = oriData[offset];
        for (j = 1; j < remainder; j++) {
		float v = oriData[offset + j];
		if (min > v)
			min = v;
		else if (max < v)
			max = v;
		}
        float valueRange = max - min;
        float radius = valueRange / 2;		
		(*radiusArray)[i] = radius;     
		(*mediusArray)[i] = min + radius;
		return nbBlocks+1;		
	}
	else
		return nbBlocks;*/
		
	return g_max - g_min;
}

size_t computeStateMedianRadius_float(float *oriData, size_t nbEle, float absErrBound, int blockSize,
                                      unsigned char *stateArray, float *medianArray, float *radiusArray) {
    size_t nbConstantBlocks = 0;
    size_t i = 0, j = 0;
    size_t nbBlocks = nbEle / blockSize;
    size_t offset = 0;

    for (i = 0; i < nbBlocks; i++) {
        float min = oriData[offset];
        float max = oriData[offset];
        for (j = 1; j < blockSize; j++) {
            float v = oriData[offset + j];
            if (min > v)
                min = v;
            else if (max < v)
                max = v;
        }
        float valueRange = max - min;
        float radius = valueRange / 2;
        float medianValue = min + radius;

        if (radius <= absErrBound) {
            stateArray[i] = 0;
            nbConstantBlocks++;
        } else
            stateArray[i] = 1;

        //stateArray[i] = radius <= absErrBound ? 0 : 1;
        medianArray[i] = medianValue;
        radiusArray[i] = radius;
        offset += blockSize;
    }

    int remainCount = nbEle % blockSize;
    if (remainCount != 0) {
        float min = oriData[offset];
        float max = oriData[offset];
        for (j = 1; j < remainCount; j++) {
            float v = oriData[offset + j];
            if (min > v)
                min = v;
            else if (max < v)
                max = v;
        }
        float valueRange = max - min;
        float radius = valueRange / 2;
        float medianValue = min + radius;
        if (radius <= absErrBound) {
            stateArray[i] = 0;
            nbConstantBlocks++;
        } else
            stateArray[i] = 1;
        medianArray[i] = medianValue;
        radiusArray[i] = radius;
    }
    return nbConstantBlocks;
}


void max_min_float(float *x, int n, float *tmp_max, float *tmp_min) {
    for (size_t i = 0; i < n; i++) {
        if (x[i] > *tmp_max) {
            *tmp_max = x[i];
        }
        if (x[i] < *tmp_min) {
            *tmp_min = x[i];
        }
    }
}

void simd_max_min_float(float *x, int n, float *tmp_max, float *tmp_min) {
    *tmp_max = x[0];
    *tmp_min = x[0];
#ifdef  __AVX512F__
    //    printf("use avx512, n=%d \n", n);
    int n16 = n & -16, i = 0, j=0;
    if (n > 16) {
        float *ptr_x = x;
        __m512 max1 = _mm512_loadu_ps(ptr_x);
//        __m512 max2 = _mm512_loadu_ps(ptr_x + 16);
        __m512 min1 = max1;
//        __m512 min2 = max2;
        __m512 tmp1;
//        __m512 tmp2;
        for (; i < n16; i += 16) {
            tmp1 = _mm512_loadu_ps(ptr_x);
            max1 = _mm512_max_ps(tmp1, max1);
            min1 = _mm512_min_ps(tmp1, min1);
//            tmp2 = _mm512_loadu_ps(ptr_x+16);
//            max2 = _mm512_max_ps(tmp2, max2);
//            min2 = _mm512_min_ps(tmp2, min2);
            ptr_x += 16;
        }
//        max1 = _mm512_max_ps(max1, max2);
//        min1 = _mm512_min_ps(min1, min2);
          __m256 max256 = _mm256_max_ps(_mm512_extractf32x8_ps(max1,0), _mm512_extractf32x8_ps(max1,1));
          __m128 max128 = _mm_max_ps(_mm256_extractf128_ps(max256,0), _mm256_extractf128_ps(max256,1));
          __m256 min256 = _mm256_min_ps(_mm512_extractf32x8_ps(min1,0), _mm512_extractf32x8_ps(min1,1));
          __m128 min128 = _mm_min_ps(_mm256_extractf128_ps(min256,0), _mm256_extractf128_ps(min256,1));
          for (j=0;j<4;j++){
            *tmp_max = *tmp_max < max128[j] ? max128[j] : *tmp_max;
            *tmp_min = *tmp_min > min128[j] ? min128[j] : *tmp_min;
          }

        if ( i < n ) {
            max_min_float(ptr_x, n - i, tmp_max, tmp_min);
        }
    } else {
        max_min_float(x, n, tmp_max, tmp_min);
    }
#elif __AVX2__
//        printf("use avx2, n=%d \n", n);
    //    fflush(stdout);
    int n16 = n & -16, i = 0;
    if (n > 16) {
        float *ptr_x = x;
        __m256 max1 = _mm256_loadu_ps(ptr_x);
        __m256 max2 = _mm256_loadu_ps(ptr_x + 8);
        __m256 min1 = max1;
        __m256 min2 = max2;
        for (; i < n16; i += 16) {
            max1 = _mm256_max_ps(_mm256_loadu_ps(ptr_x), max1);
            min1 = _mm256_min_ps(_mm256_loadu_ps(ptr_x), min1);
            max2 = _mm256_max_ps(_mm256_loadu_ps(ptr_x + 8), max2);
            min2 = _mm256_min_ps(_mm256_loadu_ps(ptr_x + 8), min2);
            ptr_x += 16;
        }
//        printf("%d %d %d\n", n, n16, i);
//        exit(0);
        max1 = _mm256_max_ps(max1, max2);
        min1 = _mm256_min_ps(min1, min2);
        for (int j = 0; j < 8; j++) {
            *tmp_max = *tmp_max < max1[j] ? max1[j] : *tmp_max;
            *tmp_min = *tmp_min > min1[j] ? min1[j] : *tmp_min;
        }
        if ( i < n ) {
            max_min_float(ptr_x, n - i, tmp_max, tmp_min);
        }
    } else {
        max_min_float(x, n, tmp_max, tmp_min);
    }
#else
    max_min_float(x, n, tmp_max, tmp_min);
#endif
}

void computeStateMedianRadius_float2(float *oriData, size_t nbEle, float absErrBound,
                                     unsigned char *state, float *median, float *radius) {
     float min = oriData[0];
     float max = oriData[0];
     simd_max_min_float(oriData, nbEle, &max, &min);

    float valueRange = max - min;
    *radius = valueRange / 2;
    *median = min + *radius;

    if (*radius <= absErrBound) {
        *state = 0;
    } else {
        *state = 1;
    }
}

float computeValueRange_float(float* oriData, size_t length, float* radius, float* medianValue)
{
	//compute value range
	float min = oriData[0];
	float max = oriData[0];
	for(size_t i=0;i<length;i++)
	{
		float v = oriData[i];
		if(min>v)
			min = v;
		else if(max<v)
			max = v;
	}
	float valueRange = max - min;
	if(radius!=NULL)
	{
		*radius = valueRange/2;
		*medianValue = min + *radius;			
	}

	return valueRange;
}

float estimateCRbasedonErrorBound_buffered_float(float errorBound, float* buffer, float* medianArray, float* radiusArray, int samplingRate, int blockSize, size_t nbEle, 
size_t *sumReqNbBytes, size_t *sum_actual_leadNumbers)
{
	size_t nbSampledBlocks = nbEle/(blockSize*samplingRate); //ignored the last remainder block
	size_t nbConstantBlocks = 0;
	*sum_actual_leadNumbers = 0;
	float metadata = 9.0*blockSize/nbEle;
	float block_cost = 33.0/8;	
	size_t i = 0, j = 0;
	*sumReqNbBytes = 0;
	for(i=0;i<nbSampledBlocks;i++) //ignored the last remainder block
	{
		size_t offset = i*blockSize;
		float medianValue = medianArray[i];
		float radius = radiusArray[i];
		float *data = &buffer[offset];
//		if(i==463)
//			printf("i=%zu\n", i);
		
		if (radius <= errorBound) {
			nbConstantBlocks++;
		}
		else //non-constant
		{
			int reqLength;
			short radExpo = getExponent_float(radius);
			computeReqLength_float(errorBound, radExpo, &reqLength, &medianValue);
			int reqBytesLength = reqLength / 8;
			int resiBitsLength = reqLength % 8;
			int rightShiftBits = 0;	

			if (resiBitsLength != 0) {
				rightShiftBits = 8 - resiBitsLength;				
				reqBytesLength++;
			}
			
			//printf("%d\n",reqBytesLength);
			*sumReqNbBytes+=	reqBytesLength;
			
			register lfloat lfBuf_pre;
			register lfloat lfBuf_cur;
			lfBuf_pre.ivalue = 0;
			register unsigned char leadingNum = 0;			
			//int leadingNum_Array[128];
			int s_actual_leadNumbers = 0;
			for(j=0;j<blockSize;j++)
			{
                lfBuf_cur.value = data[j] - medianValue;

                lfBuf_cur.ivalue = lfBuf_cur.ivalue >> rightShiftBits;

                lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;
				
				leadingNum = 0;
                if (lfBuf_pre.ivalue >> 8 == 0)
                    leadingNum = 3;
                else if (lfBuf_pre.ivalue >> 16 == 0)
                    leadingNum = 2;
                else if (lfBuf_pre.ivalue >> 24 == 0)
                    leadingNum = 1;
                   
               // leadingNum_Array[j] = leadingNum;    
				if(leadingNum >= reqBytesLength)
					s_actual_leadNumbers += reqBytesLength;
				else
					s_actual_leadNumbers += leadingNum;     
					
				lfBuf_pre = lfBuf_cur;	               
			}
			
			*sum_actual_leadNumbers += s_actual_leadNumbers;			
		}		
	}
	
	float avgReqNbBytes = 1.0*(*sumReqNbBytes)/(nbSampledBlocks - nbConstantBlocks);
	float avg_actual_lead = 1.0*(*sum_actual_leadNumbers)/(nbSampledBlocks - nbConstantBlocks);
	
	float p_lambda = 1.0*nbConstantBlocks/nbSampledBlocks;
	
	float estimatedCR = 4*blockSize/(metadata + block_cost+(1 + (0.25+avgReqNbBytes)*blockSize - avg_actual_lead)*(1 - p_lambda));
	//printf("----> sum_actual_leadNumbers=%zu, nbSampledBlocks = %zu, nbConstantBlocks = %zu, sumReqNbBytes = %zu, avgReqNbBytes=%f, avg_actual_lead=%f, p_lambda=%f\n", *sum_actual_leadNumbers, nbSampledBlocks, nbConstantBlocks, *sumReqNbBytes, avgReqNbBytes, avg_actual_lead, p_lambda);
	return estimatedCR;	
}

float estimateCRbasedonErrorBound_float(float errorBound, float* data, int blockSize, size_t nbEle)
{
	float metadata = 9.0*blockSize/nbEle;
	float block_cost = 33.0/8;
	size_t nbBlocks = nbEle/blockSize;
	size_t nbConstantBlocks = 0;
	size_t sum_actual_leadNumbers = 0;
	int samplingRate = 10; //10 means 1/10 data
	size_t nbBlocks_sampled = nbBlocks/samplingRate;
	//estimate p_lambda, avgReqNbBytes, and avg_actual_lead using sampling
	size_t i = 0, j = 0;
	size_t sumReqNbBytes = 0;
	for(i=0;i<nbBlocks;i+=samplingRate) //10 means 1/10 sampling rate
	{
		size_t offset = i*blockSize;
//		if(data[offset]!=0)
//			printf("i=%zu\n", i);
		float min = data[offset];
		float max = data[offset];
		for (j = 1; j < blockSize; j++) {
			float v = data[offset + j];	
			if (min > v)
				min = v;
			if (max < v)
				max = v;
		}
		float valueRange = max - min;
		float radius = valueRange / 2;
		float medianValue = min + radius;

		if (radius <= errorBound) {
			nbConstantBlocks++;
		}
		else //non-constant
		{
			int reqLength;
			short radExpo = getExponent_float(radius);
			computeReqLength_float(errorBound, radExpo, &reqLength, &medianValue);

			int reqBytesLength = reqLength / 8;
			int resiBitsLength = reqLength % 8;
			int rightShiftBits = 0;	

			if (resiBitsLength != 0) {
				rightShiftBits = 8 - resiBitsLength;				
				reqBytesLength++;
			}
			
			//printf("%d\n",reqBytesLength);
			sumReqNbBytes+=	reqBytesLength;
			
			register lfloat lfBuf_pre;
			register lfloat lfBuf_cur;
			lfBuf_pre.ivalue = 0;
			register unsigned char leadingNum = 0;			
			//int leadingNum_Array[128];
			int s_actual_leadNumbers = 0;
			for(j=0;j<blockSize;j++)
			{
                lfBuf_cur.value = data[offset+j] - medianValue;

                lfBuf_cur.ivalue = lfBuf_cur.ivalue >> rightShiftBits;

                lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;
				
				leadingNum = 0;
                if (lfBuf_pre.ivalue >> 8 == 0)
                    leadingNum = 3;
                else if (lfBuf_pre.ivalue >> 16 == 0)
                    leadingNum = 2;
                else if (lfBuf_pre.ivalue >> 24 == 0)
                    leadingNum = 1;
                   
               // leadingNum_Array[j] = leadingNum;    
				if(leadingNum >= reqBytesLength)
					s_actual_leadNumbers += reqBytesLength;
				else
					s_actual_leadNumbers += leadingNum;     
					
				lfBuf_pre = lfBuf_cur;	               
			}
			
			sum_actual_leadNumbers += s_actual_leadNumbers;
		}
	}
	
	float avgReqNbBytes = 1.0*sumReqNbBytes/(nbBlocks_sampled - nbConstantBlocks);
	float avg_actual_lead = 1.0*sum_actual_leadNumbers/(nbBlocks_sampled - nbConstantBlocks);
	
	float p_lambda = 1.0*nbConstantBlocks/nbBlocks_sampled;
	
	float estimatedCR = 4*blockSize/(metadata + block_cost+(1 + (0.25+avgReqNbBytes)*blockSize - avg_actual_lead)*(1 - p_lambda));
	//printf("----> sum_actual_leadNumbers=%zu, nbBlocks_sampled=%zu, nbConstantBlocks = %zu, sumReqNbBytes = %zu, avgReqNbBytes=%f, avg_actual_lead=%f, p_lambda=%f\n", sum_actual_leadNumbers, nbBlocks_sampled, nbConstantBlocks, sumReqNbBytes, avgReqNbBytes, avg_actual_lead, p_lambda);
	exit(0);
	return estimatedCR;
}

float estimateErrorBoundbasedonCR_buffered_float(float targetCompressionRatio, float tolerance, int samplingRate, float initErrorBound, int blockSize, size_t nbEle, 
float* buffer, float* medianArray, float* radiusArray)
{
	int stride = 8;
	//printf("initErrBound = %.10G\n", initErrorBound);
	int nbIterations = 10;
	float errorBound = initErrorBound;
	size_t sumReqNbBytes[20] = {0};
	size_t sum_actual_leadNumbers[20] = {0};
	float errorBounds[20];
	float ratios[20];

	int i = 0, j = 0, k = 0;	
	//search for a valid range for error bounds

	size_t sumReqNbB = 0, sum_actual_leadNum = 0;
	float CR = estimateCRbasedonErrorBound_buffered_float(errorBound, buffer, medianArray, radiusArray, samplingRate, blockSize, nbEle, &sumReqNbB, &sum_actual_leadNum);	
	sumReqNbBytes[k] = sumReqNbB;
	sum_actual_leadNumbers[k] = sum_actual_leadNum;
	errorBounds[k] = errorBound;
	ratios[k] = CR;
	k++; //k indicates the actual number of iterations
	
	if(fabs(targetCompressionRatio-CR) < tolerance)
		return errorBound;
	float left_error = 0, right_error = 0;	
	int foundFlag = 0;
	float result_error = 0;
	float targetTolerance = targetCompressionRatio*tolerance;
//	float preCR = CR;

	//leftward search
	if(targetCompressionRatio<CR)
	{
		//printf("targetCompressionRatio<CR\n");
		left_error = errorBound/stride;
		right_error = errorBound;				
		for(i=0;i<nbIterations;i++)
		{
			CR = estimateCRbasedonErrorBound_buffered_float(left_error, buffer, medianArray, radiusArray, samplingRate, blockSize, nbEle, &sumReqNbB, &sum_actual_leadNum);	
			sumReqNbBytes[k] = sumReqNbB;
			sum_actual_leadNumbers[k] = sum_actual_leadNum;
			errorBounds[k] = left_error;
			ratios[k] = CR;
			k++;
			//printf("errorBound=%.10G, CR=%f\n", left_error, CR);
			if(fabs(CR-targetCompressionRatio) < targetTolerance)// || fabs(CR-preCR)<targetTolerance/2) //stop and convergence condition
			{
				foundFlag = 1;
				result_error = left_error;
				break;
			}	
			if(targetCompressionRatio > CR)	
				break;
			right_error = left_error;
			left_error = left_error/stride;
			//preCR = CR;
		}
	}
	else if(targetCompressionRatio>CR) //rightward search
	{
		//printf("targetCompressionRatio>CR\n");
		left_error = errorBound;
		right_error = errorBound*stride;				
		for(i=0;i<nbIterations;i++)
		{	
			CR = estimateCRbasedonErrorBound_buffered_float(right_error, buffer, medianArray, radiusArray, samplingRate, blockSize, nbEle, &sumReqNbB, &sum_actual_leadNum);	
			sumReqNbBytes[k] = sumReqNbB;
			sum_actual_leadNumbers[k] = sum_actual_leadNum;
			errorBounds[k] = right_error;
			ratios[k] = CR;
			k++;			
			//printf("errorBound=%.10G, CR=%f\n", right_error, CR);			
			if(fabs(targetCompressionRatio-CR) < targetTolerance)// || fabs(CR-preCR)<targetTolerance/2) //convergence condition
			{
				foundFlag = 1;
				result_error = right_error;
				break;
			}
			if(targetCompressionRatio < CR)	
				break;
			left_error = right_error;
			right_error = right_error*stride;
		}
		//preCR = CR;
	}
	
	//printf("left_error = %.10G, right_error = %.10G\n", left_error, right_error);
	//printf("foundFlag=%d, stage 1 iterations = %d, result_error = %.10G, k = %d\n", foundFlag, i, result_error, k);	
	//binary search
	if(foundFlag==0)
	{
		for(i=0;i<nbIterations;i++)
		{
			errorBound = (left_error+right_error)/2;
			CR = estimateCRbasedonErrorBound_buffered_float(errorBound, buffer, medianArray, radiusArray, samplingRate, blockSize, nbEle, &sumReqNbB, &sum_actual_leadNum);
			
			//check if the CR appeared before to avoid pingpang effect
			int appeared = 0;
			for(j=0;j<k;j++)
			{
				if(sumReqNbBytes[j]==sumReqNbB && sum_actual_leadNumbers[j] == sum_actual_leadNum) //this compression ratio appeared before
				{
					appeared = 1;
					break;
				}	
			}
			float min_ratio_error = 10000000;
			int index = 0;
			for(j=0;j<k;j++) //use the best error bound
			{
				float r = ratios[j];
				float e = fabs(r-targetCompressionRatio);
				if(e < min_ratio_error) //this compression ratio appeared before
				{
					min_ratio_error = e;
					index = j;
				}	
			}		
			if(appeared)
			{
				result_error = errorBounds[index];
				break;
			}	
			
			//record the sumReqNBB and sum_actual_leadNum
			sumReqNbBytes[k] = sumReqNbB;
			sum_actual_leadNumbers[k] = sum_actual_leadNum;
			errorBounds[k] = errorBound;
			ratios[k] = CR;
			k++;			
			//printf("errorBound=%.10G, CR=%f\n", errorBound, CR);			
			//check if the CR meets the targetTolerance
			if(fabs(CR-targetCompressionRatio) < targetTolerance)// || fabs(CR-preCR)<targetTolerance/2) //convergence condition
			{
				foundFlag = 1;
				result_error = errorBound;
				break;
			} 
			if(targetCompressionRatio > CR)
				left_error = errorBound;
			else
				right_error = errorBound;
			//preCR = CR;	
		}
	//printf("foundFlag=%d, stage 2 iterations = %d, result_error = %.10G, k=%d\n", foundFlag, i, result_error,k);						
	}
	return result_error;
}

float estimateErrorBoundbasedonCR_float(float targetCompressionRatio, float tolerance, float* data, float initErrorBound, int blockSize, size_t nbEle)
{
	float errorBound = initErrorBound;
	//search for a valid range for error bounds
	int i = 0;
	float CR = estimateCRbasedonErrorBound_float(errorBound, data, blockSize, nbEle);	
	int nbIterations = 15;
	if(fabs(targetCompressionRatio-CR) < tolerance)
		return errorBound;
	float left_error = 0, right_error = 0;	
	int foundFlag = 0;
	float result_error = 0;
	float targetTolerance = targetCompressionRatio*tolerance;
	float preCR = CR;

	//leftward search
	if(targetCompressionRatio<CR)
	{
		//printf("targetCompressionRatio<CR\n");
		left_error = errorBound/8;
		right_error = errorBound;				
		for(i=0;i<nbIterations;i++)
		{
			CR = estimateCRbasedonErrorBound_float(left_error, data, blockSize, nbEle);
			//printf("errorBound=%.30G, CR=%f\n", left_error, CR);
			if(fabs(CR-targetCompressionRatio) < targetTolerance || fabs(CR-preCR)<targetTolerance/2)
			{
				foundFlag = 1;
				result_error = left_error;
				break;
			}	
			if(targetCompressionRatio > CR)	
				break;
			right_error = left_error;
			left_error = left_error/8;
			preCR = CR;
		}
	}
	else if(targetCompressionRatio>CR) //rightward search
	{
		//printf("targetCompressionRatio>CR\n");
		left_error = errorBound;
		right_error = errorBound*8;				
		for(i=0;i<nbIterations;i++)
		{
			CR = estimateCRbasedonErrorBound_float(right_error, data, blockSize, nbEle);
			//printf("errorBound=%.30G, CR=%f\n", right_error, CR);			
			if(fabs(targetCompressionRatio-CR) < targetTolerance || fabs(CR-preCR)<targetTolerance/2)
			{
				foundFlag = 1;
				result_error = right_error;
				break;
			}
			if(targetCompressionRatio < CR)	
				break;
			left_error = right_error;
			right_error = right_error*8;
		}
		preCR = CR;
	}
	
	//printf("foundFlag=%d, stage 1 iterations = %d\n", foundFlag, i);	
	//binary search
	if(foundFlag==0)
	{
		for(i=0;i<nbIterations;i++)
		{
			errorBound = (left_error+right_error)/2;
			CR = estimateCRbasedonErrorBound_float(errorBound, data, blockSize, nbEle);
			//printf("errorBound=%.30G, CR=%f\n", errorBound, CR);			
			if(fabs(CR-targetCompressionRatio) < targetTolerance || fabs(CR-preCR)<targetTolerance/2)
			{
				foundFlag = 1;
				result_error = errorBound;
				break;
			} 
			if(targetCompressionRatio > CR)
				left_error = errorBound;
			else
				right_error = errorBound;
			preCR = CR;	
		}
		//printf("foundFlag=%d, stage 2 iterations = %d\n", foundFlag, i);					
	}
	return result_error;
}

unsigned char *
SZ_fast_compress_args_unpredictable_blocked_fixed_rate_float(float *oriData, size_t *outSize, float compressionRatio, float tolerance, size_t nbEle,
                                                  int blockSize) 
{
    *outSize = 0;
    
    size_t maxPreservedBufferSize = sizeof(float) * nbEle/compressionRatio*(1+tolerance*2); //assume that the compressed data size would not exceed the original size
    unsigned char *outputBytes = (unsigned char *) malloc(maxPreservedBufferSize);
    //memset(outputBytes, 0, maxPreservedBufferSize);

	int samplingStride = 10;
	float* radiusArray = NULL;
	float* mediusArray = NULL;
	float* buffer = NULL;
	float approximateValueRange = computeRadiusBuffer_float(oriData, nbEle, samplingStride, blockSize, &radiusArray, &mediusArray, &buffer);
//	size_t stride = samplingStride*blockSize;
//	size_t nbBlocks = nbEle/stride;
//	int status = 0;
//	writeFloatData(buffer, nbBlocks*blockSize, "1.csv", &status);
	
	float initErrorBound = approximateValueRange*1E-2;

	float errorBound = estimateErrorBoundbasedonCR_buffered_float(compressionRatio, tolerance, samplingStride, initErrorBound, blockSize, nbEle, buffer, mediusArray, radiusArray);
	
	//unsigned char* bytes = SZ_fast_compress_args_unpredictable_blocked_float(oriData, outSize, errorBound, nbEle, blockSize);
	SZ_fast_compress_args_unpredictable_blocked_float2(oriData, outSize, outputBytes, errorBound, nbEle, blockSize);
                    
    return outputBytes;                                                                    
}

void SZ_fast_compress_args_unpredictable_blocked_fixed_rate_float2(float *oriData, size_t *outSize, unsigned char* outputBytes, float compressionRatio, float tolerance, size_t nbEle,
                                                  int blockSize) 
{
    *outSize = 0;
    
    //size_t maxPreservedBufferSize = sizeof(float) * nbEle/compressionRatio*(1+tolerance*2); //assume that the compressed data size would not exceed the original size
    //unsigned char *outputBytes = (unsigned char *) malloc(maxPreservedBufferSize);
    //memset(outputBytes, 0, maxPreservedBufferSize);

	int samplingStride = 10;
	float* radiusArray = NULL;
	float* mediusArray = NULL;
	float* buffer = NULL;
	float approximateValueRange = computeRadiusBuffer_float(oriData, nbEle, samplingStride, blockSize, &radiusArray, &mediusArray, &buffer);
//	size_t stride = samplingStride*blockSize;
//	size_t nbBlocks = nbEle/stride;
//	int status = 0;
//	writeFloatData(buffer, nbBlocks*blockSize, "1.csv", &status);
	
	float initErrorBound = approximateValueRange*1E-2;

	float errorBound = estimateErrorBoundbasedonCR_buffered_float(compressionRatio, tolerance, samplingStride, initErrorBound, blockSize, nbEle, buffer, mediusArray, radiusArray);
	
	//unsigned char* bytes = SZ_fast_compress_args_unpredictable_blocked_float(oriData, outSize, errorBound, nbEle, blockSize);
	SZ_fast_compress_args_unpredictable_blocked_float2(oriData, outSize, outputBytes, errorBound, nbEle, blockSize);
                                                                                         
}

		
size_t SZ_fast_compress_args_unpredictable_blocked_args_float(float *oriData, unsigned char* outputBytes, float absErrBound, size_t nbEle,
                                                  int blockSize)
{
    float *op = oriData;

    size_t outSize = 0;
    unsigned char *leadNumberArray_int = (unsigned char *) malloc(blockSize * sizeof(int));

    size_t i = 0;
    int oSize = 0;

    size_t nbBlocks = nbEle / blockSize;
    size_t remainCount = nbEle % blockSize;
    size_t stateNBBytes =
            remainCount == 0 ? (nbBlocks % 8 == 0 ? nbBlocks / 8 : nbBlocks / 8 + 1) : ((nbBlocks + 1) % 8 == 0 ?
                                                                                        (nbBlocks + 1) / 8 :
                                                                                        (nbBlocks + 1) / 8 + 1);
    size_t actualNBBlocks = remainCount == 0 ? nbBlocks : nbBlocks + 1;

    unsigned char *stateArray = (unsigned char *) malloc(actualNBBlocks);
    float *medianArray = (float *) malloc(actualNBBlocks * sizeof(float));
    float *radiusArray = (float *) malloc(actualNBBlocks * sizeof(float));

    size_t nbConstantBlocks = computeStateMedianRadius_float(oriData, nbEle, absErrBound, blockSize, stateArray,
                                                             medianArray, radiusArray);
                                                             
    //printf("nbConstantBlocks=%zu, %f\n", nbConstantBlocks, 1.0*nbConstantBlocks/actualNBBlocks);

    unsigned char *r = outputBytes; // + sizeof(size_t) + stateNBBytes;
    r[0] = SZx_VER_MAJOR;
    r[1] = SZx_VER_MINOR;
    r[2] = 1;
    r[3] = 0; // indicates this is not a random access version
    r[4] = (unsigned char) blockSize;
    r = r + 5; //1 byte
    sizeToBytes(r, nbConstantBlocks);
    r += sizeof(size_t); //r is the starting address of 'stateNBBytes'

    unsigned char *p = r + stateNBBytes; //p is the starting address of constant median values.
    unsigned char *q =
            p + sizeof(float) * nbConstantBlocks; //q is the starting address of the non-constant data sblocks
    //3: versions, 1: metadata: state, 1: metadata: blockSize, sizeof(size_t): nbConstantBlocks, ....
    outSize += (3 + 1 + 1 + sizeof(size_t) + stateNBBytes + sizeof(float) * nbConstantBlocks);

    //printf("nbConstantBlocks = %zu, percent = %f\n", nbConstantBlocks, 1.0f*(nbConstantBlocks*blockSize)/nbEle);

    for (i = 0; i < nbBlocks; i++, op += blockSize) {
        if (stateArray[i]) {
            SZ_fast_compress_args_unpredictable_one_block_float(op, blockSize, absErrBound, q, &oSize,
                                                                leadNumberArray_int, medianArray[i], radiusArray[i]);
            q += oSize;
            outSize += oSize;
        } else {
            floatToBytes(p, medianArray[i]);
            p += sizeof(float);
        }
    }

    if (remainCount != 0) {
        if (stateArray[i]) {
            SZ_fast_compress_args_unpredictable_one_block_float(op, remainCount, absErrBound, q, &oSize,
                                                                leadNumberArray_int, medianArray[i], radiusArray[i]);
            outSize += oSize;
        } else {
            floatToBytes(p, medianArray[i]);
        }

    }

    convertIntArray2ByteArray_fast_1b_args(stateArray, actualNBBlocks, r);
	
    free(stateArray);
    free(medianArray);	
    free(radiusArray);
    free(leadNumberArray_int);
    
    return outSize;
}		
					
unsigned char *
SZ_fast_compress_args_unpredictable_blocked_float(float *oriData, size_t *outSize, float absErrBound, size_t nbEle,
                                                  int blockSize) {
    float *op = oriData;

    *outSize = 0;
    size_t maxPreservedBufferSize =
            sizeof(float) * nbEle; //assume that the compressed data size would not exceed the original size
    unsigned char *outputBytes = (unsigned char *) malloc(maxPreservedBufferSize);
    memset(outputBytes, 0, maxPreservedBufferSize);
    unsigned char *leadNumberArray_int = (unsigned char *) malloc(blockSize * sizeof(int));

    size_t i = 0;
    int oSize = 0;

    size_t nbBlocks = nbEle / blockSize;
    size_t remainCount = nbEle % blockSize;
    size_t stateNBBytes =
            remainCount == 0 ? (nbBlocks % 8 == 0 ? nbBlocks / 8 : nbBlocks / 8 + 1) : ((nbBlocks + 1) % 8 == 0 ?
                                                                                        (nbBlocks + 1) / 8 :
                                                                                        (nbBlocks + 1) / 8 + 1);
    size_t actualNBBlocks = remainCount == 0 ? nbBlocks : nbBlocks + 1;

    unsigned char *stateArray = (unsigned char *) malloc(actualNBBlocks);
    float *medianArray = (float *) malloc(actualNBBlocks * sizeof(float));
    float *radiusArray = (float *) malloc(actualNBBlocks * sizeof(float));

    size_t nbConstantBlocks = computeStateMedianRadius_float(oriData, nbEle, absErrBound, blockSize, stateArray,
                                                             medianArray, radiusArray);
                                                            

    unsigned char *r = outputBytes; // + sizeof(size_t) + stateNBBytes;
    r[0] = SZx_VER_MAJOR;
    r[1] = SZx_VER_MINOR;
    r[2] = 1;
    r[3] = 0; // indicates this is not a random access version
    r[4] = (unsigned char) blockSize;
    r = r + 5; //1 byte
    sizeToBytes(r, nbConstantBlocks);
    r += sizeof(size_t); //r is the starting address of 'stateNBBytes'

    unsigned char *p = r + stateNBBytes; //p is the starting address of constant median values.
    unsigned char *q =
            p + sizeof(float) * nbConstantBlocks; //q is the starting address of the non-constant data sblocks
    //3: versions, 1: metadata: state, 1: metadata: blockSize, sizeof(size_t): nbConstantBlocks, ....
    *outSize += (3 + 1 + 1 + sizeof(size_t) + stateNBBytes + sizeof(float) * nbConstantBlocks);

    //printf("nbConstantBlocks = %zu, percent = %f\n", nbConstantBlocks, 1.0f*(nbConstantBlocks*blockSize)/nbEle);

    for (i = 0; i < nbBlocks; i++, op += blockSize) {
        if (stateArray[i]) {
            SZ_fast_compress_args_unpredictable_one_block_float(op, blockSize, absErrBound, q, &oSize,
                                                                leadNumberArray_int, medianArray[i], radiusArray[i]);
            q += oSize;
            *outSize += oSize;
        } else {
            floatToBytes(p, medianArray[i]);
            p += sizeof(float);
        }
    }

    if (remainCount != 0) {
        if (stateArray[i]) {
            SZ_fast_compress_args_unpredictable_one_block_float(op, remainCount, absErrBound, q, &oSize,
                                                                leadNumberArray_int, medianArray[i], radiusArray[i]);
            *outSize += oSize;
        } else {
            floatToBytes(p, medianArray[i]);
        }

    }

    convertIntArray2ByteArray_fast_1b_args(stateArray, actualNBBlocks, r);
	
    free(stateArray);
    free(medianArray);	
    free(radiusArray);
    free(leadNumberArray_int);

    return outputBytes;
}
					
										
void SZ_fast_compress_args_unpredictable_blocked_float2(float *oriData, size_t *outSize, unsigned char* outputBytes, float absErrBound, size_t nbEle,
                                                  int blockSize) {
    float *op = oriData;

    *outSize = 0;

    unsigned char *leadNumberArray_int = (unsigned char *) malloc(blockSize * sizeof(int));

    size_t i = 0;
    int oSize = 0;

    size_t nbBlocks = nbEle / blockSize;
    size_t remainCount = nbEle % blockSize;
    size_t stateNBBytes =
            remainCount == 0 ? (nbBlocks % 8 == 0 ? nbBlocks / 8 : nbBlocks / 8 + 1) : ((nbBlocks + 1) % 8 == 0 ?
                                                                                        (nbBlocks + 1) / 8 :
                                                                                        (nbBlocks + 1) / 8 + 1);
    size_t actualNBBlocks = remainCount == 0 ? nbBlocks : nbBlocks + 1;

    unsigned char *stateArray = (unsigned char *) malloc(actualNBBlocks);
    float *medianArray = (float *) malloc(actualNBBlocks * sizeof(float));
    float *radiusArray = (float *) malloc(actualNBBlocks * sizeof(float));

    size_t nbConstantBlocks = computeStateMedianRadius_float(oriData, nbEle, absErrBound, blockSize, stateArray,
                                                             medianArray, radiusArray);
                                                            

    unsigned char *r = outputBytes; // + sizeof(size_t) + stateNBBytes;
    r[0] = SZx_VER_MAJOR;
    r[1] = SZx_VER_MINOR;
    r[2] = 1;
    r[3] = 0; // indicates this is not a random access version
    r[4] = (unsigned char) blockSize;
    r = r + 5; //1 byte
    sizeToBytes(r, nbConstantBlocks);
    r += sizeof(size_t); //r is the starting address of 'stateNBBytes'

    unsigned char *p = r + stateNBBytes; //p is the starting address of constant median values.
    unsigned char *q =
            p + sizeof(float) * nbConstantBlocks; //q is the starting address of the non-constant data sblocks
    //3: versions, 1: metadata: state, 1: metadata: blockSize, sizeof(size_t): nbConstantBlocks, ....
    *outSize += (3 + 1 + 1 + sizeof(size_t) + stateNBBytes + sizeof(float) * nbConstantBlocks);

    //printf("nbConstantBlocks = %zu, percent = %f\n", nbConstantBlocks, 1.0f*(nbConstantBlocks*blockSize)/nbEle);

    for (i = 0; i < nbBlocks; i++, op += blockSize) {
        if (stateArray[i]) {
            SZ_fast_compress_args_unpredictable_one_block_float(op, blockSize, absErrBound, q, &oSize,
                                                                leadNumberArray_int, medianArray[i], radiusArray[i]);
            q += oSize;
            *outSize += oSize;
        } else {
            floatToBytes(p, medianArray[i]);
            p += sizeof(float);
        }
    }

    if (remainCount != 0) {
        if (stateArray[i]) {
            SZ_fast_compress_args_unpredictable_one_block_float(op, remainCount, absErrBound, q, &oSize,
                                                                leadNumberArray_int, medianArray[i], radiusArray[i]);
            *outSize += oSize;
        } else {
            floatToBytes(p, medianArray[i]);
        }

    }

    convertIntArray2ByteArray_fast_1b_args(stateArray, actualNBBlocks, r);
	
    free(stateArray);
    free(medianArray);	
    free(radiusArray);
    free(leadNumberArray_int);
}

unsigned char *
SZ_fast_compress_args_unpredictable_blocked_randomaccess_float_openmp(float *oriData, size_t *outSize, float absErrBound,
                                                               size_t nbEle, int blockSize) {
#ifdef _OPENMP
    printf("use openmp\n");

#ifdef __AVX512F__
    printf("use avx512\n");
#elif __AVX2__
    printf("use avx2\n");
#else
#endif
    printf("blockSize = %d\n",blockSize);
    sz_cost_start();
    float *op = oriData;

    size_t i = 0;
    size_t nbBlocks = nbEle / blockSize;
    size_t remainCount = nbEle % blockSize;
    size_t actualNBBlocks = remainCount == 0 ? nbBlocks : nbBlocks + 1;
    size_t stateNBBytes = (actualNBBlocks % 8 == 0 ? actualNBBlocks / 8 : actualNBBlocks / 8 + 1);

    unsigned char *stateArray = (unsigned char *) malloc(actualNBBlocks);
    float *medianArray = (float *) malloc(actualNBBlocks * sizeof(float));

    size_t nbNonConstantBlocks = 0;

    unsigned char *tmp_q = (unsigned char *) malloc(blockSize * sizeof(float) * actualNBBlocks);
    int *outSizes = (int *) malloc(actualNBBlocks * sizeof(int));
    size_t *outSizesAccumlate = (size_t *) malloc(actualNBBlocks * sizeof(size_t));
    int *nbNonConstantBlockAccumlate = (int *) malloc(actualNBBlocks * sizeof(int));

    (*outSize) = 0;
    size_t maxPreservedBufferSize =
    sizeof(float) * nbEle; //assume that the compressed data size would not exceed the original size
    unsigned char *outputBytes = (unsigned char *) malloc(maxPreservedBufferSize);
    memset(outputBytes, 0, maxPreservedBufferSize);
    unsigned char *r = outputBytes; // + sizeof(size_t) + stateNBBytes;
    r[0] = SZx_VER_MAJOR;
    r[1] = SZx_VER_MINOR;
    r[2] = 1;
    r[3] = 1; //support random access decompression
    r = r + 4; //4 byte

    int nbThreads = 1;
    unsigned char *leadNumberArray_int;
    size_t z0[200],z1[200];

    size_t nbConstantBlocks;
    unsigned char *R, *p, *q;
    float *pf;
    uint16_t *O;

#pragma omp parallel
{
#pragma omp single
{
    nbThreads = omp_get_num_threads();
    printf("nbThreads = %d\n", nbThreads);
    assert(nbThreads<200);
    leadNumberArray_int = (unsigned char *) malloc(blockSize * sizeof(int) * nbThreads);

    sz_cost_end_msg("sequential-1 malloc");
    sz_cost_start();
}
#pragma omp for reduction(+:nbNonConstantBlocks) schedule(static)
    for (i = 0; i < nbBlocks; i++) {
        float radius;
        computeStateMedianRadius_float2(op + i * blockSize, blockSize, absErrBound, stateArray + i, medianArray + i,
                                        &radius);
        if (stateArray[i]) {
            SZ_fast_compress_args_unpredictable_one_block_float(op + i * blockSize, blockSize, absErrBound,
                                                                tmp_q + i * blockSize * sizeof(float), outSizes + i,
                                                                leadNumberArray_int +
                                                                omp_get_thread_num() * blockSize * sizeof(int),
                                                                medianArray[i], radius);
            outSizesAccumlate[i]=outSizes[i];
            nbNonConstantBlocks += 1;
        }else{
            outSizes[i]=0;
            outSizesAccumlate[i]=0;
        }
    }
#pragma omp single
{
    sz_cost_end_msg("parallel-1 compress");
//    exit(0);
    if (remainCount != 0) {
        i = nbBlocks;
        float radius;
        computeStateMedianRadius_float2(op + i * blockSize, remainCount, absErrBound, stateArray + i, medianArray + i,
                                        &radius);
        if (stateArray[i]) {
            SZ_fast_compress_args_unpredictable_one_block_float(op + i * blockSize, remainCount, absErrBound,
                                                                tmp_q + i * blockSize * sizeof(float), outSizes + i,
                                                                leadNumberArray_int, medianArray[i], radius);
            outSizesAccumlate[i] = outSizes[i];
            nbNonConstantBlocks += 1;
        }else{
            outSizesAccumlate[i] = 0;
            outSizes[i]=0;
        }
    }

    nbConstantBlocks = actualNBBlocks - nbNonConstantBlocks;

    sizeToBytes(r, blockSize);
    r += sizeof(size_t);
    sizeToBytes(r, nbConstantBlocks);
    r += sizeof(size_t);
    O = (uint16_t*) r; //o is the starting address of 'block-size array'
    R = r + nbNonConstantBlocks * sizeof(uint16_t); //R is the starting address of the state array
    p = R + stateNBBytes; //p is the starting address of constant median values.
    pf = (float *) p;
    q = p + sizeof(float) * nbConstantBlocks; //q is the starting address of the non-constant data sblocks
    // unsigned char *q0 = q;
    // printf("%lu %lu %lu %lu\n",r-outputBytes, R-outputBytes, p-outputBytes, q-outputBytes);
    // 3: versions, 1: metadata: state, 1: metadata: blockSize, sizeof(size_t): nbConstantBlocks, ....
    *outSize = q - outputBytes;

    sz_cost_start();

}
    int tid = omp_get_thread_num();
    int lo = tid * actualNBBlocks / nbThreads;
    int hi = (tid + 1) * actualNBBlocks / nbThreads;
    int b;
    nbNonConstantBlockAccumlate[lo]=stateArray[lo];
    for (b = lo+1; b < hi; b++){
        outSizesAccumlate[b] = outSizesAccumlate[b] + outSizesAccumlate[b-1];
    }
    for (b = lo+1; b < hi; b++){
        nbNonConstantBlockAccumlate[b]=stateArray[b]+nbNonConstantBlockAccumlate[b-1];
    }
    z0[tid] = outSizesAccumlate[hi-1];
    z1[tid] = nbNonConstantBlockAccumlate[hi-1];
    size_t offset0=0, offset1=0;
#pragma omp barrier
    for (int j = 0; j < tid; j++) {
        offset0+=z0[j];
        offset1+=z1[j];
    }
    for (b = lo; b < hi; b++){
        outSizesAccumlate[b] = outSizesAccumlate[b] + offset0;
        nbNonConstantBlockAccumlate[b] = nbNonConstantBlockAccumlate[b] + offset1;
    }
#pragma omp single
{
    sz_cost_end_msg("parallel-2 prefix sum");
    sz_cost_start();
};
#pragma omp for schedule(static)
    for (i = 0; i < actualNBBlocks; i++) {
        if (stateArray[i]) {
            memcpy(q+outSizesAccumlate[i]-outSizes[i], tmp_q + i * blockSize * sizeof(float), outSizes[i]);
            O[nbNonConstantBlockAccumlate[i]-1]=outSizes[i];
        } else {
            pf[i-nbNonConstantBlockAccumlate[i]]=medianArray[i];
        }
    }
#pragma omp single
{
    sz_cost_end_msg("parallel-3 memcpy");
    sz_cost_start();

    *outSize += outSizesAccumlate[actualNBBlocks-1];

    convertIntArray2ByteArray_fast_1b_args(stateArray, actualNBBlocks, R);
    sz_cost_end_msg("sequential-2 int2byte");
    sz_cost_start();
    free(nbNonConstantBlockAccumlate);
    free(outSizesAccumlate);
    free(leadNumberArray_int);
    free(tmp_q);
    free(medianArray);
    free(stateArray);
    free(outSizes);
    sz_cost_end_msg("sequential-3 free");
    printf("blocksize = %d, actualNBBlocks = %lu\n", blockSize, actualNBBlocks);
    printf("nbConstantBlocks = %zu, percent = %f\n", nbConstantBlocks, 1.0f * (nbConstantBlocks * blockSize) / nbEle);
    printf("CR = %.3f, nbEle = %lu \n", nbEle*4.0/(*outSize), nbEle);
}
}
    return outputBytes;
#else
    return NULL;
#endif
}

unsigned char *
SZ_fast_compress_args_unpredictable_blocked_randomaccess_float(float *oriData, size_t *outSize, float absErrBound,
    size_t nbEle, int blockSize) {
    float *op = oriData;

    *outSize = 0;
    size_t maxPreservedBufferSize =
            sizeof(float) * nbEle; //assume that the compressed data size would not exceed the original size
    unsigned char *outputBytes = (unsigned char *) malloc(maxPreservedBufferSize);
    memset(outputBytes, 0, maxPreservedBufferSize);
    unsigned char *leadNumberArray_int = (unsigned char *) malloc(blockSize * sizeof(int));

    size_t i = 0;
    int oSize = 0;

    size_t nbBlocks = nbEle / blockSize;
    size_t remainCount = nbEle % blockSize;
    size_t actualNBBlocks = remainCount == 0 ? nbBlocks : nbBlocks + 1;

    size_t stateNBBytes = (actualNBBlocks % 8 == 0 ? actualNBBlocks / 8 : actualNBBlocks / 8 + 1);

    unsigned char *stateArray = (unsigned char *) malloc(actualNBBlocks);
    float *medianArray = (float *) malloc(actualNBBlocks * sizeof(float));
    float *radiusArray = (float *) malloc(actualNBBlocks * sizeof(float));

    size_t nbConstantBlocks = computeStateMedianRadius_float(oriData, nbEle, absErrBound, blockSize, stateArray,
                                                             medianArray, radiusArray);

    size_t nbNonConstantBlocks = actualNBBlocks - nbConstantBlocks;

    unsigned char *r = outputBytes; // + sizeof(size_t) + stateNBBytes;
    r[0] = SZx_VER_MAJOR;
    r[1] = SZx_VER_MINOR;
    r[2] = 1;
    r[3] = 1; //support random access decompression
    r = r + 4; //1 byte

    sizeToBytes(r, blockSize);
    r += sizeof(size_t);
    sizeToBytes(r, nbConstantBlocks);
    r += sizeof(size_t); //r is the starting address of 'block-size array'
    uint16_t *O=(uint16_t*)r;
    unsigned char *R = r + nbNonConstantBlocks*sizeof(uint16_t); //R is the starting address of the state array
    unsigned char *p = R + stateNBBytes; //p is the starting address of constant median values.
    unsigned char *q =
            p + sizeof(float) * nbConstantBlocks; //q is the starting address of the non-constant data sblocks
    //3: versions, 1: metadata: state, 1: metadata: blockSize, sizeof(size_t): nbConstantBlocks, ....
    *outSize = q-outputBytes;

    size_t nonConstantBlockID = 0;
    //printf("nbConstantBlocks = %zu, percent = %f\n", nbConstantBlocks, 1.0f*(nbConstantBlocks*blockSize)/nbEle);
    for (i = 0; i < nbBlocks; i++, op += blockSize) {
        if (stateArray[i]) {
            SZ_fast_compress_args_unpredictable_one_block_float(op, blockSize, absErrBound, q, &oSize,
                                                                leadNumberArray_int, medianArray[i], radiusArray[i]);
            q += oSize;
            *outSize += oSize;
            O[nonConstantBlockID++] = oSize;
        } else {
            floatToBytes(p, medianArray[i]);
            p += sizeof(float);
        }
    }

    if (remainCount != 0) {
        if (stateArray[i]) {
            SZ_fast_compress_args_unpredictable_one_block_float(op, remainCount, absErrBound, q, &oSize,
                                                                leadNumberArray_int, medianArray[i], radiusArray[i]);
            *outSize += oSize;
            O[nonConstantBlockID] = oSize;
        } else {
            floatToBytes(p, medianArray[i]);
        }

    }

    convertIntArray2ByteArray_fast_1b_args(stateArray, actualNBBlocks, R);

    free(leadNumberArray_int);

    return outputBytes;

}


unsigned char *
SZ_fast_compress_args_unpredictable_float(float *data, size_t *outSize, float absErrBound, size_t r5, size_t r4,
                                          size_t r3, size_t r2, size_t r1, float mValue, float radius) {
    size_t totalSize = 0;
    float medianValue = mValue;

    size_t dataLength = computeDataLength(r5, r4, r3, r2, r1);

    size_t maxPreservedBufferSize =
            sizeof(float) * dataLength; //assume that the compressed data size would not exceed the original size

    unsigned char *outputBytes = (unsigned char *) malloc(maxPreservedBufferSize);
    memset(outputBytes, 0, maxPreservedBufferSize);
    unsigned char *r = outputBytes; // + sizeof(size_t) + stateNBBytes;
    r[0] = SZx_VER_MAJOR;
    r[1] = SZx_VER_MINOR;
    r[2] = 1; //SZx_VER_SUPERFAST
    r[3] = 0; //support random access decompression

//	sz_cost_start();
    size_t i;
    int reqLength;
    short radExpo = getExponent_float(radius);

    computeReqLength_float(absErrBound, radExpo, &reqLength, &medianValue);

    int reqBytesLength = reqLength / 8;
    int resiBitsLength = reqLength % 8;
    int rightShiftBits = 0;

    size_t leadNumberArray_size = dataLength % 4 == 0 ? dataLength / 4 : dataLength / 4 + 1;

    register lfloat lfBuf_pre;
    register lfloat lfBuf_cur;
    lfBuf_pre.ivalue = 0;

    unsigned char *leadNumberArray = outputBytes + 4 + 1 + sizeof(float) + sizeof(size_t);

    unsigned char *exactMidbyteArray = leadNumberArray + leadNumberArray_size;

    if (resiBitsLength != 0) {
        rightShiftBits = 8 - resiBitsLength;
        reqBytesLength++;
    }

    register unsigned char leadingNum = 0;

    unsigned char *leadNumberArray_int = (unsigned char *) malloc(dataLength);

    size_t residualMidBytes_size = 0;
    if (sysEndianType == LITTLE_ENDIAN_SYSTEM) {
        if (reqBytesLength == 3) {
            for (i = 0; i < dataLength; i++) {
                leadingNum = 0;
                lfBuf_cur.value = data[i] - medianValue;

                lfBuf_cur.ivalue = lfBuf_cur.ivalue >> rightShiftBits;

                lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;

                if (lfBuf_pre.ivalue >> 8 == 0)
                    leadingNum = 3;
                else if (lfBuf_pre.ivalue >> 16 == 0)
                    leadingNum = 2;
                else if (lfBuf_pre.ivalue >> 24 == 0)
                    leadingNum = 1;

                leadNumberArray_int[i] = leadingNum;

                if (leadingNum == 0) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[1];
                    exactMidbyteArray[residualMidBytes_size + 1] = lfBuf_cur.byte[2];
                    exactMidbyteArray[residualMidBytes_size + 2] = lfBuf_cur.byte[3];
                    residualMidBytes_size += 3;
                } else if (leadingNum == 1) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[1];
                    exactMidbyteArray[residualMidBytes_size + 1] = lfBuf_cur.byte[2];
                    residualMidBytes_size += 2;
                } else if (leadingNum == 2) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[1];
                    residualMidBytes_size++;
                }

                lfBuf_pre = lfBuf_cur;
            }
        } else if (reqBytesLength == 2) {
            for (i = 0; i < dataLength; i++) {

                leadingNum = 0;
                lfBuf_cur.value = data[i] - medianValue;

                lfBuf_cur.ivalue = lfBuf_cur.ivalue >> rightShiftBits;

                lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;

                if (lfBuf_pre.ivalue >> 8 == 0)
                    leadingNum = 3;
                else if (lfBuf_pre.ivalue >> 16 == 0)
                    leadingNum = 2;
                else if (lfBuf_pre.ivalue >> 24 == 0)
                    leadingNum = 1;

                leadNumberArray_int[i] = leadingNum;

                if (leadingNum == 0) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[2];
                    exactMidbyteArray[residualMidBytes_size + 1] = lfBuf_cur.byte[3];
                    residualMidBytes_size += 2;
                } else if (leadingNum == 1) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[2];
                    residualMidBytes_size++;
                }

                lfBuf_pre = lfBuf_cur;
            }
        } else if (reqBytesLength == 1) {
            for (i = 0; i < dataLength; i++) {
                leadingNum = 0;
                lfBuf_cur.value = data[i] - medianValue;

                lfBuf_cur.ivalue = lfBuf_cur.ivalue >> rightShiftBits;

                lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;

                if (lfBuf_pre.ivalue >> 8 == 0)
                    leadingNum = 3;
                else if (lfBuf_pre.ivalue >> 16 == 0)
                    leadingNum = 2;
                else if (lfBuf_pre.ivalue >> 24 == 0)
                    leadingNum = 1;

                leadNumberArray_int[i] = leadingNum;

                if (leadingNum == 0) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[3];
                    residualMidBytes_size++;
                }

                lfBuf_pre = lfBuf_cur;
            }
        }else //reqBytesLength == 4
        {
            for (i = 0; i < dataLength; i++) {
                leadingNum = 0;
                lfBuf_cur.value = data[i] - medianValue;

                lfBuf_cur.ivalue = lfBuf_cur.ivalue >> rightShiftBits;

                lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;

                if (lfBuf_pre.ivalue >> 8 == 0)
                    leadingNum = 3;
                else if (lfBuf_pre.ivalue >> 16 == 0)
                    leadingNum = 2;
                else if (lfBuf_pre.ivalue >> 24 == 0)
                    leadingNum = 1;

                leadNumberArray_int[i] = leadingNum;

                if (leadingNum == 0) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[0];
                    exactMidbyteArray[residualMidBytes_size + 1] = lfBuf_cur.byte[1];
                    exactMidbyteArray[residualMidBytes_size + 2] = lfBuf_cur.byte[2];
                    exactMidbyteArray[residualMidBytes_size + 3] = lfBuf_cur.byte[3];
                    residualMidBytes_size += 4;
                } else if (leadingNum == 1) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[0];
                    exactMidbyteArray[residualMidBytes_size + 1] = lfBuf_cur.byte[1];
                    exactMidbyteArray[residualMidBytes_size + 2] = lfBuf_cur.byte[2];
                    residualMidBytes_size += 3;
                } else if (leadingNum == 2) {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[0];
                    exactMidbyteArray[residualMidBytes_size + 1] = lfBuf_cur.byte[1];
                    residualMidBytes_size += 2;
                } else //leadingNum == 3
                {
                    exactMidbyteArray[residualMidBytes_size] = lfBuf_cur.byte[0];
                    residualMidBytes_size++;
                }

                lfBuf_pre = lfBuf_cur;
            }
        }

        convertIntArray2ByteArray_fast_2b_args(leadNumberArray_int, dataLength, leadNumberArray);

        int k = 4;

        unsigned char reqLengthB = (unsigned char) reqLength;
        outputBytes[k] = reqLengthB;
        k++;
        floatToBytes(&(outputBytes[k]), medianValue);
        k += sizeof(float);
        sizeToBytes(&(outputBytes[k]), leadNumberArray_size);

        totalSize = 4 + 1 + sizeof(float) + sizeof(size_t) + leadNumberArray_size + residualMidBytes_size;
    } else {

    }

    *outSize = totalSize;

    free(leadNumberArray_int);
//	sz_cost_end();
//	printf("compression time = %f\n", sz_totalCost);

    return outputBytes;
}

unsigned char *SZ_skip_compress_float(float *data, size_t dataLength, size_t *outSize) {
    *outSize = dataLength * sizeof(float);
    unsigned char *out = (unsigned char *) malloc(dataLength * sizeof(float));
    memcpy(out, data, dataLength * sizeof(float));
    return out;
}

inline void computeReqLength_float(double realPrecision, short radExpo, int *reqLength, float *medianValue) {
    short reqExpo = getPrecisionReqLength_double(realPrecision);
    *reqLength = 9 + radExpo - reqExpo + 1; //radExpo-reqExpo == reqMantiLength
    if (*reqLength < 9)
        *reqLength = 9;
    if (*reqLength > 32) {
        *reqLength = 32;
        *medianValue = 0;
    }
}
