/**
 *  @file szx.c
 *  @author Sheng Di
 *  @date Jan, 2022
 *  @brief 
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "szx.h"
#include "szx_rw.h"
#include <szx_float.h>
#include <szx_double.h>
#include <szxd_float.h>
#include <szxd_double.h>
#include <szx_defines.h>
#include <szx_globals.h>

using namespace szx;

namespace szx{

int sumReqNbBits = 0;
int sumReqNbBytes = 0;
int sum_leadNumberArray_size = 0;
int sum_residualMidBytes_size = 0;
int sum_actual_leadNumbers = 0;

int versionNumber[4] = {SZx_VER_MAJOR,SZx_VER_MINOR,SZx_VER_BUILD,SZx_VER_REVISION};

int dataEndianType = LITTLE_ENDIAN_DATA; //*endian type of the data read from disk
int sysEndianType = LITTLE_ENDIAN_SYSTEM; //*sysEndianType is actually set automatically.

}

int SZx_computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	//printf("hello\n");
	int dimension;
	if(r1==0)
	{
		dimension = 0;
	}
	else if(r2==0)
	{
		dimension = 1;
	}
	else if(r3==0)
	{
		dimension = 2;
		//printf("hello\n");
	}
	else if(r4==0)
	{
		dimension = 3;
	}
	else if(r5==0)
	{
		dimension = 4;
	}
	else
	{
		dimension = 5;
	}
	return dimension;
}

size_t SZx_computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t dataLength;
	if(r1==0)
	{
		dataLength = 0;
	}
	else if(r2==0)
	{
		dataLength = r1;
	}
	else if(r3==0)
	{
		dataLength = r1*r2;
	}
	else if(r4==0)
	{
		dataLength = r1*r2*r3;
	}
	else if(r5==0)
	{
		dataLength = r1*r2*r3*r4;
	}
	else
	{
		dataLength = r1*r2*r3*r4*r5;
	}
	return dataLength;
}

/**
 * @brief		check dimension and correct it if needed
 * @return 	0 (didn't change dimension)
 * 					1 (dimension is changed)
 * 					2 (dimension is problematic)
 **/
int SZx_filterDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* correctedDimension)
{
	int dimensionCorrected = 0;
	int dim = SZx_computeDimension(r5, r4, r3, r2, r1);
	correctedDimension[0] = r1;
	correctedDimension[1] = r2;
	correctedDimension[2] = r3;
	correctedDimension[3] = r4;
	correctedDimension[4] = r5;
	size_t* c = correctedDimension;
	if(dim==1)
	{
		if(r1<1)
			return 2;
	}
	else if(dim==2)
	{
		if(r2==1)
		{
			c[1]= 0;
			dimensionCorrected = 1;
		}	
		if(r1==1) //remove this dimension
		{
			c[0] = c[1]; 
			c[1] = c[2];
			dimensionCorrected = 1;
		}
	}
	else if(dim==3)
	{
		if(r3==1)
		{
			c[2] = 0;
			dimensionCorrected = 1;
		}	
		if(r2==1)
		{
			c[1] = c[2];
			c[2] = c[3];
			dimensionCorrected = 1;
		}
		if(r1==1)
		{
			c[0] = c[1];
			c[1] = c[2];
			c[2] = c[3];
			dimensionCorrected = 1;
		}
	}
	else if(dim==4)
	{
		if(r4==1)
		{
			c[3] = 0;
			dimensionCorrected = 1;
		}
		if(r3==1)
		{
			c[2] = c[3];
			c[3] = c[4];
			dimensionCorrected = 1;
		}
		if(r2==1)
		{
			c[1] = c[2];
			c[2] = c[3];
			c[3] = c[4];
			dimensionCorrected = 1;
		}
		if(r1==1)
		{
			c[0] = c[1];
			c[1] = c[2];
			c[2] = c[3];
			c[3] = c[4];
			dimensionCorrected = 1;
		}
	}
	else if(dim==5)
	{
		if(r5==1)
		{
			c[4] = 0;
			dimensionCorrected = 1;
		}
		if(r4==1)
		{
			c[3] = c[4];
			c[4] = 0;
			dimensionCorrected = 1;
		}
		if(r3==1)
		{
			c[2] = c[3];
			c[3] = c[4];
			c[4] = 0;
			dimensionCorrected = 1;
		}
		if(r2==1)
		{
			c[1] = c[2];
			c[2] = c[3];
			c[3] = c[4];
			c[4] = 0;
			dimensionCorrected = 1;
		}
		if(r1==1)
		{
			c[0] = c[1];
			c[1] = c[2];
			c[2] = c[3];
			c[3] = c[4];
			c[4] = 0;
			dimensionCorrected = 1;
		}
	}
	
	return dimensionCorrected;
	
}

unsigned char* SZx_fast_compress_args(int fastMode, int dataType, void *data, size_t *outSize, int errBoundMode, float absErrBound,
float relBoundRatio, float compressionRatio, float tolerance, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	int blockSize = 128;	
	unsigned char*  bytes = NULL;
	size_t length = SZx_computeDataLength(r5, r4, r3, r2, r1);
	
	if(dataType == SZx_FLOAT)
	{
		if(errBoundMode==FXR)
		{
			bytes = SZx_fast_compress_args_unpredictable_blocked_fixed_rate_float((float*)data, outSize, compressionRatio, tolerance, length, blockSize); 			
		}
		else //not fix rate mode
		{
			if(fastMode == SZx_WITH_BLOCK_FAST_CMPR || fastMode == SZx_RANDOMACCESS_FAST_CMPR || fastMode == SZx_OPENMP_FAST_CMPR)
			{
				float realPrecision = absErrBound;
				if(errBoundMode==REL)
				{
					float valueRange = computeValueRange_float((float*)data, length, NULL, NULL);
					realPrecision = valueRange*relBoundRatio;
				}

				if (fastMode == SZx_RANDOMACCESS_FAST_CMPR) {
					bytes = SZx_fast_compress_args_unpredictable_blocked_randomaccess_float((float*)data, outSize, realPrecision, length, blockSize);
				} 
				else if(fastMode == SZx_OPENMP_FAST_CMPR)
				{
					#ifdef _OPENMP
					bytes = SZx_fast_compress_args_unpredictable_blocked_randomaccess_float_openmp((float*)data, outSize, realPrecision, length,
																								  blockSize);
					#else
					bytes = SZx_fast_compress_args_unpredictable_blocked_randomaccess_float((float*)data, outSize, realPrecision, length, blockSize);
					printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the compression is degraded to serial version automatically.\n");
					#endif
				}
				else {
					bytes = SZx_fast_compress_args_unpredictable_blocked_float((float*)data, outSize, realPrecision, length, blockSize);
				}
				return bytes;
			}
			else
			{
				//compute value range
				
				float radius = 0;
				float medianValue = 0;
				float valueRange = computeValueRange_float((float*)data, length, &radius, &medianValue);

				float realPrecision = 0;
				if(errBoundMode==ABS)
					realPrecision = absErrBound;
				else if(errBoundMode==REL)
					realPrecision = valueRange*relBoundRatio;

				bytes = SZx_fast_compress_args_unpredictable_float((float*)data, outSize, realPrecision, r5, r4, r3, r2, r1, medianValue, radius);		
			}			
		}
		

	}
	else if(dataType == SZx_DOUBLE)
	{
		if(errBoundMode==FXR)
		{
			bytes = SZx_fast_compress_args_unpredictable_blocked_fixed_rate_double((double*)data, outSize, compressionRatio, tolerance, length, blockSize); 			
		}
		else //not fix rate mode
		{
			if(fastMode == SZx_WITH_BLOCK_FAST_CMPR || fastMode == SZx_RANDOMACCESS_FAST_CMPR || fastMode == SZx_OPENMP_FAST_CMPR)
			{
				float realPrecision = absErrBound;
				if(errBoundMode==REL)
				{
					double valueRange = computeValueRange_double((double*)data, length, NULL, NULL);
					realPrecision = valueRange*relBoundRatio;
				}

				int blockSize = 128;
				if (fastMode == SZx_RANDOMACCESS_FAST_CMPR) {
					bytes = SZx_fast_compress_args_unpredictable_blocked_randomaccess_double((double*)data, outSize, realPrecision, length, blockSize);
				} 
				else if(fastMode == SZx_OPENMP_FAST_CMPR)
				{
					#ifdef _OPENMP
					bytes = SZx_fast_compress_args_unpredictable_blocked_randomaccess_double_openmp((double*)data, outSize, realPrecision, length,
																								blockSize);
					#else
					bytes = SZx_fast_compress_args_unpredictable_blocked_randomaccess_double((double*)data, outSize, realPrecision, length, blockSize);
					printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the compression is degraded to serial version automatically.\n");
					#endif
				}
				else {
					bytes = SZx_fast_compress_args_unpredictable_blocked_double((double*)data, outSize, realPrecision, length, blockSize);
				}
				return bytes;
			}
			else
			{
				//compute value range
				float radius = 0;
				float medianValue = 0;
				double valueRange = computeValueRange_double((double*)data, length, &radius, &medianValue);

				float realPrecision = 0;
				if(errBoundMode==ABS)
					realPrecision = absErrBound;
				else if(errBoundMode==REL)
					realPrecision = valueRange*relBoundRatio;

				bytes = SZx_fast_compress_args_unpredictable_double((double*)data, outSize, realPrecision, r5, r4, r3, r2, r1, medianValue, radius);		
			}
		}		
	}
    return bytes;

}

/**
 * @deprecated
 * */
void* SZx_fast_decompress_pred(int dataType, float* preData, unsigned char *curBytes, size_t byteLength, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    int x = 1;
    char *y = (char*)&x;
    if(*y==1)
        sysEndianType = LITTLE_ENDIAN_SYSTEM;
    else //=0
        sysEndianType = BIG_ENDIAN_SYSTEM;

    if(dataType == SZx_FLOAT)
    {
        float* newFloatData = NULL;
        SZx_fast_decompress_args_with_prediction_float(&newFloatData, preData, r5, r4, r3, r2, r1, curBytes, byteLength);
        return newFloatData;
    }
    else if(dataType == SZx_DOUBLE)
    {
        double* newDoubleData = NULL;
        //SZx_fast_decompress_args_unpredictable_float(&newDoubleData, r5, r4, r3, r2, r1, bytes, byteLength, 0, NULL);
        return newDoubleData;
    }

    return NULL;
}

void* SZx_fast_decompress(int fastMode, int dataType, unsigned char *bytes, size_t byteLength, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t nbEle = SZx_computeDataLength(r5, r4, r3, r2, r1);
    int x = 1;
    char *y = (char*)&x;
    if(*y==1)
        sysEndianType = LITTLE_ENDIAN_SYSTEM;
    else //=0
        sysEndianType = BIG_ENDIAN_SYSTEM;

    if(dataType == SZx_FLOAT)
    {
        float* newFloatData = NULL;
        if(fastMode == SZx_NO_BLOCK_FAST_CMPR)
            SZx_fast_decompress_args_unpredictable_float(&newFloatData, r5, r4, r3, r2, r1, bytes, byteLength);
		else if(fastMode == SZx_WITH_BLOCK_FAST_CMPR)
			SZx_fast_decompress_args_unpredictable_blocked_float(&newFloatData, nbEle, bytes);            
        else if(fastMode == SZx_RANDOMACCESS_FAST_CMPR)
			newFloatData = (float*)SZx_fast_decompress_args_unpredictable_blocked_randomaccess_float(nbEle, bytes);
        else //SZx_openmp
        {
#ifdef _OPENMP
                newFloatData = (float*)SZx_fast_decompress_args_unpredictable_blocked_randomaccess_float_openmp(nbEle, bytes);
#else
                SZx_fast_decompress_args_unpredictable_blocked_float(&newFloatData, nbEle, bytes);
                printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the decompression is degraded to serial version automatically.\n");
#endif
        }
        return newFloatData;
    }
    else if(dataType == SZx_DOUBLE)
    {
        double* newFloatData = NULL;
        if(fastMode == SZx_NO_BLOCK_FAST_CMPR)
            SZx_fast_decompress_args_unpredictable_double(&newFloatData, r5, r4, r3, r2, r1, bytes, byteLength);
		else if(fastMode == SZx_WITH_BLOCK_FAST_CMPR)
			SZx_fast_decompress_args_unpredictable_blocked_double(&newFloatData, nbEle, bytes);            
        else if(fastMode == SZx_RANDOMACCESS_FAST_CMPR)
			newFloatData = (double*)SZx_fast_decompress_args_unpredictable_blocked_randomaccess_double(nbEle, bytes);
        else //SZx_openmp
        {
#ifdef _OPENMP
                newFloatData = SZx_fast_decompress_args_unpredictable_blocked_randomaccess_double_openmp(nbEle, bytes);
#else
                SZx_fast_decompress_args_unpredictable_blocked_double(&newFloatData, nbEle, bytes);
                printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the decompression is degraded to serial version automatically.\n");
#endif
        }
        return newFloatData;
    }

    return NULL;
}

double SZx_computeValueRange(void* oriData, int dataType, size_t length, float* radius, float* medianValue)
{
	if(dataType == SZx_FLOAT)
		return computeValueRange_float((float*)oriData, length, radius, medianValue);
	else //SZx_DOUBLE
		return computeValueRange_double((double*)oriData, length, radius, medianValue);
}


unsigned char *
SZx_fast_compress_args_unpredictable_blocked_randomaccess_openmp(void *oriData, int dataType, size_t *outSize, float absErrBound,
                                                               size_t nbEle, int blockSize) {
	if(dataType == SZx_FLOAT)
		return SZx_fast_compress_args_unpredictable_blocked_randomaccess_float_openmp((float*)oriData, outSize, absErrBound, nbEle, blockSize);
	else //SZx_DOUBLE
		return SZx_fast_compress_args_unpredictable_blocked_randomaccess_double_openmp((double*)oriData, outSize, absErrBound, nbEle, blockSize);
}

void* SZx_fast_decompress_args_unpredictable_blocked_randomaccess_openmp(int dataType, size_t nbEle, unsigned char* cmpBytes)
{
	if(dataType == SZx_FLOAT)
		return SZx_fast_decompress_args_unpredictable_blocked_randomaccess_float_openmp(nbEle, cmpBytes);
	else //SZx_DOUBLE
		return SZx_fast_decompress_args_unpredictable_blocked_randomaccess_double_openmp(nbEle, cmpBytes);
}

float SZx_estimateErrorBoundbasedonCR(int dataType, float targetCompressionRatio, float tolerance, void* data, float initErrorBound, int blockSize, size_t nbEle)
{
	if(dataType == SZx_FLOAT)
		return estimateErrorBoundbasedonCR_float(targetCompressionRatio, tolerance, (float*)data, initErrorBound, blockSize, nbEle);
	else//dataType==SZx_DOUBLE
		return estimateErrorBoundbasedonCR_double(targetCompressionRatio, tolerance, (double*)data, initErrorBound, blockSize, nbEle);
}

unsigned char *
SZx_fast_compress_args_unpredictable_blocked(void *oriData, int dataType, size_t *outSize, float absErrBound, size_t nbEle, int blockSize) 
{
	if(dataType == SZx_FLOAT)
		return SZx_fast_compress_args_unpredictable_blocked_float((float*)oriData, outSize, absErrBound, nbEle, blockSize);
	else //dataType == SZx_DOUBLE
		return SZx_fast_compress_args_unpredictable_blocked_double((double*)oriData, outSize, absErrBound, nbEle, blockSize);
}

float computeRadiusBuffer(void* oriData, int dataType, size_t nbEle, int samplingRate, int blockSize, float** radiusArray, float** mediusArray, void* buffer)
{
	if(dataType==SZx_FLOAT)
	{
		return computeRadiusBuffer_float((float*)oriData, nbEle, samplingRate, blockSize, radiusArray, mediusArray, (float**)buffer);
	}
	else //SZx_DOUBLE
	{
		return computeRadiusBuffer_double((double*)oriData, nbEle, samplingRate, blockSize, radiusArray, mediusArray, (double**)buffer);
	}
}

float estimateCRbasedonErrorBound_buffered(int dataType, float errorBound, void* buffer, float* medianArray, float* radiusArray, int samplingRate, int blockSize, size_t nbEle, 
size_t *sumReqNbBytes, size_t *sum_actual_leadNumbers)
{
	if(dataType==SZx_FLOAT)
	{
		return estimateCRbasedonErrorBound_buffered_float(errorBound, (float*)buffer, medianArray, radiusArray, samplingRate, blockSize, nbEle, sumReqNbBytes, sum_actual_leadNumbers);
	}
	else
	{
		return estimateCRbasedonErrorBound_buffered_double(errorBound, (double*)buffer, medianArray, radiusArray, samplingRate, blockSize, nbEle, sumReqNbBytes, sum_actual_leadNumbers);
	}
}


size_t SZx_fast_compress_args_unpredictable_blocked_args(void *oriData, int dataType, unsigned char* outputBytes, float absErrBound, size_t nbEle, int blockSize)
{
	if(dataType==SZx_FLOAT)
	{
		return SZx_fast_compress_args_unpredictable_blocked_args_float((float*)oriData, outputBytes, absErrBound, nbEle, blockSize);
	}
	else
	{
		return 0;  //not implemented yet
		//return SZx_fast_compress_args_unpredictable_blocked_args_double((double*)oriData, outputBytes, absErrBound, nbEle, blockSize);
	}
}

void SZx_fast_decompress_args_unpredictable_blocked(void* newData, int dataType, size_t nbEle, unsigned char* cmpBytes)
{
	if(dataType==SZx_FLOAT)
	{
		SZx_fast_decompress_args_unpredictable_blocked_float((float**)newData, nbEle, cmpBytes);
	}
	else
	{
		//not implemented yet?
	}
}



