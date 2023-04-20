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

int sumReqNbBits = 0;
int sumReqNbBytes = 0;
int sum_leadNumberArray_size = 0;
int sum_residualMidBytes_size = 0;
int sum_actual_leadNumbers = 0;

int versionNumber[4] = {SZx_VER_MAJOR,SZx_VER_MINOR,SZx_VER_BUILD,SZx_VER_REVISION};

int dataEndianType = LITTLE_ENDIAN_DATA; //*endian type of the data read from disk
int sysEndianType = LITTLE_ENDIAN_SYSTEM; //*sysEndianType is actually set automatically.

int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
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

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
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
int filterDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* correctedDimension)
{
	int dimensionCorrected = 0;
	int dim = computeDimension(r5, r4, r3, r2, r1);
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

unsigned char* SZ_fast_compress_args(int fastMode, int dataType, void *data, size_t *outSize, int errBoundMode, float absErrBound,
float relBoundRatio, float compressionRatio, float tolerance, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	int blockSize = 128;	
	unsigned char*  bytes = NULL;
	size_t length = computeDataLength(r5, r4, r3, r2, r1);
	
	if(dataType == SZ_FLOAT)
	{
		if(errBoundMode==FXR)
		{
			bytes = SZ_fast_compress_args_unpredictable_blocked_fixed_rate_float(data, outSize, compressionRatio, tolerance, length, blockSize); 			
		}
		else //not fix rate mode
		{
			if(fastMode == SZx_WITH_BLOCK_FAST_CMPR || fastMode == SZx_RANDOMACCESS_FAST_CMPR || fastMode == SZx_OPENMP_FAST_CMPR)
			{
				float realPrecision = absErrBound;
				if(errBoundMode==REL)
				{
					float valueRange = computeValueRange_float(data, length, NULL, NULL);
					realPrecision = valueRange*relBoundRatio;
				}

				if (fastMode == SZx_RANDOMACCESS_FAST_CMPR) {
					bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_float(data, outSize, realPrecision, length, blockSize);
				} 
				else if(fastMode == SZx_OPENMP_FAST_CMPR)
				{
					#ifdef _OPENMP
					bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_float_openmp(data, outSize, realPrecision, length,
																								  blockSize);
					#else
					bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_float(data, outSize, realPrecision, length, blockSize);
					printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the compression is degraded to serial version automatically.\n");
					#endif
				}
				else {
					bytes = SZ_fast_compress_args_unpredictable_blocked_float(data, outSize, realPrecision, length, blockSize);
				}
				return bytes;
			}
			else
			{
				//compute value range
				
				float radius = 0;
				float medianValue = 0;
				float valueRange = computeValueRange_float(data, length, &radius, &medianValue);

				float realPrecision = 0;
				if(errBoundMode==ABS)
					realPrecision = absErrBound;
				else if(errBoundMode==REL)
					realPrecision = valueRange*relBoundRatio;

				bytes = SZ_fast_compress_args_unpredictable_float(data, outSize, realPrecision, r5, r4, r3, r2, r1, medianValue, radius);		
			}			
		}
		

	}
	else if(dataType == SZ_DOUBLE)
	{
		if(fastMode == SZx_WITH_BLOCK_FAST_CMPR || fastMode == SZx_RANDOMACCESS_FAST_CMPR || fastMode == SZx_OPENMP_FAST_CMPR)
		{
			float realPrecision = absErrBound;
			if(errBoundMode==REL)
			{
				double valueRange = computeValueRange_double(data, length, NULL, NULL);
				realPrecision = valueRange*relBoundRatio;
			}

			int blockSize = 128;
			if (fastMode == SZx_RANDOMACCESS_FAST_CMPR) {
				bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_double(data, outSize, realPrecision, length, blockSize);
			} 
			else if(fastMode == SZx_OPENMP_FAST_CMPR)
			{
				#ifdef _OPENMP
				bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_double_openmp(data, outSize, realPrecision, length,
																							  blockSize);
				#else
				bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_double(data, outSize, realPrecision, length, blockSize);
				printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the compression is degraded to serial version automatically.\n");
				#endif
			}
			else {
				bytes = SZ_fast_compress_args_unpredictable_blocked_double(data, outSize, realPrecision, length, blockSize);
			}
			return bytes;
		}
		else
		{
			//compute value range
			float radius = 0;
			float medianValue = 0;
			double valueRange = computeValueRange_double(data, length, &radius, &medianValue);

			float realPrecision = 0;
			if(errBoundMode==ABS)
				realPrecision = absErrBound;
			else if(errBoundMode==REL)
				realPrecision = valueRange*relBoundRatio;

			bytes = SZ_fast_compress_args_unpredictable_double(data, outSize, realPrecision, r5, r4, r3, r2, r1, medianValue, radius);		
		}		
	}

    return bytes;

}

unsigned char* SZ_fast_compress_args_test(int fastMode, int dataType, void *data, size_t *outSize, int errBoundMode, float absErrBound,
float relBoundRatio, float compressionRatio, float tolerance, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	int blockSize = 128;	
	unsigned char*  bytes = NULL;
	size_t length = computeDataLength(r5, r4, r3, r2, r1);
	
	if(dataType == SZ_FLOAT)
	{
		if(errBoundMode==FXR)
		{
			bytes = SZ_fast_compress_args_unpredictable_blocked_fixed_rate_float(data, outSize, compressionRatio, tolerance, length, blockSize); 			
		}
		else //not fix rate mode
		{
			if(fastMode == SZx_WITH_BLOCK_FAST_CMPR || fastMode == SZx_RANDOMACCESS_FAST_CMPR || fastMode == SZx_OPENMP_FAST_CMPR)
			{
				float realPrecision = absErrBound;
				if(errBoundMode==REL)
				{
					float valueRange = computeValueRange_float(data, length, NULL, NULL);
					realPrecision = valueRange*relBoundRatio;
				}

				if (fastMode == SZx_RANDOMACCESS_FAST_CMPR) {
					bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_float(data, outSize, realPrecision, length, blockSize);
				} 
				else if(fastMode == SZx_OPENMP_FAST_CMPR)
				{
					#ifdef _OPENMP
					bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_float_openmp_test(data, outSize, realPrecision, length,
																								  blockSize);
					#else
					bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_float(data, outSize, realPrecision, length, blockSize);
					printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the compression is degraded to serial version automatically.\n");
					#endif
				}
				else {
					bytes = SZ_fast_compress_args_unpredictable_blocked_float_test(data, outSize, realPrecision, length, blockSize);
				}
				return bytes;
			}
			else
			{
				//compute value range
				
				float radius = 0;
				float medianValue = 0;
				float valueRange = computeValueRange_float(data, length, &radius, &medianValue);

				float realPrecision = 0;
				if(errBoundMode==ABS)
					realPrecision = absErrBound;
				else if(errBoundMode==REL)
					realPrecision = valueRange*relBoundRatio;

				bytes = SZ_fast_compress_args_unpredictable_float(data, outSize, realPrecision, r5, r4, r3, r2, r1, medianValue, radius);		
			}			
		}
		

	}
	else if(dataType == SZ_DOUBLE)
	{
		if(fastMode == SZx_WITH_BLOCK_FAST_CMPR || fastMode == SZx_RANDOMACCESS_FAST_CMPR || fastMode == SZx_OPENMP_FAST_CMPR)
		{
			float realPrecision = absErrBound;
			if(errBoundMode==REL)
			{
				double valueRange = computeValueRange_double(data, length, NULL, NULL);
				realPrecision = valueRange*relBoundRatio;
			}

			int blockSize = 128;
			if (fastMode == SZx_RANDOMACCESS_FAST_CMPR) {
				bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_double(data, outSize, realPrecision, length, blockSize);
			} 
			else if(fastMode == SZx_OPENMP_FAST_CMPR)
			{
				#ifdef _OPENMP
				bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_double_openmp(data, outSize, realPrecision, length,
																							  blockSize);
				#else
				bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_double(data, outSize, realPrecision, length, blockSize);
				printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the compression is degraded to serial version automatically.\n");
				#endif
			}
			else {
				bytes = SZ_fast_compress_args_unpredictable_blocked_double(data, outSize, realPrecision, length, blockSize);
			}
			return bytes;
		}
		else
		{
			//compute value range
			float radius = 0;
			float medianValue = 0;
			double valueRange = computeValueRange_double(data, length, &radius, &medianValue);

			float realPrecision = 0;
			if(errBoundMode==ABS)
				realPrecision = absErrBound;
			else if(errBoundMode==REL)
				realPrecision = valueRange*relBoundRatio;

			bytes = SZ_fast_compress_args_unpredictable_double(data, outSize, realPrecision, r5, r4, r3, r2, r1, medianValue, radius);		
		}		
	}

    return bytes;

}

unsigned char* SZ_fast_compress_args2(int fastMode, int dataType, void *data, size_t *outSize, unsigned char* outputBytes, int errBoundMode, float absErrBound,
float relBoundRatio, float compressionRatio, float tolerance, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	int blockSize = 128;
	unsigned char*  bytes = NULL;	
	size_t length = computeDataLength(r5, r4, r3, r2, r1);
	
	if(dataType == SZ_FLOAT)
	{
		if(errBoundMode==FXR)
		{
			if(fastMode == SZx_OPENMP_FAST_CMPR)
			{
				SZ_fast_compress_args_unpredictable_blocked_fixed_rate_float2_openmp(data, outSize, outputBytes, compressionRatio, tolerance, length, blockSize);
			}
			else
			{
				SZ_fast_compress_args_unpredictable_blocked_fixed_rate_float2(data, outSize, outputBytes, compressionRatio, tolerance, length, blockSize); 	
			}
		}
		else //not fix rate mode
		{
			if(fastMode == SZx_WITH_BLOCK_FAST_CMPR || fastMode == SZx_RANDOMACCESS_FAST_CMPR || fastMode == SZx_OPENMP_FAST_CMPR)
			{
				float realPrecision = absErrBound;
				if(errBoundMode==REL)
				{
					float valueRange = computeValueRange_float(data, length, NULL, NULL);
					realPrecision = valueRange*relBoundRatio;
				}

				if (fastMode == SZx_RANDOMACCESS_FAST_CMPR) {
					SZ_fast_compress_args_unpredictable_blocked_float2(data, outSize, outputBytes, realPrecision, length, blockSize);
				} 
				else if(fastMode == SZx_OPENMP_FAST_CMPR)
				{
					#ifdef _OPENMP
					SZ_fast_compress_args_unpredictable_blocked_randomaccess_float2_openmp(data, outSize, outputBytes, realPrecision, length,
																								  blockSize);
					#else
					SZ_fast_compress_args_unpredictable_blocked_float2(data, outSize, outputBytes, realPrecision, length, blockSize);
					printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the compression is degraded to serial version automatically.\n");
					#endif
				}
				else {
					SZ_fast_compress_args_unpredictable_blocked_float2(data, outSize, outputBytes, realPrecision, length, blockSize);
				}
			}
			else
			{
				//compute value range
				
				float radius = 0;
				float medianValue = 0;
				float valueRange = computeValueRange_float(data, length, &radius, &medianValue);

				float realPrecision = 0;
				if(errBoundMode==ABS)
					realPrecision = absErrBound;
				else if(errBoundMode==REL)
					realPrecision = valueRange*relBoundRatio;

				// SZ_fast_compress_args_unpredictable_float2(data, outSize, outputBytes, realPrecision, r5, r4, r3, r2, r1, medianValue, radius);		
			}			
		}
		

	}
	else if(dataType == SZ_DOUBLE)
	{
		if(fastMode == SZx_WITH_BLOCK_FAST_CMPR || fastMode == SZx_RANDOMACCESS_FAST_CMPR || fastMode == SZx_OPENMP_FAST_CMPR)
		{
			float realPrecision = absErrBound;
			if(errBoundMode==REL)
			{
				double valueRange = computeValueRange_double(data, length, NULL, NULL);
				realPrecision = valueRange*relBoundRatio;
			}

			int blockSize = 128;
			if (fastMode == SZx_RANDOMACCESS_FAST_CMPR) {
				bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_double(data, outSize, realPrecision, length, blockSize);
			} 
			else if(fastMode == SZx_OPENMP_FAST_CMPR)
			{
				#ifdef _OPENMP
				bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_double_openmp(data, outSize, realPrecision, length,
																							  blockSize);
				#else
				bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_double(data, outSize, realPrecision, length, blockSize);
				printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the compression is degraded to serial version automatically.\n");
				#endif
			}
			else {
				bytes = SZ_fast_compress_args_unpredictable_blocked_double(data, outSize, realPrecision, length, blockSize);
			}
			return bytes;
		}
		else
		{
			//compute value range
			float radius = 0;
			float medianValue = 0;
			double valueRange = computeValueRange_double(data, length, &radius, &medianValue);

			float realPrecision = 0;
			if(errBoundMode==ABS)
				realPrecision = absErrBound;
			else if(errBoundMode==REL)
				realPrecision = valueRange*relBoundRatio;

			bytes = SZ_fast_compress_args_unpredictable_double(data, outSize, realPrecision, r5, r4, r3, r2, r1, medianValue, radius);		
		}		
		return bytes;
	}

}

unsigned char* SZ_fast_compress_args2_split(int fastMode, int dataType, void *data, size_t *outSize, unsigned char* outputBytes, unsigned char* chunk_arr, size_t chunk_iter,
int errBoundMode, float absErrBound, float relBoundRatio, float compressionRatio, float tolerance, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	int blockSize = 128;
	unsigned char*  bytes = NULL;	
	size_t length = computeDataLength(r5, r4, r3, r2, r1);
	
	if(dataType == SZ_FLOAT)
	{
		// if(errBoundMode==FXR)
		// {
		// 	if(fastMode == SZx_OPENMP_FAST_CMPR)
		// 	{
		// 		SZ_fast_compress_args_unpredictable_blocked_fixed_rate_float2_openmp(data, outSize, outputBytes, compressionRatio, tolerance, length, blockSize);
		// 	}
		// 	else
		// 	{
		// 		SZ_fast_compress_args_unpredictable_blocked_fixed_rate_float2(data, outSize, outputBytes, compressionRatio, tolerance, length, blockSize); 	
		// 	}
		// }
		// else //not fix rate mode
		// {
			if(fastMode == SZx_WITH_BLOCK_FAST_CMPR || fastMode == SZx_RANDOMACCESS_FAST_CMPR || fastMode == SZx_OPENMP_FAST_CMPR)
			{
				float realPrecision = absErrBound;
				if(errBoundMode==REL)
				{
					float valueRange = computeValueRange_float(data, length, NULL, NULL);
					realPrecision = valueRange*relBoundRatio;
				}

				// if (fastMode == SZx_RANDOMACCESS_FAST_CMPR) {
				// 	SZ_fast_compress_args_unpredictable_blocked_float2(data, outSize, outputBytes, realPrecision, length, blockSize);
				// } 
				if(fastMode == SZx_OPENMP_FAST_CMPR)
				{
					#ifdef _OPENMP
					SZ_fast_compress_args_unpredictable_blocked_randomaccess_float2_openmp_split(data, outSize, outputBytes, 
    												chunk_arr,chunk_iter, realPrecision,length,
                                                  blockSize);
					#else
					SZ_fast_compress_args_unpredictable_blocked_float2_split(data, outSize, outputBytes, 
    												chunk_arr,chunk_iter, realPrecision,length,
                                                  blockSize);
					printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the compression is degraded to serial version automatically.\n");
					#endif
				}
				else {
					SZ_fast_compress_args_unpredictable_blocked_float2_split(data, outSize, outputBytes, 
    												chunk_arr,chunk_iter, realPrecision,length,
                                                  blockSize);
					// SZ_fast_compress_args_unpredictable_blocked_float2(data, outSize, outputBytes, realPrecision, length, blockSize);
				}
			}
			else
			{
				//compute value range
				
				float radius = 0;
				float medianValue = 0;
				float valueRange = computeValueRange_float(data, length, &radius, &medianValue);

				float realPrecision = 0;
				if(errBoundMode==ABS)
					realPrecision = absErrBound;
				else if(errBoundMode==REL)
					realPrecision = valueRange*relBoundRatio;

				// SZ_fast_compress_args_unpredictable_float2(data, outSize, outputBytes, realPrecision, r5, r4, r3, r2, r1, medianValue, radius);		
			}			
		// }
		

	}
	// else if(dataType == SZ_DOUBLE)
	// {
	// 	if(fastMode == SZx_WITH_BLOCK_FAST_CMPR || fastMode == SZx_RANDOMACCESS_FAST_CMPR || fastMode == SZx_OPENMP_FAST_CMPR)
	// 	{
	// 		float realPrecision = absErrBound;
	// 		if(errBoundMode==REL)
	// 		{
	// 			double valueRange = computeValueRange_double(data, length, NULL, NULL);
	// 			realPrecision = valueRange*relBoundRatio;
	// 		}

	// 		int blockSize = 128;
	// 		if (fastMode == SZx_RANDOMACCESS_FAST_CMPR) {
	// 			bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_double(data, outSize, realPrecision, length, blockSize);
	// 		} 
	// 		else if(fastMode == SZx_OPENMP_FAST_CMPR)
	// 		{
	// 			#ifdef _OPENMP
	// 			bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_double_openmp(data, outSize, realPrecision, length,
	// 																						  blockSize);
	// 			#else
	// 			bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_double(data, outSize, realPrecision, length, blockSize);
	// 			printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the compression is degraded to serial version automatically.\n");
	// 			#endif
	// 		}
	// 		else {
	// 			bytes = SZ_fast_compress_args_unpredictable_blocked_double(data, outSize, realPrecision, length, blockSize);
	// 		}
	// 		return bytes;
	// 	}
	// 	else
	// 	{
	// 		//compute value range
	// 		float radius = 0;
	// 		float medianValue = 0;
	// 		double valueRange = computeValueRange_double(data, length, &radius, &medianValue);

	// 		float realPrecision = 0;
	// 		if(errBoundMode==ABS)
	// 			realPrecision = absErrBound;
	// 		else if(errBoundMode==REL)
	// 			realPrecision = valueRange*relBoundRatio;

	// 		bytes = SZ_fast_compress_args_unpredictable_double(data, outSize, realPrecision, r5, r4, r3, r2, r1, medianValue, radius);		
	// 	}		
	// 	return bytes;
	// }

}

/**
 * @deprecated
 * */
void* SZ_fast_decompress_pred(int dataType, float* preData, unsigned char *curBytes, size_t byteLength, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    int x = 1;
    char *y = (char*)&x;
    if(*y==1)
        sysEndianType = LITTLE_ENDIAN_SYSTEM;
    else //=0
        sysEndianType = BIG_ENDIAN_SYSTEM;

    if(dataType == SZ_FLOAT)
    {
        float* newFloatData = NULL;
        SZ_fast_decompress_args_with_prediction_float(&newFloatData, preData, r5, r4, r3, r2, r1, curBytes, byteLength);
        return newFloatData;
    }
    else if(dataType == SZ_DOUBLE)
    {
        double* newDoubleData = NULL;
        //SZ_fast_decompress_args_unpredictable_float(&newDoubleData, r5, r4, r3, r2, r1, bytes, byteLength, 0, NULL);
        return newDoubleData;
    }

    return NULL;
}

void* SZ_fast_decompress(int fastMode, int dataType, unsigned char *bytes, size_t byteLength, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t nbEle = computeDataLength(r5, r4, r3, r2, r1);
    int x = 1;
    char *y = (char*)&x;
    if(*y==1)
        sysEndianType = LITTLE_ENDIAN_SYSTEM;
    else //=0
        sysEndianType = BIG_ENDIAN_SYSTEM;

    if(dataType == SZ_FLOAT)
    {
        float* newFloatData = NULL;
        if(fastMode == SZx_NO_BLOCK_FAST_CMPR)
            SZ_fast_decompress_args_unpredictable_float(&newFloatData, r5, r4, r3, r2, r1, bytes, byteLength);
		else if(fastMode == SZx_WITH_BLOCK_FAST_CMPR)
			SZ_fast_decompress_args_unpredictable_blocked_float(&newFloatData, nbEle, bytes);            
        else if(fastMode == SZx_RANDOMACCESS_FAST_CMPR)
			SZ_fast_decompress_args_unpredictable_blocked_randomaccess_float(&newFloatData, nbEle, bytes);
        else //SZx_openmp
        {
#ifdef _OPENMP
                SZ_fast_decompress_args_unpredictable_blocked_randomaccess_float_openmp(&newFloatData, nbEle, bytes);
#else
                SZ_fast_decompress_args_unpredictable_blocked_float(&newFloatData, nbEle, bytes);
                printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the decompression is degraded to serial version automatically.\n");
#endif
        }
        return newFloatData;
    }
    else if(dataType == SZ_DOUBLE)
    {
        double* newFloatData = NULL;
        if(fastMode == SZx_NO_BLOCK_FAST_CMPR)
            SZ_fast_decompress_args_unpredictable_double(&newFloatData, r5, r4, r3, r2, r1, bytes, byteLength);
		else if(fastMode == SZx_WITH_BLOCK_FAST_CMPR)
			SZ_fast_decompress_args_unpredictable_blocked_double(&newFloatData, nbEle, bytes);            
        else if(fastMode == SZx_RANDOMACCESS_FAST_CMPR)
			SZ_fast_decompress_args_unpredictable_blocked_randomaccess_double(&newFloatData, nbEle, bytes);
        else //SZx_openmp
        {
#ifdef _OPENMP
                SZ_fast_decompress_args_unpredictable_blocked_randomaccess_double_openmp(&newFloatData, nbEle, bytes);
#else
                SZ_fast_decompress_args_unpredictable_blocked_double(&newFloatData, nbEle, bytes);
                printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the decompression is degraded to serial version automatically.\n");
#endif
        }
        return newFloatData;
    }

    return NULL;
}

void* SZ_fast_decompress_test(int fastMode, int dataType, unsigned char *bytes, size_t byteLength, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t nbEle = computeDataLength(r5, r4, r3, r2, r1);
    int x = 1;
    char *y = (char*)&x;
    if(*y==1)
        sysEndianType = LITTLE_ENDIAN_SYSTEM;
    else //=0
        sysEndianType = BIG_ENDIAN_SYSTEM;

    if(dataType == SZ_FLOAT)
    {
        float* newFloatData = NULL;
        if(fastMode == SZx_NO_BLOCK_FAST_CMPR)
            SZ_fast_decompress_args_unpredictable_float(&newFloatData, r5, r4, r3, r2, r1, bytes, byteLength);
		else if(fastMode == SZx_WITH_BLOCK_FAST_CMPR)
			SZ_fast_decompress_args_unpredictable_blocked_float_test(&newFloatData, nbEle, bytes);            
        else if(fastMode == SZx_RANDOMACCESS_FAST_CMPR)
			SZ_fast_decompress_args_unpredictable_blocked_randomaccess_float(&newFloatData, nbEle, bytes);
        else //SZx_openmp
        {
#ifdef _OPENMP
                SZ_fast_decompress_args_unpredictable_blocked_randomaccess_float_openmp_test(&newFloatData, nbEle, bytes);
#else
                SZ_fast_decompress_args_unpredictable_blocked_float_test(&newFloatData, nbEle, bytes);
                printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the decompression is degraded to serial version automatically.\n");
#endif
        }
        return newFloatData;
    }
    else if(dataType == SZ_DOUBLE)
    {
        double* newFloatData = NULL;
        if(fastMode == SZx_NO_BLOCK_FAST_CMPR)
            SZ_fast_decompress_args_unpredictable_double(&newFloatData, r5, r4, r3, r2, r1, bytes, byteLength);
		else if(fastMode == SZx_WITH_BLOCK_FAST_CMPR)
			SZ_fast_decompress_args_unpredictable_blocked_double(&newFloatData, nbEle, bytes);            
        else if(fastMode == SZx_RANDOMACCESS_FAST_CMPR)
			SZ_fast_decompress_args_unpredictable_blocked_randomaccess_double(&newFloatData, nbEle, bytes);
        else //SZx_openmp
        {
#ifdef _OPENMP
                SZ_fast_decompress_args_unpredictable_blocked_randomaccess_double_openmp(&newFloatData, nbEle, bytes);
#else
                SZ_fast_decompress_args_unpredictable_blocked_double(&newFloatData, nbEle, bytes);
                printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the decompression is degraded to serial version automatically.\n");
#endif
        }
        return newFloatData;
    }

    return NULL;
}

void* SZ_fast_decompress_split(int fastMode, int dataType, float *newData, unsigned char *bytes, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t nbEle = computeDataLength(r5, r4, r3, r2, r1);
    int x = 1;
    char *y = (char*)&x;
    if(*y==1)
        sysEndianType = LITTLE_ENDIAN_SYSTEM;
    else //=0
        sysEndianType = BIG_ENDIAN_SYSTEM;

    if(dataType == SZ_FLOAT)
    {
        // float* newFloatData = NULL;
        // if(fastMode == SZx_NO_BLOCK_FAST_CMPR)
            // SZ_fast_decompress_args_unpredictable_float(&newFloatData, r5, r4, r3, r2, r1, bytes, byteLength);
		if(fastMode == SZx_WITH_BLOCK_FAST_CMPR) //{previous} else if(fastMode == SZx_WITH_BLOCK_FAST_CMPR)
		{
			// printf("Decompress start\n");
			SZ_fast_decompress_args_unpredictable_blocked_float2_split(newData, nbEle, bytes);
			// printf("Decompress end\n");
		}
			            
        // else if(fastMode == SZx_RANDOMACCESS_FAST_CMPR)
			// SZ_fast_decompress_args_unpredictable_blocked_randomaccess_float(&newFloatData, nbEle, bytes);
        else //SZx_openmp
        {
#ifdef _OPENMP
                // SZ_fast_decompress_args_unpredictable_blocked_randomaccess_float_openmp(&newFloatData, nbEle, bytes);
				SZ_fast_decompress_args_unpredictable_blocked_randomaccess_float2_openmp(newData, nbEle, bytes);
#else
                SZ_fast_decompress_args_unpredictable_blocked_float2_split(newData, nbEle, bytes); // The split does not have any special optimization compared with float2
                printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the decompression is degraded to serial version automatically.\n");
#endif
        }
    }
//     else if(dataType == SZ_DOUBLE)
//     {
//         double* newFloatData = NULL;
//         if(fastMode == SZx_NO_BLOCK_FAST_CMPR)
//             SZ_fast_decompress_args_unpredictable_double(&newFloatData, r5, r4, r3, r2, r1, bytes, byteLength);
// 		else if(fastMode == SZx_WITH_BLOCK_FAST_CMPR)
// 			SZ_fast_decompress_args_unpredictable_blocked_double(&newFloatData, nbEle, bytes);            
//         else if(fastMode == SZx_RANDOMACCESS_FAST_CMPR)
// 			SZ_fast_decompress_args_unpredictable_blocked_randomaccess_double(&newFloatData, nbEle, bytes);
//         else //SZx_openmp
//         {
// #ifdef _OPENMP
//                 SZ_fast_decompress_args_unpredictable_blocked_randomaccess_double_openmp(&newFloatData, nbEle, bytes);
// #else
//                 SZ_fast_decompress_args_unpredictable_blocked_double(&newFloatData, nbEle, bytes);
//                 printf("WARNING: It seems that you want to run the code with openmp mode but you didn't compile the code in openmp mode.\nSo, the decompression is degraded to serial version automatically.\n");
// #endif
//         }
//         return newFloatData;
//     }

    return NULL;
}
