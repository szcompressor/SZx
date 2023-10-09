/**
 *  @file szx_double.h
 *  @author Sheng Di
 *  @date July, 2017
 *  @brief Header file for the sz_double.c.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <szx_float.h>

#ifndef _SZx_Double_H
#define _SZx_Double_H

namespace szx{	
	
void SZx_fast_compress_args_unpredictable_one_block_double(double *oriData, size_t nbEle, float absErrBound,
                                                                unsigned char *outputBytes, int *outSize,
                                                                unsigned char *leadNumberArray_int, float medianValue,
                                                                float radius);
                                                                
size_t computeStateMedianRadius_double(double *oriData, size_t nbEle, float absErrBound, int blockSize,
                                      unsigned char *stateArray, float *medianArray, float *radiusArray) ;
                                      
void max_min_double(double *x, int n, double *tmp_max, double *tmp_min);

void simd_max_min_double(double *x, int n, double *tmp_max, double *tmp_min);

void computeStateMedianRadius_double2(double *oriData, size_t nbEle, float absErrBound,
                                     unsigned char *state, float *median, float *radius) ;
                                     
double computeValueRange_double(double* oriData, size_t length, float* radius, float* medianValue);
                                     
// added

float computeRadiusBuffer_double(double *oriData, size_t nbEle, int samplingRate, int blockSize, float** radiusArray, float** mediusArray, double** buffer);

float estimateCRbasedonErrorBound_buffered_double(float errorBound, double* buffer, float* medianArray, float* radiusArray, int samplingRate, int blockSize, size_t nbEle, 
size_t *sumReqNbBytes, size_t *sum_actual_leadNumbers);

float estimateCRbasedonErrorBound_double(float errorBound, double* data, int blockSize, size_t nbEle);                                     

float estimateErrorBoundbasedonCR_buffered_double(float targetCompressionRatio, float tolerance, int samplingRate, float initErrorBound, int blockSize, size_t nbEle, 
double* buffer, float* medianArray, float* radiusArray);

float estimateErrorBoundbasedonCR_double(float targetCompressionRatio, float tolerance, double* data, float initErrorBound, int blockSize, size_t nbEle);
unsigned char *
SZx_fast_compress_args_unpredictable_blocked_fixed_rate_double(double *oriData, size_t *outSize, float compressionRatio, float tolerance, size_t nbEle,
                                                  int blockSize) ;

void SZx_fast_compress_args_unpredictable_blocked_fixed_rate_double2(double *oriData, size_t *outSize, unsigned char* outputBytes, float compressionRatio, float tolerance, size_t nbEle,
                                                  int blockSize) ;

void SZx_fast_compress_args_unpredictable_blocked_double2(double *oriData, size_t *outSize, unsigned char* outputBytes, float absErrBound, size_t nbEle,
                                                  int blockSize) ;   

// added

unsigned char *
SZx_fast_compress_args_unpredictable_blocked_double(double *oriData, size_t *outSize, float absErrBound, size_t nbEle,
                                                  int blockSize) ;
                                                  
unsigned char *
SZx_fast_compress_args_unpredictable_blocked_randomaccess_double_openmp(double *oriData, size_t *outSize, 
									float absErrBound, size_t nbEle, int blockSize) ;
                                                               
                                                               
unsigned char *
SZx_fast_compress_args_unpredictable_blocked_randomaccess_double(double *oriData, size_t *outSize, 
								float absErrBound, size_t nbEle, int blockSize) ;
    
unsigned char *
SZx_fast_compress_args_unpredictable_double(double *data, size_t *outSize, float absErrBound, size_t r5, size_t r4,
                                          size_t r3, size_t r2, size_t r1, float mValue, float radius);
                                          
unsigned char *SZx_skip_compress_double(double *data, size_t dataLength, size_t *outSize) ;

void computeReqLength_double(float realPrecision, short radExpo, int *reqLength, float *medianValue) ;

}

#endif /* ----- #ifndef _SZx_Double_H  ----- */

