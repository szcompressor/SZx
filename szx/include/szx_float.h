/**
 *  @file sz_float.h
 *  @author Sheng Di
 *  @date July, 2017
 *  @brief Header file for the sz_float.c.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _SZ_Float_H
#define _SZ_Float_H

#ifdef __cplusplus
extern "C" {
#endif

unsigned char * SZ_fast_compress_args_with_prediction_float(float *pred, float *data, size_t *outSize, float absErrBound, size_t r5,
                                            size_t r4, size_t r3, size_t r2, size_t r1, float medianValue, float radius);

void SZ_fast_compress_args_unpredictable_one_block_float(float *oriData, size_t nbEle, float absErrBound,
                                                                unsigned char *outputBytes, int *outSize,
                                                                unsigned char *leadNumberArray_int, float medianValue,
                                                                float radius);
                  
float computeRadiusBuffer_float(float *oriData, size_t nbEle, int samplingRate, int blockSize, float** radiusArray, float** mediusArray, float** buffer);
size_t computeStateMedianRadius_float(float *oriData, size_t nbEle, float absErrBound, int blockSize,
                                      unsigned char *stateArray, float *medianArray, float *radiusArray) ;
                                      
void max_min_float(float *x, int n, float *tmp_max, float *tmp_min);

void simd_max_min_float(float *x, int n, float *tmp_max, float *tmp_min);

void computeStateMedianRadius_float2(float *oriData, size_t nbEle, float absErrBound,
                                     unsigned char *state, float *median, float *radius) ;

float computeValueRange_float(float* oriData, size_t length, float* radius, float* medianValue);

float estimateCRbasedonErrorBound_buffered_float(float errorBound, float* buffer, float* medianArray, float* radiusArray, int samplingRate, int blockSize, size_t nbEle, 
size_t *sumReqNbBytes, size_t *sum_actual_leadNumbers);

float estimateCRbasedonErrorBound_float(float errorBound, float* data, int blockSize, size_t nbEle);                                     

float estimateErrorBoundbasedonCR_buffered_float(float targetCompressionRatio, float tolerance, int samplingRate, float initErrorBound, int blockSize, size_t nbEle, 
float* buffer, float* medianArray, float* radiusArray);

float estimateErrorBoundbasedonCR_float(float targetCompressionRatio, float tolerance, float* data, float initErrorBound, int blockSize, size_t nbEle);
unsigned char *
SZ_fast_compress_args_unpredictable_blocked_fixed_rate_float(float *oriData, size_t *outSize, float compressionRatio, float tolerance, size_t nbEle,
                                                  int blockSize) ;

void SZ_fast_compress_args_unpredictable_blocked_fixed_rate_float2(float *oriData, size_t *outSize, unsigned char* outputBytes, float compressionRatio, float tolerance, size_t nbEle,
                                                  int blockSize) ;

void SZ_fast_compress_args_unpredictable_blocked_fixed_rate_float2_openmp(float *oriData, size_t *outSize, unsigned char* outputBytes, float compressionRatio, float tolerance, size_t nbEle,
                                                  int blockSize) ;

size_t SZ_fast_compress_args_unpredictable_blocked_args_float(float *oriData, unsigned char* outputBytes, float absErrBound, size_t nbEle,
                                                  int blockSize);

unsigned char *
SZ_fast_compress_args_unpredictable_blocked_float(float *oriData, size_t *outSize, float absErrBound, size_t nbEle,
                                                  int blockSize) ;

unsigned char *
SZ_fast_compress_args_unpredictable_blocked_float_test(float *oriData, size_t *outSize, float absErrBound, size_t nbEle,
                                                  int blockSize) ;                                            
                                                  
void SZ_fast_compress_args_unpredictable_blocked_float2(float *oriData, size_t *outSize, unsigned char* outputBytes, float absErrBound, size_t nbEle,
                                                  int blockSize) ; 

void SZ_fast_compress_args_unpredictable_blocked_float2_split(float *oriData, size_t *outSize, unsigned char* outputBytes, 
    unsigned char* chunk_arr, size_t chunk_iter, float absErrBound, size_t nbEle,
                                                  int blockSize);                                                                                        
                                                  
unsigned char *
SZ_fast_compress_args_unpredictable_blocked_randomaccess_float_openmp(float *oriData, size_t *outSize, float absErrBound,
                                                               size_t nbEle, int blockSize) ;

unsigned char *
SZ_fast_compress_args_unpredictable_blocked_randomaccess_float_openmp_test(float *oriData, size_t *outSize, float absErrBound,
                                                               size_t nbEle, int blockSize) ;

unsigned char *
SZ_fast_compress_args_unpredictable_blocked_randomaccess_float2_openmp(float *oriData, size_t *outSize, unsigned char* outputBytes, float absErrBound,
                                                               size_t nbEle, int blockSize);

unsigned char *
SZ_fast_compress_args_unpredictable_blocked_randomaccess_float(float *oriData, size_t *outSize, float absErrBound,
    size_t nbEle, int blockSize) ;
    
unsigned char *
SZ_fast_compress_args_unpredictable_float(float *data, size_t *outSize, float absErrBound, size_t r5, size_t r4,
                                          size_t r3, size_t r2, size_t r1, float mValue, float radius);
                                          
unsigned char *SZ_skip_compress_float(float *data, size_t dataLength, size_t *outSize) ;

void computeReqLength_float(double realPrecision, short radExpo, int *reqLength, float *medianValue) ;



#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _SZ_Float_H  ----- */

