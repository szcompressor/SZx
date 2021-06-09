#ifndef CUSZX_ENTRY_H
#define CUSZX_ENTRY_H

#include <stdio.h>
#include "cuSZx_compress_float.h"
#include "cuSZx_decompress_float.h"

unsigned char* cuSZx_fast_compress_args_unpredictable_blocked_float(float *oriData, size_t *outSize, float absErrBound, size_t nbEle, int blockSize, unsigned char *test);
void cuSZx_fast_decompress_args_unpredictable_blocked_float(float** newData, size_t nbEle, unsigned char* cmpBytes, float* m);

#endif /* ----- #ifndef CUSZX_ENTRY_H  ----- */
