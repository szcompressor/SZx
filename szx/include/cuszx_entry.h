#ifndef CUSZX_ENTRY_H
#define CUSZX_ENTRY_H

#include <stdio.h>
#include "cuszx_float.h"
#include "cuszxd_float.h"

#define GPU

#ifdef __cplusplus
extern "C" {
#endif

unsigned char* cuSZx_fast_compress_args_unpredictable_blocked_float(float *oriData, size_t *outSize, float absErrBound, size_t nbEle, int blockSize);
void cuSZx_fast_decompress_args_unpredictable_blocked_float(float** newData, size_t nbEle, unsigned char* cmpBytes);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef CUSZX_ENTRY_H  ----- */
