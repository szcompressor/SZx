// #include <immintrin.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
typedef float data_type;
#define ALL_REDUCE_TAG 14
#define MPIR_ALLGATHERV_TAG 8
#define MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE 32768
#define fast_mode 2
#define cpr_mode ABS
// #define fast_mode 2
// #define cpr_mode FXR
#define absErrBound (1e-3)
#define relBoundRatio (1e-5)
#define tag_base 100