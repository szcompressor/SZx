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
#define MPIR_BCAST_TAG 2
#define MPIR_SCATTER_TAG 5
#define MPI_IN_PLACE (void *) -1
// #define MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE 32768
//for error analysis
#define MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE 42347060
// #define MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE 327680
#define MPIR_CVAR_ALLGATHERV_CPR_PIPELINE_MSG_SIZE 983040
#define SINGLETHREAD_CHUNK_SIZE 5120
#define MULTITHREAD_CHUNK_SIZE 256000
#define fast_mode 2
#define cpr_mode ABS
// #define fast_mode 2
// #define cpr_mode FXR
// #define absErrBound (1e-2)
#define relBoundRatio (1e-5)
#define SAVE_CPR_RESULT 0
#define PRINT_DETAILS 0
#define PRINT_EXPLANATION 0
#define PRINT_EXPERIMENTS 1
#define tag_base 100

