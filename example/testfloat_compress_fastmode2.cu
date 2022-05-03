/**
 *  @file test_compress.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "szx.h"
#include "szx_rw.h"
#include "cuszx_entry.h"

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;


void cost_start()
{
	totalCost = 0;
	gettimeofday(&costStart, NULL);
}

void cost_end()
{
	double elapsed;
	struct timeval costEnd;
	gettimeofday(&costEnd, NULL);
	elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
	totalCost += elapsed;
}


int main(int argc, char * argv[])
{
    char oriFilePath[640], outputFilePath[645];
    if(argc < 3)
    {
		printf("Usage: testfloat_compress_fastmode2 [srcFilePath] [block size] [err bound] [--cuda]\n");
		printf("Example: testfloat_compress_fastmode2 testfloat_8_8_128.dat 64 1E-3 --cuda\n");
		exit(0);
    }

    sprintf(oriFilePath, "%s", argv[1]);
    int blockSize = atoi(argv[2]);
    float errBound = atof(argv[3]);
    bool withGPU = false;
    if (argc > 4 && !strcmp(argv[4], "--cuda")) withGPU = true;

    sprintf(outputFilePath, "%s.szx", oriFilePath);

    int status = 1;
    size_t nbEle;
    float *data = readFloatData(oriFilePath, &nbEle, &status);
    if(status != SZ_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
    //float *revValue = (float *)malloc(sizeof(float));
    //*revValue = 1.0E36;

    size_t outSize;
    unsigned char* bytes = NULL;
    if (withGPU)
    {
        bytes = cuSZx_fast_compress_args_unpredictable_blocked_float(data, &outSize, errBound, nbEle, blockSize);
    }else{
        cost_start();
        bytes = SZ_fast_compress_args_unpredictable_blocked_float(data, &outSize, errBound, nbEle, blockSize);
    //    unsigned char* bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_float(data, &outSize, errBound, nbEle, blockSize);

    //    unsigned char* bytes =  SZ_fast_compress_args(SZ_WITH_BLOCK_FAST_CMPR, SZ_FLOAT, data, &outSize, ABS, errBound, errBound, 0, 0, 0, 0, 0, nbEle);
        cost_end();
        printf("\ntimecost=%f, total fastmode2\n",totalCost);
    }
    printf("compression size = %zu, CR = %f\n", outSize, 1.0f*nbEle*sizeof(float)/outSize);
    writeByteData(bytes, outSize, outputFilePath, &status);
    if(status != SZ_SCES)
    {
        printf("Error: data file %s cannot be written!\n", outputFilePath);
        exit(0);
    }

    printf("done\n");
    free(bytes);
    free(data);

    return 0;
}
