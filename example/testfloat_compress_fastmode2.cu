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
#include "sz.h"
#include "rw.h"
#include "cuSZx_entry.h"

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
    char *cfgFile;
    
    if(argc < 4)
    {
		printf("Usage: testfloat_compress_fastmode2 [config_file] [srcFilePath] [block size] [err bound]\n");
		printf("Example: testfloat_compress_fastmode2 sz.config testfloat_8_8_128.dat 64 1E-3\n");
		exit(0);
    }
   
    cfgFile=argv[1];
    sprintf(oriFilePath, "%s", argv[2]);
    int blockSize = atoi(argv[3]);
    float errBound = atof(argv[4]);
    printf("cfgFile=%s\n", cfgFile); 
    int status = SZ_Init(cfgFile);
    if(status == SZ_NSCS)
	exit(0);
    sprintf(outputFilePath, "%s.sz", oriFilePath);
   
    size_t nbEle;
    float *data = readFloatData(oriFilePath, &nbEle, &status);
    if(status != SZ_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
    //float *revValue = (float *)malloc(sizeof(float));
    //*revValue = 1.0E36;
    unsigned char *test_meta = (unsigned char*)malloc(sizeof(float) * nbEle * sizeof(unsigned char));
    memset(test_meta, 0, sizeof(float) * nbEle * sizeof(unsigned char));
    //int *test_meta = (int*)malloc(nbEle/blockSize*sizeof(int));
    //memset(test_meta, 0, nbEle/blockSize*sizeof(int));
   
    size_t outSize; 
    cost_start();
    unsigned char* bytes = SZ_fast_compress_args_unpredictable_blocked_float(data, &outSize, errBound, nbEle, blockSize, test_meta);
    cuSZx_fast_compress_args_unpredictable_blocked_float(data, &outSize, errBound, nbEle, blockSize, test_meta);
    //unsigned char* bytes =  SZ_fast_compress_args(SZ_WITH_BLOCK_FAST_CMPR, SZ_FLOAT, data, &outSize, ABS, errBound, 0.001, 0, 0, 0, 0, 0, nbEle);
    cost_end();
    printf("timecost=%f, %d\n",totalCost, bytes[0]); 
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
    SZ_Finalize();
    
    return 0;}
