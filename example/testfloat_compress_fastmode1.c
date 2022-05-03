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
#include "szx.h"
#include "szx_rw.h"

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
    
    if(argc < 2)
    {
		printf("Usage: testfloat_compress_fastmode1 [srcFilePath] [err bound]\n");
		printf("Example: testfloat_compress_fastmode1 testfloat_8_8_128.dat 1E-3\n");
		exit(0);
    }
   
    sprintf(oriFilePath, "%s", argv[1]);
    float errorBound = atof(argv[2]); 
    sprintf(outputFilePath, "%s.szx", oriFilePath);
  
    int status = 0; 
    size_t nbEle;
    float *data = readFloatData(oriFilePath, &nbEle, &status);
    if(status != SZ_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
   
    size_t outSize; 
    cost_start();
    unsigned char* bytes =  SZ_fast_compress_args(SZx_NO_BLOCK_FAST_CMPR, SZ_FLOAT, data, &outSize, ABS, errorBound, 0.001, 0, 0, 0, 0, nbEle);
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
    
    return 0;
}
