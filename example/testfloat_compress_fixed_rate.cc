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
 
    if(argc < 3)
    {
		printf("Usage: testfloat_compress_fixed_rate [srcFilePath] [target ratio] [tolerance]\n");
		printf("Example: testfloat_compress_fixed_rate testfloat_8_8_128.dat 10 0.1\n");
		exit(0);
    }
   
    sprintf(oriFilePath, "%s", argv[1]);
    float targetCompressionRatio = atof(argv[2]); 
    float tolerance = atof(argv[3]);
    sprintf(outputFilePath, "%s.szx", oriFilePath);
  
    int status = 0; 
    size_t nbEle;
    float *data = SZx_readFloatData(oriFilePath, &nbEle, &status);
    if(status != SZx_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
   
    int blockSize = 128;
    float radius = 0, medianValue = 0;
    float valueRange = SZx_computeValueRange(data, SZx_FLOAT, nbEle, &radius, &medianValue);
    float initErrorBound = valueRange*1E-3;
    cost_start();

    float errorBound = SZx_estimateErrorBoundbasedonCR(SZx_FLOAT, targetCompressionRatio, tolerance, data, initErrorBound, blockSize, nbEle);
    cost_end();
    printf("searched error bound = %.30G\n", errorBound);
    size_t outSize = 0;
    unsigned char* bytes = SZx_fast_compress_args_unpredictable_blocked(data, SZx_FLOAT, &outSize, errorBound, nbEle, blockSize);

    printf("timecost=%f, %d\n",totalCost, bytes[0]); 
    printf("compression size = %zu, CR = %f\n", outSize, 1.0f*nbEle*sizeof(float)/outSize);
    SZx_writeByteData(bytes, outSize, outputFilePath, &status);
    if(status != SZx_SCES)
    {
        printf("Error: data file %s cannot be written!\n", outputFilePath);
        exit(0);
    }

//    printf("done\n");
    free(bytes); 
    free(data);
    
    return 0;
}
