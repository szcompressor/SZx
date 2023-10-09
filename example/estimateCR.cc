/**
 *  @file estimateCR.c
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
#include "szx_float.h"

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;


void szx_cost_start()
{
	totalCost = 0;
	gettimeofday(&costStart, NULL);
}

void szx_cost_end()
{
	double elapsed;
	struct timeval costEnd;
	gettimeofday(&costEnd, NULL);
	elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
	totalCost += elapsed;
}


int main(int argc, char * argv[])
{
    char oriFilePath[640];
    
    if(argc < 2)
    {
		printf("Usage: estimateCR [srcFilePath] [err bound]\n");
		printf("Example: estimateCR testfloat_8_8_128.dat 1E-3\n");
		exit(0);
    }
   
    sprintf(oriFilePath, "%s", argv[1]);
    float errorBound = atof(argv[2]); 
  
    int status = 0; 
    size_t nbEle;
    float *data = SZx_readFloatData(oriFilePath, &nbEle, &status);
    if(status != SZx_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
   
    szx_cost_start();
    int samplingStride = 10;
    int blockSize = 128;
    float* radiusArray = NULL;
    float* mediusArray = NULL;
    float* buffer = NULL;
    //float approximateValueRange = computeRadiusBuffer(data, SZx_FLOAT, nbEle, samplingStride, blockSize, &radiusArray, &mediusArray, &buffer);    
    size_t sumReqNbB = 0, sum_actual_leadNum = 0;
    float CR = estimateCRbasedonErrorBound_buffered(SZx_FLOAT, errorBound, buffer, mediusArray, radiusArray, samplingStride, blockSize, nbEle, &sumReqNbB, &sum_actual_leadNum);	
    szx_cost_end();
    printf("timecost=%f\n",totalCost); 
    printf("estimated compression ratio = %f\n", CR);

    printf("done\n");
    free(data);
    
    return 0;
}
