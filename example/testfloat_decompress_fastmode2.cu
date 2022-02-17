/**
 *  @file test_decompress.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using Decompression interface.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
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
    size_t nbEle, totalNbEle;
    char zipFilePath[640], outputFilePath[645];
    if(argc < 2)
    {
		printf("Usage: testfloat_decompress_fastmode2 [srcFilePath] [nbEle] [--cuda]\n");
		printf("Example: testfloat_decompress_fastmode2 testfloat_8_8_128.dat.sz 8192 --cuda\n");
		exit(0);
	}

    sprintf(zipFilePath, "%s", argv[1]);
//    nbEle = atoi(argv[2]);
    nbEle = strtoimax(argv[2], NULL, 10);
    bool withGPU = false;
    if (argc > 3 && !strcmp(argv[3], "--cuda")) withGPU = true;

    sprintf(outputFilePath, "%s.out", zipFilePath);

    size_t byteLength;
    int status;
    unsigned char *bytes = readByteData(zipFilePath, &byteLength, &status);
    if(status!=SZ_SCES)
    {
        printf("Error: %s cannot be read!\n", zipFilePath);
        exit(0);
    }


    float *data = NULL;
    if (withGPU)
    {
        cuSZx_fast_decompress_args_unpredictable_blocked_float(&data, nbEle, bytes);
    }else{
        cost_start();
        SZ_fast_decompress_args_unpredictable_blocked_float(&data, nbEle, bytes);
//    SZ_fast_decompress_args_unpredictable_blocked_randomaccess_float(&data, nbEle, bytes);
        cost_end();
        printf("timecost=%f\n",totalCost);
    }

    free(bytes);
    writeFloatData_inBytes(data, nbEle, outputFilePath, &status);
    if(status!=SZ_SCES)
    {
	printf("Error: %s cannot be written!\n", outputFilePath);
	exit(0);
    }
    printf("done\n");

    char oriFilePath[645];
    strcpy(oriFilePath, zipFilePath);
    oriFilePath[strlen(zipFilePath)-4] = '\0';
    float *ori_data = readFloatData(oriFilePath, &totalNbEle, &status);
    if(status!=SZ_SCES)
    {
        printf("Error: %s cannot be read!\n", oriFilePath);
        exit(0);
    }

    size_t i = 0;
    float Max = 0, Min = 0, diffMax = 0;
    Max = ori_data[0];
    Min = ori_data[0];
    diffMax = fabs(data[0] - ori_data[0]);
    double sum1 = 0, sum2 = 0;
    for (i = 0; i < nbEle; i++)
    {
        sum1 += ori_data[i];
		sum2 += data[i];
    }
    double mean1 = sum1/nbEle;
    double mean2 = sum2/nbEle;

    double sum3 = 0, sum4 = 0;
    double sum = 0, prodSum = 0, relerr = 0;

    double maxpw_relerr = 0;
    for (i = 0; i < nbEle; i++)
    {
        if (Max < ori_data[i]) Max = ori_data[i];
        if (Min > ori_data[i]) Min = ori_data[i];

        float err = fabs(data[i] - ori_data[i]);
	if(ori_data[i]!=0)
	{
		if(fabs(ori_data[i])>1)
			relerr = err/ori_data[i];
		else
			relerr = err;
		if(maxpw_relerr<relerr)
			maxpw_relerr = relerr;
        }

	/*if(err > 1600000)
	{
		printf("i=%zu, ori=%f, dec=%f, diff=%f\n", i, ori_data[i], data[i], err);
		exit(0);
	}*/
	if (diffMax < err)
		diffMax = err;
        prodSum += (ori_data[i]-mean1)*(data[i]-mean2);
        sum3 += (ori_data[i] - mean1)*(ori_data[i]-mean1);
        sum4 += (data[i] - mean2)*(data[i]-mean2);
	sum += err*err;
    }
    double std1 = sqrt(sum3/nbEle);
    double std2 = sqrt(sum4/nbEle);
    double ee = prodSum/nbEle;
    double acEff = ee/std1/std2;

    double mse = sum/nbEle;
    double range = Max - Min;
    double psnr = 20*log10(range)-10*log10(mse);
    double nrmse = sqrt(mse)/range;

    double compressionRatio = 1.0*nbEle*sizeof(float)/byteLength;

    printf ("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
    printf ("Max absolute error = %.10f\n", diffMax);
    printf ("Max relative error = %f\n", diffMax/(Max-Min));
    printf ("Max pw relative error = %f\n", maxpw_relerr);
    printf ("PSNR = %f, NRMSE= %.20G\n", psnr,nrmse);
    printf ("acEff=%f\n", acEff);
    printf ("compressionRatio = %f\n", compressionRatio);
    free(data);
    free(ori_data);
    return 0;
}
