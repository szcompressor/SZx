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
#include "omp.h"
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

void computeStateMedianRadius_float2(float *oriData, size_t nbEle, float absErrBound, int blockSize,
                                     unsigned char *stateArray, float *medianArray, float *radiusArray,
                                     size_t *nbConstantBlocks, size_t *constantBlocks,
                                     size_t *nbNonConstantBlocks, size_t *nonConstantBlocks) {
//    size_t nbConstantBlocks = 0;
    size_t i = 0, j = 0;
    size_t nbBlocks = nbEle / blockSize;
    size_t offset = 0;

    for (i = 0; i < nbBlocks; i++) {
        float min = oriData[offset];
        float max = oriData[offset];
        for (j = 1; j < blockSize; j++) {
            float v = oriData[offset + j];
            if (min > v)
                min = v;
            else if (max < v)
                max = v;
        }
        float valueRange = max - min;
        float radius = valueRange / 2;
        float medianValue = min + radius;

        if (radius <= absErrBound) {
            stateArray[i] = 0;
            constantBlocks[(*nbConstantBlocks)++] = i;
        } else {
            stateArray[i] = 1;
            nonConstantBlocks[(*nbNonConstantBlocks)++] = i;
        }
        stateArray[i] = radius <= absErrBound ? 0 : 1;
        medianArray[i] = medianValue;
        radiusArray[i] = radius;
        offset += blockSize;
    }

    int remainCount = nbEle % blockSize;
    if (remainCount != 0) {
        float min = oriData[offset];
        float max = oriData[offset];
        for (j = 1; j < remainCount; j++) {
            float v = oriData[offset + j];
            if (min > v)
                min = v;
            else if (max < v)
                max = v;
        }
        float valueRange = max - min;
        float radius = valueRange / 2;
        float medianValue = min + radius;
        if (radius <= absErrBound) {
            stateArray[i] = 0;
            constantBlocks[(*nbConstantBlocks)++] = i;
        } else {
            stateArray[i] = 1;
            nonConstantBlocks[(*nbNonConstantBlocks)++] = i;
        }

        medianArray[i] = medianValue;
        radiusArray[i] = radius;
    }
//    return nbConstantBlocks;
}

unsigned char *
SZ_fast_compress_args_unpredictable_blocked_randomaccess_float_2(float *oriData, size_t *outSize, float absErrBound,
                                                               size_t nbEle, int blockSize) {
//#ifdef _OPENMP
    printf("use openmp\n");
    float *op = oriData;

    *outSize = 0;
    size_t maxPreservedBufferSize =
            sizeof(float) * nbEle; //assume that the compressed data size would not exceed the original size
    unsigned char *outputBytes = (unsigned char *) malloc(maxPreservedBufferSize);
    memset(outputBytes, 0, maxPreservedBufferSize);

    size_t i = 0;

    size_t nbBlocks = nbEle / blockSize;
    size_t remainCount = nbEle % blockSize;
    size_t stateNBBytes =
            remainCount == 0 ? (nbBlocks % 8 == 0 ? nbBlocks / 8 : nbBlocks / 8 + 1) : ((nbBlocks + 1) % 8 == 0 ?
                                                                                        (nbBlocks + 1) / 8 :
                                                                                        (nbBlocks + 1) / 8 + 1);
    size_t actualNBBlocks = remainCount == 0 ? nbBlocks : nbBlocks + 1;

    unsigned char *stateArray = (unsigned char *) malloc(nbBlocks);
    float *medianArray = (float *) malloc(actualNBBlocks * sizeof(float));
    float *radiusArray = (float *) malloc(actualNBBlocks * sizeof(float));

    size_t nbConstantBlocks = 0, nbNonConstantBlocks = 0;
    size_t *constantBlocks = (size_t *) malloc(actualNBBlocks * sizeof(size_t));
    size_t *nonConstantBlocks = (size_t *) malloc(actualNBBlocks * sizeof(size_t));
    computeStateMedianRadius_float2(oriData, nbEle, absErrBound, blockSize, stateArray, medianArray, radiusArray,
                                    &nbConstantBlocks, constantBlocks, &nbNonConstantBlocks, nonConstantBlocks);

    assert(nbConstantBlocks + nbNonConstantBlocks == actualNBBlocks);
    unsigned char *r = outputBytes; // + sizeof(size_t) + stateNBBytes;
    r[0] = SZ_VER_MAJOR;
    r[1] = SZ_VER_MINOR;
    r[2] = SZ_VER_SUPERFAST;
    r[3] = 1; //support random access decompression
    r[4] = (unsigned char) blockSize;
    r = r + 5; //1 byte
    sizeToBytes(r, nbConstantBlocks);
    r += sizeof(size_t); //r is the starting address of 'block-size array'

    unsigned char *R = r + nbNonConstantBlocks; //R is the starting address of the state array
    unsigned char *p = R + stateNBBytes; //p is the starting address of constant median values.
    unsigned char *q =
            p + sizeof(float) * nbConstantBlocks; //q is the starting address of the non-constant data sblocks
    unsigned char *q0 = q;
    //3: versions, 1: metadata: state, 1: metadata: blockSize, sizeof(size_t): nbConstantBlocks, ....
    *outSize = (3 + 1 + 1 + sizeof(size_t) + nbNonConstantBlocks + stateNBBytes + sizeof(float) * nbConstantBlocks);

    size_t nonConstantBlockID = 0;
    //printf("nbConstantBlocks = %zu, percent = %f\n", nbConstantBlocks, 1.0f*(nbConstantBlocks*blockSize)/nbEle);
    unsigned char *leadNumberArray_int = (unsigned char *) malloc(blockSize * sizeof(int) * nbNonConstantBlocks);
    unsigned char *tmp_q = (unsigned char *) malloc(blockSize * sizeof(float) * nbNonConstantBlocks);

#pragma omp parallel for
    for (i = 0; i < nbNonConstantBlocks; i++) {
//        printf(" Thread %d: %d\n", omp_get_thread_num(), i);
        int oSize = 0;
        SZ_fast_compress_args_unpredictable_one_block_float(op + nonConstantBlocks[i] * blockSize, blockSize,
                                                            absErrBound,
                                                            tmp_q + i * blockSize, &oSize,
                                                            leadNumberArray_int + i * blockSize,
                                                            medianArray[nonConstantBlocks[i]],
                                                            radiusArray[nonConstantBlocks[i]]);
#pragma omp critical
        {
            memcpy(q, tmp_q + i * blockSize, oSize);
            q += oSize;
            r[i] = oSize;
        }
    }
    *outSize += q - q0;

#pragma omp parallel for num_threads(4)
    for (i = 0; i < nbConstantBlocks; i++) {
        floatToBytes(p + i * sizeof(float), medianArray[constantBlocks[i]]);
    }
    p += nbConstantBlocks * sizeof(float);

    if (remainCount != 0) {
        int oSize = 0;
        if (stateArray[i]) {
            SZ_fast_compress_args_unpredictable_one_block_float(op, remainCount, absErrBound, q, &oSize,
                                                                leadNumberArray_int, medianArray[i], radiusArray[i]);
            *outSize += oSize;
            r[nbNonConstantBlocks] = oSize;
        } else {
            floatToBytes(p, medianArray[i]);
        }

    }

    convertIntArray2ByteArray_fast_1b_args(stateArray, actualNBBlocks, R);

    free(leadNumberArray_int);
    free(constantBlocks);
    free(nonConstantBlocks);
    free(tmp_q);
    return outputBytes;
}


int main(int argc, char * argv[])
{
    char oriFilePath[640], outputFilePath[645];
    char *cfgFile;
#ifdef _OPENMP
    omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel for
    for(int n = 0; n < 10; n++) {
        printf(" ThreadID = %d, Total= %d\n", omp_get_thread_num(), omp_get_num_threads());
    }
#endif
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

    size_t outSize;
    cost_start();
//    unsigned char* bytes = SZ_fast_compress_args_unpredictable_blocked_float(data, &outSize, errBound, nbEle, blockSize);
    unsigned char* bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_float_2(data, &outSize, errBound, nbEle, blockSize);

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
