/**
 *  @file test_compress.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "sz.h"
#include "rw.h"

#ifdef _OPENMP

#include "omp.h"

#endif
struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;


void cost_start() {
    totalCost = 0;
    gettimeofday(&costStart, NULL);
}

void cost_end() {
    double elapsed;
    struct timeval costEnd;
    gettimeofday(&costEnd, NULL);
    elapsed = ((costEnd.tv_sec * 1000000 + costEnd.tv_usec) - (costStart.tv_sec * 1000000 + costStart.tv_usec)) /
              1000000.0;
    totalCost += elapsed;
}

template<typename Type>
void writefile(const char *file, Type *data, size_t num_elements) {
    std::ofstream fout(file, std::ios::binary);
    fout.write(reinterpret_cast<const char *>(&data[0]), num_elements * sizeof(Type));
    fout.close();
}

int main(int argc, char *argv[]) {
    char oriFilePath[640], outputFilePath[645];
    char *cfgFile;
    if (argc < 4) {
        printf("Usage: testfloat_compress_fastmode6 [config_file] [srcFilePath] [block size] [reb] [abs]\n");
        printf("Example: testfloat_compress_fastmode6 sz.config testfloat_8_8_128.dat 500 1E-3 1E-3\n");
        exit(0);
    }

    cfgFile = argv[1];
    sprintf(oriFilePath, "%s", argv[2]);
    int blockSize = atoi(argv[3]);
    float reb = atof(argv[4]);
    float abs = atof(argv[5]);

    if (argc >= 7) {
        omp_set_num_threads(atoi(argv[6]));
    }
    printf("cfgFile=%s\n", cfgFile);

    sprintf(outputFilePath, "aramco-reb%s-aeb%s.b%s.stack.out", argv[4], argv[5], argv[3]);

    size_t nbEle = 449 * 449 * 235;
    size_t allocBufSize = nbEle * 300;
    float buf[allocBufSize];
    float stack[nbEle];
    std::fill_n(stack, nbEle, 0);

    std::ifstream is(oriFilePath);
    size_t snapshot = 0;
    do {
        is.read((char *) &buf[0], sizeof(float) * allocBufSize);
        size_t bufSize = is.gcount() / sizeof(float);

        for (size_t l = 0; l < bufSize; l += nbEle) {
            printf("\n\n************* snapshot %04lu ***********\n", snapshot++);
            float *data = buf + l;
            float max = data[0], min = data[0];
            for (size_t i = 0; i < nbEle; i++) {
                if (data[i] > max) {
                    max = data[i];
                }
                if (data[i] < min) {
                    min = data[i];
                }
            }
            float errBound = (max - min) * reb;
            if (errBound > abs) {
                errBound = abs;
            }

            size_t outSize;
            int status = SZ_Init(cfgFile);
            if (status == SZ_NSCS)
                exit(0);
            cost_start();
#ifdef _OPENMP
            unsigned char *bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_float_openmp(
                    data, &outSize, errBound, nbEle, blockSize);
#else
            unsigned char* bytes = SZ_fast_compress_args_unpredictable_blocked_randomaccess_float(data, &outSize, errBound, nbEle, blockSize);
#endif

            cost_end();
            printf("\ntimecost=%f, total fastmode2\n", totalCost);
            printf("compression size = %zu, CR = %f\n", outSize, 1.0f * nbEle * sizeof(float) / outSize);

            SZ_Finalize();


            float *dec_data;
#ifdef _OPENMP
            SZ_fast_decompress_args_unpredictable_blocked_randomaccess_float_openmp(&dec_data, nbEle, bytes);
#else
            SZ_fast_decompress_args_unpredictable_blocked_randomaccess_float(&data, nbEle, bytes);
#endif
            cost_end();

            free(bytes);
            printf("timecost=%f\n", totalCost);

            size_t i = 0;
            float Max = 0, Min = 0, diffMax = 0;
            Max = data[0];
            Min = data[0];
            diffMax = fabs(dec_data[0] - data[0]);
            double sum1 = 0, sum2 = 0;
            for (i = 0; i < nbEle; i++) {
                sum1 += data[i];
                sum2 += dec_data[i];
            }
            double mean1 = sum1 / nbEle;
            double mean2 = sum2 / nbEle;

            double sum3 = 0, sum4 = 0;
            double sum = 0, prodSum = 0, relerr = 0;

            double maxpw_relerr = 0;
            for (i = 0; i < nbEle; i++) {
                stack[i] += dec_data[i];
                if (Max < data[i]) Max = data[i];
                if (Min > data[i]) Min = data[i];

                float err = fabs(dec_data[i] - data[i]);
                if (data[i] != 0) {
                    if (fabs(data[i]) > 1)
                        relerr = err / data[i];
                    else
                        relerr = err;
                    if (maxpw_relerr < relerr)
                        maxpw_relerr = relerr;
                }

                if (diffMax < err)
                    diffMax = err;
                prodSum += (data[i] - mean1) * (dec_data[i] - mean2);
                sum3 += (data[i] - mean1) * (data[i] - mean1);
                sum4 += (dec_data[i] - mean2) * (dec_data[i] - mean2);
                sum += err * err;
            }
            double std1 = sqrt(sum3 / nbEle);
            double std2 = sqrt(sum4 / nbEle);
            double ee = prodSum / nbEle;
            double acEff = ee / std1 / std2;

            double mse = sum / nbEle;
            double range = Max - Min;
            double psnr = 20 * log10(range) - 10 * log10(mse);
            double nrmse = sqrt(mse) / range;

            double compressionRatio = 1.0 * nbEle * sizeof(float) / outSize;

            printf("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
            printf("Max absolute error = %.10f\n", diffMax);
            printf("Max relative error = %f\n", diffMax / (Max - Min));
            printf("Max pw relative error = %f\n", maxpw_relerr);
            printf("PSNR = %f, NRMSE= %.20G\n", psnr, nrmse);
            printf("acEff=%f\n", acEff);
            printf("compressionRatio = %f\n", compressionRatio);

            free(dec_data);

        }
    } while (is);

    writefile(outputFilePath, stack, nbEle);
    return 0;
}
