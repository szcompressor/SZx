#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "szx.h"
#include "szx_rw.h"

double calculate_mse(float *ori_data, float *dec_data, size_t nbEle)
{
    size_t i = 0;
    double Max = 0, Min = 0, diffMax = 0;
    Max = ori_data[0];
    Min = ori_data[0];
    diffMax = dec_data[0] > ori_data[0] ? dec_data[0] - ori_data[0] : ori_data[0] - dec_data[0];

    // diffMax = fabs(data[0] - ori_data[0]);
    double sum1 = 0, sum2 = 0, sum22 = 0;

    for (i = 0; i < nbEle; i++)
    {
        sum1 += ori_data[i];
        sum2 += dec_data[i];
        sum22 += dec_data[i] * dec_data[i];
    }
    double mean1 = sum1 / nbEle;
    double mean2 = sum2 / nbEle;

    double sum3 = 0, sum4 = 0;
    double sum = 0, prodSum = 0, relerr = 0;

    double maxpw_relerr = 0;
    for (i = 0; i < nbEle; i++)
    {
        if (Max < ori_data[i])
            Max = ori_data[i];
        if (Min > ori_data[i])
            Min = ori_data[i];
        float err_ = dec_data[i] - ori_data[i];
        // if (printError)
        //     errors[i] = err_;
        float err = fabs(err_);
        if (ori_data[i] != 0)
        {
            relerr = err / fabs(ori_data[i]);
            if (maxpw_relerr < relerr)
                maxpw_relerr = relerr;
        }

        if (diffMax < err)
            diffMax = err;
        prodSum += (ori_data[i] - mean1) * (dec_data[i] - mean2);
        sum3 += (ori_data[i] - mean1) * (ori_data[i] - mean1);
        sum4 += (dec_data[i] - mean2) * (dec_data[i] - mean2);
        sum += err * err;
    }
    double std1 = sqrt(sum3 / nbEle);
    double std2 = sqrt(sum4 / nbEle);
    double ee = prodSum / nbEle;
    double acEff = ee / std1 / std2;

    double mse = sum / nbEle;
    double range = Max - Min;
    // printf("the range is %.20G, the mse is %.20G\n", range, mse);
    double psnr = 20 * log10(range) - 10 * log10(mse);
    double normErr = sqrt(sum);
    double normErr_norm = normErr / sqrt(sum22);
    double nrmse = sqrt(mse) / range;
    return mse;
}

int main(int argc, char *argv[])
{
    int rank, iteration;
    size_t nbEle, decomnbEle;
    

    printf("Reduce-scatter start\n");
    for (iteration = 0; iteration < 15; iteration++)
    {
        for (rank = 0; rank < 16; rank++)
        {
            int status;
            char originalFilePath[256];
            char decomFilePath[256];
            char cmpPath[256] = "/lcrc/project/sbi-fair/jiajun/MPI-Coll-SZx-data/16_allreduce_midsteps";
            sprintf(originalFilePath, "%s/%d_%d.reduce-scatter", cmpPath, rank, iteration);
            sprintf(decomFilePath, "%s/zfp_FXR%d_%d.reduce-scatter", cmpPath, rank, iteration);
            float *ori_data = readFloatData(originalFilePath, &nbEle, &status);
            float *dec_data = readFloatData(decomFilePath, &decomnbEle, &status);
            if (nbEle != decomnbEle)
            {
                printf("Error: number of elements is not consistent\n");
                exit(0);
            }
            double mse = calculate_mse(ori_data, dec_data, nbEle);
            // printf("for iteration %d, rank %d the mse is %.20G, \n", iteration, rank, mse);
            printf("%.20G, ", mse);
        }
        printf("\n");
    }

    printf("Allgather start\n");
    for (iteration = 0; iteration < 15; iteration++)
    {
        for (rank = 0; rank < 16; rank++)
        {
            int status;
            char originalFilePath[256];
            char decomFilePath[256];
            char cmpPath[256] = "/lcrc/project/sbi-fair/jiajun/MPI-Coll-SZx-data/16_allreduce_midsteps";
            sprintf(originalFilePath, "%s/%d_%d.allagther", cmpPath, rank, iteration);
            sprintf(decomFilePath, "%s/zfp_FXR%d_%d.allagther", cmpPath, rank, iteration);
            float *ori_data = readFloatData(originalFilePath, &nbEle, &status);
            float *dec_data = readFloatData(decomFilePath, &decomnbEle, &status);
            // printf("\nFinished_reading\n");
            if (nbEle != decomnbEle)
            {
                printf("Error: number of elements is not consistent\n");
                exit(0);
            }
            double mse = calculate_mse(ori_data, dec_data, nbEle);
            // printf("for iteration %d, rank %d the mse is %f, \n", iteration, rank, mse);
            printf("%.20G, ", mse);
        }
        printf("\n");
    }
}
