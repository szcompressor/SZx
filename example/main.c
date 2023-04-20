#include <stdio.h>
#include "mpi.h"
#include "recursive_doubling.h"
#include "ring.h"
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
// #include "./utils.h"

#define ITERATIONS_LARGE 100
#define LARGE_MESSAGE_SIZE 1024 * 1024 // This is in bytes
#define MIN_MESSAGE_LENGTH 1           // This is in length
// #define compressionRatio 20
#define tolerance 0.08
#define chunk_size 10240
// #define data_size 520000
typedef float data_type;

int main()
{
    int validation = 1;
    int i = 0;
    int provided;
    double compressionRatio = 1E-3;
    double absErrBound = compressionRatio;
    // MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
    MPI_Init_thread(NULL, NULL, MPI_THREAD_SINGLE, &provided);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // data_type *invec = NULL;
    // srand(time(NULL));
    // invec = create_rand_nums(size);
    // int size = 10000;
    // invec = create_rand_nums(size);

    // if (validation && i == 0)
    // {
    // Allocate memory for the testing array
    // data_type *our_array = NULL;
    // our_array = inilize_arr(size);
    // data_type *MPI_array = NULL;
    // MPI_array = inilize_arr(size);
    // MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Allreduce(invec, MPI_array, size, MPI_FLOAT,
    //               MPI_SUM, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);
    // MPIR_Allreduce_intra_ring(invec, our_array, size, MPI_FLOAT,
    //                           MPI_SUM, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);
    // if (world_rank == 0)
    // {
    //     printf("Verification begins\n");
    // }
    // if (!verify_arrays(our_array, MPI_array, size))
    // {
    //     printf("Oops!\n");
    //     exit(1);
    // }
    // else if (world_rank == 0)
    // {
    //     printf("Verification passed!\n");
    // }

    int extent = sizeof(MPI_FLOAT);
    MPI_Request reqs[2]; /* one send and one recv per transfer */
    MPI_Status stas[2];

    char oriFilePath[645] = "/lcrc/project/ECP-EZ/shdi/RtmLab/examples/overthrust_model/data/snpData_0600_849X849X235.dat";
    int status = 0;
    size_t nbEle;
    float *recvbuf = readFloatData(oriFilePath, &nbEle, &status);
    int data_size = nbEle;

    int outSize = 0;
    int chunk_num = data_size / chunk_size;
    int chunk_remainder_size = data_size % chunk_size;
    printf("data_size = %d chunk_num = %d chunk_reminder_size = %d\n", data_size, chunk_num, chunk_remainder_size);

    int flag, iter;
    size_t chunk_out_size = 0;
    unsigned char *outputBytes = (unsigned char *)malloc(data_size * extent + sizeof(size_t) * (chunk_num + 1));
    float *newData = (float *)malloc(data_size * extent);
    int first_out_size, second_out_size, third_out_size, sum_out_size = 0;
    int sum_iter = 12;
    for (iter = 0; iter < chunk_num; iter++)
    {
        // MPI_Test(&reqs[0], &flag, &stas[0]);
        // SZ_fast_compress_args2(fast_mode, SZ_FLOAT, (char *)recvbuf + iter * chunk_size * extent,
        //                        &chunk_out_size, outputBytes + outSize,
        //                        cpr_mode, absErrBound, relBoundRatio, compressionRatio, tolerance, 0, 0, 0, 0, chunk_size);
        SZ_fast_compress_args2_split(fast_mode, SZ_FLOAT, (char *)recvbuf + iter * chunk_size * extent,
                                     &chunk_out_size, outputBytes + outSize + sizeof(size_t) * (chunk_num + 1), outputBytes, iter,
                                     cpr_mode, absErrBound, relBoundRatio, compressionRatio, tolerance, 0, 0, 0, 0, chunk_size);
        outSize += chunk_out_size;
        // if (iter == 0)
        // {
        //     first_out_size = chunk_out_size;
        // }
        // if (iter == 1)
        // {
        //     second_out_size = chunk_out_size;
        // }
        // if (iter == 2)
        // {
        //     third_out_size = chunk_out_size;
        // }
        if (iter < sum_iter)
        {
            sum_out_size += chunk_out_size;
            printf("At iter: %d chunk_out_size = %d\n", iter, chunk_out_size);
        }
        // printf("outSize = %d\n", outSize);
    }
    printf("outSize = %d\n", outSize);
    printf("Main iterations finished\n");

    // Handling the remainder size
    if (chunk_remainder_size != 0)
    {
        // MPI_Test(&reqs[0], &flag, &stas[0]);
        // SZ_fast_compress_args2(fast_mode, SZ_FLOAT, (char *)recvbuf + (data_size - chunk_remainder_size) * extent,
        //                        &chunk_out_size, outputBytes + outSize,
        //                        cpr_mode, absErrBound, relBoundRatio, compressionRatio, tolerance, 0, 0, 0, 0, chunk_remainder_size);
        SZ_fast_compress_args2_split(fast_mode, SZ_FLOAT, (char *)recvbuf + iter * chunk_size * extent,
                                     &chunk_out_size, outputBytes + outSize + sizeof(size_t) * (chunk_num + 1), outputBytes, iter,
                                     cpr_mode, absErrBound, relBoundRatio, compressionRatio, tolerance, 0, 0, 0, 0, chunk_remainder_size);
        outSize += chunk_out_size;
        printf("outSize = %d\n", outSize);
    }
    printf("Remainder iterations finished\n");

    size_t *bytes = (size_t *)outputBytes;
    // unsigned char *first_data = SZ_fast_decompress(fast_mode, SZ_FLOAT, bytes, first_out_size/extent, 0, 0, 0, 0, chunk_size);
    // unsigned char *second_data = SZ_fast_decompress(fast_mode, SZ_FLOAT, bytes+first_out_size, second_out_size/extent, 0, 0, 0, 0, chunk_size);
    // unsigned char *sum_data = SZ_fast_decompress(fast_mode, SZ_FLOAT, bytes, (sum_out_size) / extent, 0, 0, 0, 0, sum_iter * chunk_size);
    // unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, bytes, outSize/extent, 0, 0, 0, 0, (data_size - chunk_remainder_size));
    // unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, bytes, outSize/extent, 0, 0, 0, 0, data_size);

    outSize = 0;
    chunk_num = data_size / chunk_size;
    chunk_remainder_size = data_size % chunk_size;
    chunk_out_size = 0;
    printf("data_size = %d chunk_num = %d chunk_reminder_size = %d\n", data_size, chunk_num, chunk_remainder_size);
    int decom_offset = 0;
    int decom_out_offset = 0;
    // for (iter = 0; iter < 3; iter++)
    // {
    //     printf("iter:%d Compressed size:%d Compressed size address:%d\n", iter, (int)bytes[iter], &bytes[iter]);
    //     unsigned char *cmpBytes = (unsigned char *)bytes + sizeof(size_t) * (chunk_num + 1) + decom_offset;
    //     unsigned char *r = (unsigned char *)cmpBytes;
    //     r += 4;
    //     int blockSize = r[0]; // get block size
    //     printf("blockSize = %d Block size address: %d\n", r[0], &r[0]);
    //     printf("cmpBytes address = %d\n", cmpBytes);
    //     decom_offset += bytes[iter];
    //     decom_out_offset += chunk_size * extent;
    // }

    for (iter = 0; iter < chunk_num; iter++)
    {
        SZ_fast_decompress_split(fast_mode, SZ_FLOAT, (unsigned char *)newData + decom_out_offset, (unsigned char *)bytes + sizeof(size_t) * (chunk_num + 1) + decom_offset, 0, 0, 0, 0, chunk_size);
        decom_offset += bytes[iter];
        decom_out_offset += chunk_size * extent;
        // printf("iter%d\n",iter);
    }
    // Handling the remainder size
    if (chunk_remainder_size != 0)
    {
        // MPI_Test(&reqs[0], &flag, &stas[0]);
        SZ_fast_decompress_split(fast_mode, SZ_FLOAT, (unsigned char *)newData + decom_out_offset, (unsigned char *)bytes + sizeof(size_t) * (chunk_num + 1) + decom_offset, 0, 0, 0, 0, chunk_remainder_size);
    }

    free(newData);
    // free(data);
    free(outputBytes);
    free(recvbuf);

    // free(our_array);
    // free(MPI_array);
    MPI_Barrier(MPI_COMM_WORLD);
    // }
}
