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
#include "./include/libs.h"
// #define data_size 520000
typedef float data_type;
#define MPI_THREAD_MODE MPI_THREAD_FUNNELED

int main()
{
    int validation = 1;
    int i = 0;
    int provided;
    double compressionRatio = 1E-3;
    double absErrBound = compressionRatio;
    // MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MODE, &provided);
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
    int flag, iter, mpi_errno;
    size_t chunk_out_size = 0;
    unsigned char *outputBytes = (unsigned char *)malloc(data_size * extent);
    float *newData = (float *)malloc(data_size * extent);
    unsigned char *rbuf = (unsigned char *)malloc(data_size * extent);
    float *sbuf = (float *)malloc(data_size * extent);
    double CPR_timer = 0.0, MPI_timer = 0.0;
    unsigned char *bytes;
    for (iter = 0; iter < 10; iter++)
    {
        SZ_fast_compress_args2(2, SZ_FLOAT, recvbuf, &outSize, outputBytes, cpr_mode, absErrBound, relBoundRatio,
                               compressionRatio, tolerance, 0, 0, 0, 0, data_size);
        bytes = outputBytes;

        // mpi_errno = MPI_Sendrecv((void *)bytes, outSize / sizeof(MPI_FLOAT), MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
        //                          rbuf, data_size, MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
        //                          MPI_COMM_WORLD, &status);
        mpi_errno = MPI_Sendrecv(recvbuf, data_size/20, MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
                                 rbuf, data_size/20, MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
                                 MPI_COMM_WORLD, &status);
        if (mpi_errno)
        {
            exit(-1);
        }
    }
    for (iter = 0; iter < 10; iter++)
    {
        CPR_timer -= MPI_Wtime();
        SZ_fast_compress_args2(2, SZ_FLOAT, recvbuf, &outSize, outputBytes, cpr_mode, absErrBound, relBoundRatio,
                               compressionRatio, tolerance, 0, 0, 0, 0, data_size);
        bytes = outputBytes;
        CPR_timer += MPI_Wtime();
        MPI_timer -= MPI_Wtime();
        // mpi_errno = MPI_Sendrecv((void *)bytes, outSize / sizeof(MPI_FLOAT), MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
        //                          rbuf, data_size, MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
        //                          MPI_COMM_WORLD, &status);
        mpi_errno = MPI_Sendrecv(recvbuf, data_size/20, MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
                                 rbuf, data_size/20, MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
                                 MPI_COMM_WORLD, &status);
        if (mpi_errno)
        {
            exit(-1);
        }
        MPI_timer += MPI_Wtime();
    }

    printf("For process %d Single thread comm time is %f, cpr time is %f\n", world_rank, MPI_timer, CPR_timer);

    MPI_timer = 0.0;
    CPR_timer = 0.0;
     for (iter = 0; iter < 10; iter++)
    {
        SZ_fast_compress_args2(4, SZ_FLOAT, recvbuf, &outSize, outputBytes, cpr_mode, absErrBound, relBoundRatio,
                               compressionRatio, tolerance, 0, 0, 0, 0, data_size);
        bytes = outputBytes;

        // mpi_errno = MPI_Sendrecv((void *)bytes, outSize / sizeof(MPI_FLOAT), MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
        //                          rbuf, data_size, MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
        //                          MPI_COMM_WORLD, &status);
        mpi_errno = MPI_Sendrecv(recvbuf, data_size/20, MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
                                 rbuf, data_size/20, MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
                                 MPI_COMM_WORLD, &status);
        if (mpi_errno)
        {
            exit(-1);
        }
    }

    for (iter = 0; iter < 10; iter++)
    {
        CPR_timer -= MPI_Wtime();
        SZ_fast_compress_args2(4, SZ_FLOAT, recvbuf, &outSize, outputBytes, cpr_mode, absErrBound, relBoundRatio,
                               compressionRatio, tolerance, 0, 0, 0, 0, data_size);
        bytes = outputBytes;
        CPR_timer += MPI_Wtime();
        MPI_timer -= MPI_Wtime();
        // mpi_errno = MPI_Sendrecv((void *)bytes, outSize / sizeof(MPI_FLOAT), MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
        //                          rbuf, data_size, MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
        //                          MPI_COMM_WORLD, &status);
        mpi_errno = MPI_Sendrecv(recvbuf, data_size/20, MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
                                 rbuf, data_size/20, MPI_FLOAT, !world_rank, MPIR_ALLGATHERV_TAG,
                                 MPI_COMM_WORLD, &status);
        if (mpi_errno)
        {
            exit(-1);
        }
        MPI_timer += MPI_Wtime();
    }
    printf("For process %d Multi thread comm time is %f, cpr time is %f\n", world_rank, MPI_timer, CPR_timer);
    free(newData);
    // free(data);
    free(outputBytes);
    free(recvbuf);
    free(sbuf);
    free(rbuf);
    // free(our_array);
    // free(MPI_array);
    MPI_Barrier(MPI_COMM_WORLD);
    // }
}
