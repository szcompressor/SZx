#include <stdio.h>
#include "mpi.h"
#include "recursive_doubling.h"
#include "ring.h"
#include "ring2.h"
#include "ring2_overlap.h"
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
// #include "./utils.h"

#define ITERATIONS_LARGE 100
#define LARGE_MESSAGE_SIZE 1024 * 1024 // This is in bytes
#define MIN_MESSAGE_LENGTH 1           // This is in length
#define compressionRatio 20
#define tolerance 0.08
#define MPI_THREAD_MODE MPI_THREAD_SINGLE
// #define MPI_THREAD_MODE MPI_THREAD_FUNNELED
typedef float data_type;

int main(int argc, char *argv[])
{
        // multiple proce on one node, multi-node from one byte to 4 MB OSU benchmarks
        int provided;
        MPI_Init_thread(NULL, NULL, MPI_THREAD_MODE, &provided);
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        // check if the process size is large enough
        if (world_size < 2)
        {
                if (world_rank == 0)
                {
                        fprintf(stderr, "This test requires at least two processes\n");
                }
                MPI_Finalize();
                exit(EXIT_FAILURE);
        }
        if (world_rank == 0)
        {
                printf("Welcome to our ANL_allreduce_benchmark\n");
        }
        int opt;
        int warm_up = 0;
        int num_trials = 1000;
        int validation = 0;
        int select = 1;
        int minimal_size = 4 / sizeof(data_type);
        int maximal_size = 4 * 1024 * 1024 / sizeof(data_type);
        int large_size = LARGE_MESSAGE_SIZE / sizeof(data_type);
        while ((opt = getopt(argc, argv, "i:w:v:s:l:k:")) != EOF)
        {
                switch (opt)
                {
                case 'i':
                        num_trials = atoi(optarg);
                        break;
                case 'w':
                        warm_up = atoi(optarg);
                        break;
                case 'v':
                        validation = atoi(optarg);
                        break;
                case 's':
                        minimal_size = atoi(optarg) / sizeof(data_type);
                        break;
                case 'l':
                        maximal_size = atoi(optarg) / sizeof(data_type);
                        break;
                case 'k':
                        select = atoi(optarg);
                        break;
                case '?':
                        if (world_rank == 0)
                        {
                                printf("usage is: ./ANL_benchmarks\n"
                                       "-i <number of iterations> the default value is 1000\n"
                                       "-w <number of warmups> the default value is 0\n"
                                       "-v <enable validation> the default value is 0\n"
                                       "-s <minimal data size in bytes> the default value is 4 bytes\n"
                                       "-l <maximal data size in bytes> the default value is 4 MB\n"
                                       "-k <select the kernel> 0 for original allreduce, 1 for our allreduce, 2 for SZx allreduce default is 1\n"
                                       "-? printf this message\n");
                                break;
                        }

                default:
                        exit(1);
                }
        }        
        // char oriFilePath[645];
        //ML files
        // char oriFileDire[640] = "/lcrc/project/sbi-fair/shared/gradients_saved/";
        // sprintf(oriFilePath, "%srank0_step_%d_layer_15.bin", oriFileDire, world_rank+50);
        //small aramco files 47376235 189504940 bytes
        // char oriFileDire[645] = "/lcrc/project/ECP-EZ/shdi/RtmLab-small/examples/overthrust_model/data/";
        // sprintf(oriFilePath, "%saramco-snapshot-00%d.f32", oriFileDire, world_rank+10);
        // char oriFilePath[645] = "/lcrc/project/ECP-EZ/shdi/RtmLab-small/examples/overthrust_model/data/aramco-snapshot-0001.f32";
        //big snpData files 849X849X235=169388235 677552940 bytes
        // char oriFileDire[640] = "/lcrc/project/ECP-EZ/shdi/RtmLab/examples/overthrust_model/data/";
        // sprintf(oriFilePath, "%ssnpData_0%d_849X849X235.dat", oriFileDire, world_rank*50+100);
        char oriFilePath[645] = "/lcrc/project/ECP-EZ/shdi/RtmLab/examples/overthrust_model/data/snpData_0020_849X849X235.dat";

        int status = 0;
        size_t nbEle;
        float *numbers = readFloatData(oriFilePath, &nbEle, &status);
        if (status != SZ_SCES)
        {
                printf("Error: data file %s cannot be read!\n", oriFilePath);
                exit(0);
        }
        if (!world_rank)
                printf("Original file size: %d\n", nbEle);
        if (num_trials <= 0)
        {
                printf("Please select a valid number of iterations.\n");
                exit(1);
        }
        if (warm_up < 0)
        {
                printf("Please select a valid number of warm_up.\n");
                exit(1);
        }
        if (validation != 0 && validation != 1)
        {
                printf("Please select a valid status of validation.\n");
                exit(1);
        }
        if (minimal_size < MIN_MESSAGE_LENGTH)
                minimal_size = MIN_MESSAGE_LENGTH;
        if (maximal_size < MIN_MESSAGE_LENGTH)
                maximal_size = MIN_MESSAGE_LENGTH;
        if (maximal_size > nbEle)
                maximal_size = nbEle;
        // if (select != 0 && select != 1 && select != 2 && select != 3 && select != 4)
        // {
        //         printf("Please select a valid kernel.\n");
        //         exit(1);
        // }
        if (world_rank == 0)
        {
                printf("The settings are: %d iterations, %d warmups, validation: %d, "
                       "minimal data size: %ld bytes, maximal data size: %ld bytes, kernel: %d\n",
                       num_trials, warm_up, validation, minimal_size * sizeof(data_type), maximal_size * sizeof(data_type), select);
        }

        int size, iterations, i;

        for (size = minimal_size; size <=
                                  maximal_size;
             size *= 2)
        {
                iterations = num_trials;
                // if (size > large_size)
                // {
                //         iterations = ITERATIONS_LARGE;
                // }
                // Allocate memory for the input array
                data_type *invec = NULL;
                // srand(time(NULL));
                // invec = create_rand_nums(size);
                invec = inilize_arr_withoutset(size);
                memcpy(invec, numbers, sizeof(data_type) * size);
                // Allocate memory for the output array
                data_type *inoutvec = NULL;
                srand(time(NULL));
                inoutvec = inilize_arr(size);

                MPI_Barrier(MPI_COMM_WORLD);

                for (i = 0; i < warm_up; i++)
                {
                        if (select == 0)
                        {
                                MPI_Allreduce(invec, inoutvec, size, MPI_FLOAT, MPI_SUM,
                                              MPI_COMM_WORLD);
                        }
                        else if (select == 1)
                        {
                                MPI_Allreduce_compre_RB(invec, inoutvec, size, MPI_FLOAT, MPI_SUM,
                                                     MPI_COMM_WORLD);
                        }
                        else if (select == 2)
                        {
                                MPI_Allreduce_SZx_FXR_RB(invec, inoutvec, compressionRatio, tolerance, size, MPI_FLOAT, MPI_SUM,
                                                  MPI_COMM_WORLD);
                        }
                        else if (select == 3)
                        {
                                MPI_Allreduce_SZx_FXR_RB_record(invec, inoutvec, compressionRatio, tolerance, size, MPI_FLOAT, MPI_SUM,
                                                      MPI_COMM_WORLD);
                        }
                        else if (select == 4)
                        {
                                MPIR_Allreduce_intra_ring(invec, inoutvec, size, MPI_FLOAT, MPI_SUM,
                                                      MPI_COMM_WORLD);
                        }
                        else if (select == 5)
                        {
                                MPI_Allreduce_SZx_FXR_RI2(invec, inoutvec, compressionRatio, tolerance, size, MPI_FLOAT, MPI_SUM,
                                                      MPI_COMM_WORLD);
                        }
                        else if (select == 6)
                        {
                                MPI_Allreduce_SZx_FXR_RI2_record(invec, inoutvec, compressionRatio, tolerance, size, MPI_FLOAT, MPI_SUM,
                                                      MPI_COMM_WORLD);
                        }
                }
                MPI_Barrier(MPI_COMM_WORLD);
                double MPI_timer = 0.0;
                for (i = 0; i < iterations; i++)
                {
                        if (validation && i == 0)
                        {
                                // Allocate memory for the testing array
                                data_type *our_array = NULL;
                                our_array = inilize_arr(size);
                                data_type *MPI_array = NULL;
                                MPI_array = inilize_arr(size);
                                MPI_Barrier(MPI_COMM_WORLD);
                                MPI_Allreduce(invec, MPI_array, size, MPI_FLOAT,
                                              MPI_SUM, MPI_COMM_WORLD);
                                MPI_Barrier(MPI_COMM_WORLD);
                                MPIR_Allreduce_intra_ring(invec, inoutvec, size, MPI_FLOAT,
                                              MPI_SUM, MPI_COMM_WORLD);
                                MPI_Barrier(MPI_COMM_WORLD);
                                if (world_rank == 0)
                                {
                                        printf("Verification begins\n");
                                }
                                if (!verify_arrays(our_array, MPI_array, size))
                                {
                                        printf("Oops!\n");
                                        exit(1);
                                }
                                else if (world_rank == 0)
                                {
                                        printf("Verification passed!\n");
                                }
                                free(our_array);
                                free(MPI_array);
                                MPI_Barrier(MPI_COMM_WORLD);
                        }
                        if (select == 0)
                        {
                                MPI_Barrier(MPI_COMM_WORLD);
                                MPI_timer -= MPI_Wtime();
                                MPI_Allreduce(invec, inoutvec, size, MPI_FLOAT, MPI_SUM,
                                              MPI_COMM_WORLD);
                                MPI_timer += MPI_Wtime();
                        }
                        else if (select == 1)
                        {
                                MPI_Barrier(MPI_COMM_WORLD);
                                MPI_timer -= MPI_Wtime();
                                MPI_Allreduce_compre_RB(invec, inoutvec, size, MPI_FLOAT, MPI_SUM,
                                                     MPI_COMM_WORLD);
                                MPI_timer += MPI_Wtime();
                        }
                        else if (select == 2)
                        {
                                MPI_Barrier(MPI_COMM_WORLD);
                                MPI_timer -= MPI_Wtime();
                                MPI_Allreduce_SZx_FXR_RB(invec, inoutvec, compressionRatio, tolerance, size, MPI_FLOAT, MPI_SUM,
                                                  MPI_COMM_WORLD);
                                MPI_timer += MPI_Wtime();
                        }
                        else if (select == 3)
                        {
                                MPI_Barrier(MPI_COMM_WORLD);
                                MPI_timer -= MPI_Wtime();
                                MPI_Allreduce_SZx_FXR_RB_record(invec, inoutvec, compressionRatio, tolerance, size, MPI_FLOAT, MPI_SUM,
                                                      MPI_COMM_WORLD);
                                MPI_timer += MPI_Wtime();
                        }
                        else if (select == 4)
                        {
                                MPI_Barrier(MPI_COMM_WORLD);
                                MPI_timer -= MPI_Wtime();
                                MPIR_Allreduce_intra_ring(invec, inoutvec, size, MPI_FLOAT, MPI_SUM,
                                                      MPI_COMM_WORLD);
                                MPI_timer += MPI_Wtime();
                        }
                        else if (select == 5)
                        {
                                MPI_Barrier(MPI_COMM_WORLD);
                                MPI_timer -= MPI_Wtime();
                                MPI_Allreduce_SZx_FXR_RI2(invec, inoutvec, compressionRatio, tolerance, size, MPI_FLOAT, MPI_SUM,
                                                      MPI_COMM_WORLD);
                                MPI_timer += MPI_Wtime();
                        }
                        else if (select == 6)
                        {
                                MPI_Barrier(MPI_COMM_WORLD);
                                MPI_timer -= MPI_Wtime();
                                MPI_Allreduce_SZx_FXR_RI2_record(invec, inoutvec, compressionRatio, tolerance, size, MPI_FLOAT, MPI_SUM,
                                                      MPI_COMM_WORLD);
                                MPI_timer += MPI_Wtime();
                        }
                        else if (select == 7)
                        {
                                MPI_Barrier(MPI_COMM_WORLD);
                                MPI_timer -= MPI_Wtime();
                                MPI_Allreduce_SZx_FXR_RI2_op_record(invec, inoutvec, compressionRatio, tolerance, size, MPI_FLOAT, MPI_SUM,
                                                      MPI_COMM_WORLD);
                                MPI_timer += MPI_Wtime();
                        }
                        
                        MPI_Barrier(MPI_COMM_WORLD);
                }
                double latency = (double)(MPI_timer * 1e6) / iterations;
                double min_time = 0.0;
                double max_time = 0.0;
                double avg_time = 0.0;
                MPI_Reduce(&latency, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0,
                           MPI_COMM_WORLD);
                MPI_Reduce(&latency, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                           MPI_COMM_WORLD);
                MPI_Reduce(&latency, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0,
                           MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
                avg_time = avg_time / world_size;
                if (world_rank == 0)
                {
                        printf("Routine:%d For datasize: %ld bytes, the avg_time is %f us, the max_time is %f us, the min_time is %f us\n",
                               select, size * sizeof(data_type), avg_time, max_time, min_time);
                        if (select == 0)
                        {
                                FILE *myWrite = fopen("data/allreduce_mpi.txt", "a");
                                if (myWrite == NULL)
                                {
                                        return 0;
                                }
                                fprintf(myWrite, "%f, ", avg_time);
                                fclose(myWrite);
                        }
                        if (select == 1)
                        {
                                FILE *myWrite = fopen("data/allreduce_com.txt", "a");
                                if (myWrite == NULL)
                                {
                                        return 0;
                                }
                                fprintf(myWrite, "%f, ", avg_time);
                                fclose(myWrite);
                        }
                        if (select == 2)
                        {
                                FILE *myWrite = fopen("data/allreduce_szx_fxr.txt", "a");
                                if (myWrite == NULL)
                                {
                                        return 0;
                                }
                                fprintf(myWrite, "%f, ", avg_time);
                                fclose(myWrite);
                        }
                }
                MPI_Barrier(MPI_COMM_WORLD);
                free(invec);
                free(inoutvec);
                MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
}
