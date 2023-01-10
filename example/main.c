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
#define compressionRatio 20
#define tolerance 0.08
typedef float data_type;

int main()
{
    int validation = 1;
    int i = 0;
    int provided;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    data_type *invec = NULL;
    // srand(time(NULL));
    // invec = create_rand_nums(size);
    int size = 10000;
    invec = create_rand_nums(size);

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
        MPIR_Allreduce_intra_ring(invec, our_array, size, MPI_FLOAT,
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
}
