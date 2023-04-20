/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */
// This is for the binomial tree broadcast
#include <stdio.h>
#include "mpi.h"
// #include "./utils.h"
#include "szx.h"
#include "szx_rw.h"

#include "./include/libs.h"

/* Algorithm: Binomial bcast
 *
 * For short messages, we use a binomial tree algorithm.
 * Cost = lgp.alpha + n.lgp.beta
 */
int MPIR_Bcast_intra_binomial(void *buffer,
                              MPI_Aint count,
                              MPI_Datatype datatype,
                              int root, MPI_Comm comm)
{
    double MPI_total_timer = 0.0;
    MPI_total_timer -= MPI_Wtime();
    int rank, comm_size, src, dst;
    int relative_rank, mask;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;
    MPI_Aint nbytes = 0;
    MPI_Status status;

    int is_contig;
    MPI_Aint type_size;
    void *tmp_buf = NULL;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    type_size = sizeof(data_type);

    nbytes = type_size * count;
    if (nbytes == 0)
        goto fn_exit; /* nothing to do */

    relative_rank = (rank >= root) ? rank - root : rank - root + comm_size;

    /* Use short message algorithm, namely, binomial tree */

    /* Algorithm:
     * This uses a fairly basic recursive subdivision algorithm.
     * The root sends to the process comm_size/2 away; the receiver becomes
     * a root for a subtree and applies the same process.
     *
     * So that the new root can easily identify the size of its
     * subtree, the (subtree) roots are all powers of two (relative
     * to the root) If m = the first power of 2 such that 2^m >= the
     * size of the communicator, then the subtree at root at 2^(m-k)
     * has size 2^k (with special handling for subtrees that aren't
     * a power of two in size).
     *
     * Do subdivision.  There are two phases:
     * 1. Wait for arrival of data.  Because of the power of two nature
     * of the subtree roots, the source of this message is always the
     * process whose relative rank has the least significant 1 bit CLEARED.
     * That is, process 4 (100) receives from process 0, process 7 (111)
     * from process 6 (110), etc.
     * 2. Forward to my subtree
     *
     * Note that the process that is the tree root is handled automatically
     * by this code, since it has no bits set.  */

    mask = 0x1;
    while (mask < comm_size)
    {
        if (relative_rank & mask)
        {
            src = rank - mask;
            if (src < 0)
                src += comm_size;
            mpi_errno = MPI_Recv(buffer, count, datatype, src,
                                 MPIR_BCAST_TAG, comm, &status);
            if (mpi_errno)
            {
                exit(-1);
            }
            break;
        }
        mask <<= 1;
    }

    /* This process is responsible for all processes that have bits
     * set from the LSB up to (but not including) mask.  Because of
     * the "not including", we start by shifting mask back down one.
     *
     * We can easily change to a different algorithm at any power of two
     * by changing the test (mask > 1) to (mask > block_size)
     *
     * One such version would use non-blocking operations for the last 2-4
     * steps (this also bounds the number of MPI_Requests that would
     * be needed).  */

    mask >>= 1;
    while (mask > 0)
    {
        if (relative_rank + mask < comm_size)
        {
            dst = rank + mask;
            if (dst >= comm_size)
                dst -= comm_size;
            mpi_errno = MPI_Send(buffer, count, datatype, dst,
                                 MPIR_BCAST_TAG, comm);
            if (mpi_errno)
            {
                exit(-1);
            }
        }
        mask >>= 1;
    }
    MPI_total_timer += MPI_Wtime();
    double MPI_total_timer_all = 0.0;
    MPI_Reduce(&MPI_total_timer, &MPI_total_timer_all, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    if (rank == 0)
    {
        if (PRINT_EXPERIMENTS)
        {
            printf("%f\n",
                   MPI_total_timer_all / comm_size * 1000000);
        }
    }
fn_exit:
    return mpi_errno;
fn_fail:
    goto fn_exit;
}

int MPI_Bcast_SZx_FXR_RI_record(void *buffer,
                                float compressionRatio,
                                float tolerance,
                                MPI_Aint count,
                                MPI_Datatype datatype,
                                int root, MPI_Comm comm)
{
    double MPI_total_timer = 0.0;
    MPI_total_timer -= MPI_Wtime();
    int rank, comm_size, src, dst;
    int relative_rank, mask;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;
    MPI_Aint nbytes = 0;
    MPI_Status status;
    int tmp_len = 0;
    double absErrBound = compressionRatio;
    if (cpr_mode == FXR)
    {
        tmp_len = (count - 1) / (compressionRatio * (1 - tolerance)) + 1;
    }
    if (cpr_mode == ABS)
    {
        tmp_len = count;
    }
    void *tmpbuf;
    MPI_Request reqs[2]; /* one send and one recv per transfer */
    MPI_Status stas[2];
    int is_contig;
    MPI_Aint type_size;
    void *tmp_buf = NULL;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    // Compression variables
    size_t outSize;
    size_t byteLength;

    type_size = sizeof(data_type);

    nbytes = type_size * count;
    if (nbytes == 0)
        goto fn_exit; /* nothing to do */

    relative_rank = (rank >= root) ? rank - root : rank - root + comm_size;

    /* Use short message algorithm, namely, binomial tree */

    /* Algorithm:
     * This uses a fairly basic recursive subdivision algorithm.
     * The root sends to the process comm_size/2 away; the receiver becomes
     * a root for a subtree and applies the same process.
     *
     * So that the new root can easily identify the size of its
     * subtree, the (subtree) roots are all powers of two (relative
     * to the root) If m = the first power of 2 such that 2^m >= the
     * size of the communicator, then the subtree at root at 2^(m-k)
     * has size 2^k (with special handling for subtrees that aren't
     * a power of two in size).
     *
     * Do subdivision.  There are two phases:
     * 1. Wait for arrival of data.  Because of the power of two nature
     * of the subtree roots, the source of this message is always the
     * process whose relative rank has the least significant 1 bit CLEARED.
     * That is, process 4 (100) receives from process 0, process 7 (111)
     * from process 6 (110), etc.
     * 2. Forward to my subtree
     *
     * Note that the process that is the tree root is handled automatically
     * by this code, since it has no bits set.  */

    mask = 0x1;
    while (mask < comm_size)
    {
        if (relative_rank & mask)
        {
            src = rank - mask;
            if (src < 0)
                src += comm_size;

            mpi_errno = MPI_Recv(buffer, count, datatype, src,
                                 MPIR_BCAST_TAG, comm, &status);
            MPI_Get_count(&status, MPI_FLOAT, &byteLength);
            unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, buffer, byteLength, 0, 0, 0, 0, count);
            memcpy(buffer, data, count * type_size);
            free(data);
            if (mpi_errno)
            {
                exit(-1);
            }
            break;
        }
        mask <<= 1;
    }

    /* This process is responsible for all processes that have bits
     * set from the LSB up to (but not including) mask.  Because of
     * the "not including", we start by shifting mask back down one.
     *
     * We can easily change to a different algorithm at any power of two
     * by changing the test (mask > 1) to (mask > block_size)
     *
     * One such version would use non-blocking operations for the last 2-4
     * steps (this also bounds the number of MPI_Requests that would
     * be needed).  */

    mask >>= 1;
    while (mask > 0)
    {
        if (relative_rank + mask < comm_size)
        {
            dst = rank + mask;
            if (dst >= comm_size)
                dst -= comm_size;
            unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, buffer,
                                                         &outSize, cpr_mode, absErrBound, relBoundRatio, compressionRatio, tolerance, 0, 0, 0, 0, count);
            mpi_errno = MPI_Send(bytes, outSize, MPI_BYTE, dst,
                                 MPIR_BCAST_TAG, comm);
            free(bytes);
            if (mpi_errno)
            {
                exit(-1);
            }
        }
        mask >>= 1;
    }
    MPI_total_timer += MPI_Wtime();
    double MPI_total_timer_all = 0.0;
    MPI_Reduce(&MPI_total_timer, &MPI_total_timer_all, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    if (rank == 0)
    {
        if (PRINT_EXPERIMENTS)
        {
            printf("%f\n",
                   MPI_total_timer_all / comm_size * 1000000);
        }
    }
fn_exit:
    return mpi_errno;
fn_fail:
    goto fn_exit;
}

int MPI_Bcast_SZx_FXR_RI2_oa_record(void *buffer,
                                float compressionRatio,
                                float tolerance,
                                MPI_Aint count,
                                MPI_Datatype datatype,
                                int root, MPI_Comm comm)
{
    double MPI_total_timer = 0.0;
    MPI_total_timer -= MPI_Wtime();
    int rank, comm_size, src, dst;
    int relative_rank, mask;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;
    MPI_Aint nbytes = 0;
    MPI_Status status;
    int tmp_len = 0;
    double absErrBound = compressionRatio;
    if (cpr_mode == FXR)
    {
        tmp_len = (count - 1) / (compressionRatio * (1 - tolerance)) + 1;
    }
    if (cpr_mode == ABS)
    {
        tmp_len = count;
    }

    void *tmpbuf;
    MPI_Request reqs[2]; /* one send and one recv per transfer */
    MPI_Status stas[2];
    int is_contig;
    MPI_Aint type_size;
    void *tmp_buf = NULL;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    type_size = sizeof(data_type);

    // Compression variables
    size_t outSize;
    size_t byteLength;
    //  The compressed size buffer
    int *compressed_sizes = (int *)malloc(comm_size * sizeof(int));
    // For comparess buffer
    unsigned char *outputBytes = (unsigned char *)malloc(tmp_len * type_size);
    if (rank == root)
    {
        SZ_fast_compress_args2(fast_mode, SZ_FLOAT, buffer,
                               &outSize, outputBytes,
                               cpr_mode, absErrBound,
                               relBoundRatio, compressionRatio,
                               tolerance, 0, 0, 0, 0, count);
    }
    MPI_Bcast(&outSize, 1, MPI_INT, root,
              comm);

    nbytes = type_size * count;
    if (nbytes == 0)
        goto fn_exit; /* nothing to do */

    relative_rank = (rank >= root) ? rank - root : rank - root + comm_size;

    /* Use short message algorithm, namely, binomial tree */

    /* Algorithm:
     * This uses a fairly basic recursive subdivision algorithm.
     * The root sends to the process comm_size/2 away; the receiver becomes
     * a root for a subtree and applies the same process.
     *
     * So that the new root can easily identify the size of its
     * subtree, the (subtree) roots are all powers of two (relative
     * to the root) If m = the first power of 2 such that 2^m >= the
     * size of the communicator, then the subtree at root at 2^(m-k)
     * has size 2^k (with special handling for subtrees that aren't
     * a power of two in size).
     *
     * Do subdivision.  There are two phases:
     * 1. Wait for arrival of data.  Because of the power of two nature
     * of the subtree roots, the source of this message is always the
     * process whose relative rank has the least significant 1 bit CLEARED.
     * That is, process 4 (100) receives from process 0, process 7 (111)
     * from process 6 (110), etc.
     * 2. Forward to my subtree
     *
     * Note that the process that is the tree root is handled automatically
     * by this code, since it has no bits set.  */

    mask = 0x1;
    while (mask < comm_size)
    {
        if (relative_rank & mask)
        {
            src = rank - mask;
            if (src < 0)
                src += comm_size;

            mpi_errno = MPI_Recv(outputBytes, outSize, MPI_BYTE, src,
                                 MPIR_BCAST_TAG, comm, &status);
            if (mpi_errno)
            {
                exit(-1);
            }
            break;
        }
        mask <<= 1;
    }

    /* This process is responsible for all processes that have bits
     * set from the LSB up to (but not including) mask.  Because of
     * the "not including", we start by shifting mask back down one.
     *
     * We can easily change to a different algorithm at any power of two
     * by changing the test (mask > 1) to (mask > block_size)
     *
     * One such version would use non-blocking operations for the last 2-4
     * steps (this also bounds the number of MPI_Requests that would
     * be needed).  */

    mask >>= 1;
    while (mask > 0)
    {
        if (relative_rank + mask < comm_size)
        {
            dst = rank + mask;
            if (dst >= comm_size)
                dst -= comm_size;
            mpi_errno = MPI_Send(outputBytes, outSize, MPI_BYTE, dst,
                                 MPIR_BCAST_TAG, comm);
            if (mpi_errno)
            {
                exit(-1);
            }
        }
        mask >>= 1;
    }

    if(rank != root)
    {
        SZ_fast_decompress_split(fast_mode, SZ_FLOAT, buffer, 
        (char *)outputBytes, 0, 0, 0, 0, count);
    }
    free(outputBytes);
    MPI_total_timer += MPI_Wtime();
    double MPI_total_timer_all = 0.0;
    MPI_Reduce(&MPI_total_timer, &MPI_total_timer_all, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    if (rank == 0)
    {
        if (PRINT_EXPERIMENTS)
        {
            printf("%f\n",
                   MPI_total_timer_all / comm_size * 1000000);
        }
    }
fn_exit:
    return mpi_errno;
fn_fail:
    goto fn_exit;
}

int MPI_Bcast_SZx_FXR_RI2_oa_mt_record(void *buffer,
                                float compressionRatio,
                                float tolerance,
                                MPI_Aint count,
                                MPI_Datatype datatype,
                                int root, MPI_Comm comm)
{
    double MPI_total_timer = 0.0;
    MPI_total_timer -= MPI_Wtime();
    int rank, comm_size, src, dst;
    int relative_rank, mask;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;
    MPI_Aint nbytes = 0;
    MPI_Status status;
    int tmp_len = 0;
    double absErrBound = compressionRatio;
    if (cpr_mode == FXR)
    {
        tmp_len = (count - 1) / (compressionRatio * (1 - tolerance)) + 1;
    }
    if (cpr_mode == ABS)
    {
        tmp_len = count;
    }

    void *tmpbuf;
    MPI_Request reqs[2]; /* one send and one recv per transfer */
    MPI_Status stas[2];
    int is_contig;
    MPI_Aint type_size;
    void *tmp_buf = NULL;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    type_size = sizeof(data_type);

    // Compression variables
    size_t outSize;
    size_t byteLength;
    //  The compressed size buffer
    int *compressed_sizes = (int *)malloc(comm_size * sizeof(int));
    // For comparess buffer
    unsigned char *outputBytes = (unsigned char *)malloc(tmp_len * type_size);
    if (rank == root)
    {
        SZ_fast_compress_args2(4, SZ_FLOAT, buffer,
                               &outSize, outputBytes,
                               cpr_mode, absErrBound,
                               relBoundRatio, compressionRatio,
                               tolerance, 0, 0, 0, 0, count);
    }
    MPI_Bcast(&outSize, 1, MPI_INT, root,
              comm);

    nbytes = type_size * count;
    if (nbytes == 0)
        goto fn_exit; /* nothing to do */

    relative_rank = (rank >= root) ? rank - root : rank - root + comm_size;

    /* Use short message algorithm, namely, binomial tree */

    /* Algorithm:
     * This uses a fairly basic recursive subdivision algorithm.
     * The root sends to the process comm_size/2 away; the receiver becomes
     * a root for a subtree and applies the same process.
     *
     * So that the new root can easily identify the size of its
     * subtree, the (subtree) roots are all powers of two (relative
     * to the root) If m = the first power of 2 such that 2^m >= the
     * size of the communicator, then the subtree at root at 2^(m-k)
     * has size 2^k (with special handling for subtrees that aren't
     * a power of two in size).
     *
     * Do subdivision.  There are two phases:
     * 1. Wait for arrival of data.  Because of the power of two nature
     * of the subtree roots, the source of this message is always the
     * process whose relative rank has the least significant 1 bit CLEARED.
     * That is, process 4 (100) receives from process 0, process 7 (111)
     * from process 6 (110), etc.
     * 2. Forward to my subtree
     *
     * Note that the process that is the tree root is handled automatically
     * by this code, since it has no bits set.  */

    mask = 0x1;
    while (mask < comm_size)
    {
        if (relative_rank & mask)
        {
            src = rank - mask;
            if (src < 0)
                src += comm_size;

            mpi_errno = MPI_Recv(outputBytes, outSize, MPI_BYTE, src,
                                 MPIR_BCAST_TAG, comm, &status);
            if (mpi_errno)
            {
                exit(-1);
            }
            break;
        }
        mask <<= 1;
    }

    /* This process is responsible for all processes that have bits
     * set from the LSB up to (but not including) mask.  Because of
     * the "not including", we start by shifting mask back down one.
     *
     * We can easily change to a different algorithm at any power of two
     * by changing the test (mask > 1) to (mask > block_size)
     *
     * One such version would use non-blocking operations for the last 2-4
     * steps (this also bounds the number of MPI_Requests that would
     * be needed).  */

    mask >>= 1;
    while (mask > 0)
    {
        if (relative_rank + mask < comm_size)
        {
            dst = rank + mask;
            if (dst >= comm_size)
                dst -= comm_size;
            mpi_errno = MPI_Send(outputBytes, outSize, MPI_BYTE, dst,
                                 MPIR_BCAST_TAG, comm);
            if (mpi_errno)
            {
                exit(-1);
            }
        }
        mask >>= 1;
    }

    if(rank != root)
    {
        SZ_fast_decompress_split(4, SZ_FLOAT, buffer, 
        (char *)outputBytes, 0, 0, 0, 0, count);
    }
    free(outputBytes);
    MPI_total_timer += MPI_Wtime();
    double MPI_total_timer_all = 0.0;
    MPI_Reduce(&MPI_total_timer, &MPI_total_timer_all, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    if (rank == 0)
    {
        if (PRINT_EXPERIMENTS)
        {
            printf("%f\n",
                   MPI_total_timer_all / comm_size * 1000000);
        }
    }
fn_exit:
    return mpi_errno;
fn_fail:
    goto fn_exit;
}