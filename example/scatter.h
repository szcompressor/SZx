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

/* This is the machine-independent implementation of scatter. The algorithm is:

   Algorithm: Binomial

   We use a binomial tree algorithm for both short and
   long messages. At nodes other than leaf nodes we need to allocate
   a temporary buffer to store the incoming message. If the root is
   not rank 0, we reorder the sendbuf in order of relative ranks by
   copying it into a temporary buffer, so that all the sends from the
   root are contiguous and in the right order.

   Cost = lgp.alpha + n.((p-1)/p).beta
   where n is the total size of the data to be scattered from the root.

   Possible improvements:

   End Algorithm: MPI_Scatter
*/

/* not declared static because a machine-specific function may call this one in some cases */
int MPIR_Scatter_intra_binomial(const void *sendbuf, MPI_Aint sendcount, MPI_Datatype sendtype,
                                void *recvbuf, MPI_Aint recvcount, MPI_Datatype recvtype, int root,
                                MPI_Comm comm)
{
    double MPI_total_timer = 0.0;
    MPI_total_timer -= MPI_Wtime();
    MPI_Status status;
    MPI_Aint extent = 0;
    int rank, comm_size;
    int relative_rank;
    MPI_Aint curr_cnt, send_subtree_cnt;
    int mask, src, dst;
    MPI_Aint tmp_buf_size = 0;
    void *tmp_buf = NULL;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    if (rank == root)
        extent = sizeof(sendtype);

    relative_rank = (rank >= root) ? rank - root : rank - root + comm_size;

    MPI_Aint nbytes;
    if (rank == root)
    {
        /* We separate the two cases (root and non-root) because
         * in the event of recvbuf=MPI_IN_PLACE on the root,
         * recvcount and recvtype are not valid */
        MPI_Aint sendtype_size;
        sendtype_size = sizeof(sendtype);
        nbytes = sendtype_size * sendcount;
    }
    else
    {
        MPI_Aint recvtype_size;
        recvtype_size = sizeof(recvtype);
        nbytes = recvtype_size * recvcount;
    }

    curr_cnt = 0;

    /* all even nodes other than root need a temporary buffer to
     * receive data of max size (nbytes*comm_size)/2 */
    if (relative_rank && !(relative_rank % 2))
    {
        tmp_buf_size = (nbytes * comm_size) / 2;
        tmp_buf = (void *)malloc(tmp_buf_size);
    }

    /* if the root is not rank 0, we reorder the sendbuf in order of
     * relative ranks and copy it into a temporary buffer, so that
     * all the sends from the root are contiguous and in the right
     * order. */
    if (rank == root)
    {
        if (root != 0)
        {
            tmp_buf_size = nbytes * comm_size;
            tmp_buf = (void *)malloc(tmp_buf_size);

            if (recvbuf != MPI_IN_PLACE)
                memcpy(tmp_buf, ((char *)sendbuf + extent * sendcount * rank),
                       sendcount * (comm_size - rank) * sizeof(sendtype));
            else
                memcpy((char *)tmp_buf + nbytes, ((char *)sendbuf + extent * sendcount * (rank + 1)),
                       sendcount * (comm_size - rank - 1) * sizeof(sendtype));

            memcpy(((char *)tmp_buf + nbytes * (comm_size - rank)), sendbuf,
                   sendcount * rank * sizeof(sendtype));

            curr_cnt = nbytes * comm_size;
        }
        else
            curr_cnt = sendcount * comm_size;
    }

    /* root has all the data; others have zero so far */

    mask = 0x1;
    while (mask < comm_size)
    {
        if (relative_rank & mask)
        {
            src = rank - mask;
            if (src < 0)
                src += comm_size;

            /* The leaf nodes receive directly into recvbuf because
             * they don't have to forward data to anyone. Others
             * receive data into a temporary buffer. */
            if (relative_rank % 2)
            {
                mpi_errno = MPI_Recv(recvbuf, recvcount, recvtype,
                                     src, MPIR_SCATTER_TAG, comm, &status);
                if (mpi_errno)
                {
                    exit(-1);
                }
            }
            else
            {
                mpi_errno = MPI_Recv(tmp_buf, tmp_buf_size, MPI_BYTE, src,
                                     MPIR_SCATTER_TAG, comm, &status);
                if (mpi_errno)
                {
                    exit(-1);
                }
                else
                    /* the recv size is larger than what may be sent in
                     * some cases. query amount of data actually received */
                    MPI_Get_count(&status, MPI_BYTE, &curr_cnt);
            }
            break;
        }
        mask <<= 1;
    }

    /* This process is responsible for all processes that have bits
     * set from the LSB up to (but not including) mask.  Because of
     * the "not including", we start by shifting mask back down
     * one. */

    mask >>= 1;
    while (mask > 0)
    {
        if (relative_rank + mask < comm_size)
        {
            dst = rank + mask;
            if (dst >= comm_size)
                dst -= comm_size;

            if ((rank == root) && (root == 0))
            {
                send_subtree_cnt = curr_cnt - sendcount * mask;
                /* mask is also the size of this process's subtree */
                mpi_errno = MPI_Send(((char *)sendbuf +
                                      extent * sendcount * mask),
                                     send_subtree_cnt,
                                     sendtype, dst,
                                     MPIR_SCATTER_TAG, comm);
            }
            else
            {
                /* non-zero root and others */
                send_subtree_cnt = curr_cnt - nbytes * mask;
                /* mask is also the size of this process's subtree */
                mpi_errno = MPI_Send(((char *)tmp_buf + nbytes * mask),
                                     send_subtree_cnt,
                                     MPI_BYTE, dst, MPIR_SCATTER_TAG, comm);
            }
            if (mpi_errno)
            {
                exit(-1);
            }
            curr_cnt -= send_subtree_cnt;
        }
        mask >>= 1;
    }

    if ((rank == root) && (root == 0) && (recvbuf != MPI_IN_PLACE))
    {
        /* for root=0, put root's data in recvbuf if not MPI_IN_PLACE */
        memcpy(recvbuf, sendbuf, sendcount * sizeof(sendtype));
    }
    else if (!(relative_rank % 2) && (recvbuf != MPI_IN_PLACE))
    {
        /* for non-zero root and non-leaf nodes, copy from tmp_buf
         * into recvbuf */
        memcpy(recvbuf, tmp_buf, nbytes);
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
    free(tmp_buf);
    if (mpi_errno_ret)
        mpi_errno = mpi_errno_ret;
    return mpi_errno;
fn_fail:
    goto fn_exit;
}

int MPIR_Scatter_SZx_FXR_RI_record(const void *sendbuf,
                                   float compressionRatio,
                                   float tolerance,
                                   MPI_Aint sendcount, MPI_Datatype sendtype,
                                   void *recvbuf, MPI_Aint recvcount, MPI_Datatype recvtype, int root,
                                   MPI_Comm comm)
{
    double MPI_total_timer = 0.0;
    MPI_total_timer -= MPI_Wtime();
    MPI_Status status;
    MPI_Aint extent = 0;
    int rank, comm_size;
    int relative_rank;
    MPI_Aint curr_cnt, send_subtree_cnt;
    int mask, src, dst;
    MPI_Aint tmp_buf_size = 0;
    void *tmp_buf = NULL;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    if (rank == root)
        extent = sizeof(sendtype);

    int tmp_len = 0;
    double absErrBound = compressionRatio;
    if (cpr_mode == FXR)
    {
        tmp_len = (sendcount - 1) / (compressionRatio * (1 - tolerance)) + 1;
    }
    if (cpr_mode == ABS)
    {
        tmp_len = sendcount;
    }
    void *tmpbuf;
    MPI_Request reqs[2]; /* one send and one recv per transfer */
    MPI_Status stas[2];
    // Compression variables
    size_t outSize;
    size_t byteLength;

    relative_rank = (rank >= root) ? rank - root : rank - root + comm_size;

    MPI_Aint nbytes;
    if (rank == root)
    {
        /* We separate the two cases (root and non-root) because
         * in the event of recvbuf=MPI_IN_PLACE on the root,
         * recvcount and recvtype are not valid */
        MPI_Aint sendtype_size;
        sendtype_size = sizeof(sendtype);
        nbytes = sendtype_size * sendcount;
    }
    else
    {
        MPI_Aint recvtype_size;
        recvtype_size = sizeof(recvtype);
        nbytes = recvtype_size * recvcount;
    }

    curr_cnt = 0;

    /* all even nodes other than root need a temporary buffer to
     * receive data of max size (nbytes*comm_size)/2 */
    if (relative_rank && !(relative_rank % 2))
    {
        tmp_buf_size = (nbytes * comm_size) / 2;
        tmp_buf = (void *)malloc(tmp_buf_size);
    }

    /* if the root is not rank 0, we reorder the sendbuf in order of
     * relative ranks and copy it into a temporary buffer, so that
     * all the sends from the root are contiguous and in the right
     * order. */
    if (rank == root)
    {
        if (root != 0)
        {
            tmp_buf_size = nbytes * comm_size;
            tmp_buf = (void *)malloc(tmp_buf_size);

            if (recvbuf != MPI_IN_PLACE)
                memcpy(tmp_buf, ((char *)sendbuf + extent * sendcount * rank),
                       sendcount * (comm_size - rank) * sizeof(sendtype));
            else
                memcpy((char *)tmp_buf + nbytes, ((char *)sendbuf + extent * sendcount * (rank + 1)),
                       sendcount * (comm_size - rank - 1) * sizeof(sendtype));

            memcpy(((char *)tmp_buf + nbytes * (comm_size - rank)), sendbuf,
                   sendcount * rank * sizeof(sendtype));

            curr_cnt = nbytes * comm_size;
        }
        else
            curr_cnt = sendcount * comm_size;
    }

    /* root has all the data; others have zero so far */

    mask = 0x1;
    while (mask < comm_size)
    {
        if (relative_rank & mask)
        {
            src = rank - mask;
            if (src < 0)
                src += comm_size;

            /* The leaf nodes receive directly into recvbuf because
             * they don't have to forward data to anyone. Others
             * receive data into a temporary buffer. */
            if (relative_rank % 2)
            {
                int original_size;
                MPI_Recv(&original_size, 1, MPI_INT, src,
                         MPIR_SCATTER_TAG, comm, &status);
                mpi_errno = MPI_Recv(recvbuf, recvcount, recvtype,
                                     src, MPIR_SCATTER_TAG, comm, &status);
                MPI_Get_count(&status, MPI_FLOAT, &byteLength);
                unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, recvbuf, byteLength, 0, 0, 0, 0, recvcount);
                memcpy(recvbuf, data, recvcount * sizeof(recvtype));
                free(data);
                if (mpi_errno)
                {
                    exit(-1);
                }
            }
            else
            {
                int original_size;
                MPI_Recv(&original_size, 1, MPI_INT, src,
                         MPIR_SCATTER_TAG, comm, &status);
                // printf("rank %d finished recv data size\n", rank);
                mpi_errno = MPI_Recv(tmp_buf, tmp_buf_size, MPI_BYTE, src,
                                     MPIR_SCATTER_TAG, comm, &status);
                MPI_Get_count(&status, MPI_FLOAT, &byteLength);
                unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, tmp_buf, byteLength, 0, 0, 0, 0, original_size / sizeof(recvtype));
                memcpy(tmp_buf, data, original_size);
                free(data);
                if (mpi_errno)
                {
                    exit(-1);
                }
                else
                    /* the recv size is larger than what may be sent in
                     * some cases. query amount of data actually received */
                    // MPI_Get_count(&status, MPI_BYTE, &curr_cnt);
                    curr_cnt = original_size;
            }
            break;
        }
        mask <<= 1;
    }

    /* This process is responsible for all processes that have bits
     * set from the LSB up to (but not including) mask.  Because of
     * the "not including", we start by shifting mask back down
     * one. */

    mask >>= 1;
    while (mask > 0)
    {
        if (relative_rank + mask < comm_size)
        {
            dst = rank + mask;
            if (dst >= comm_size)
                dst -= comm_size;

            if ((rank == root) && (root == 0))
            {
                send_subtree_cnt = curr_cnt - sendcount * mask;
                /* mask is also the size of this process's subtree */
                unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, ((char *)sendbuf + extent * sendcount * mask),
                                                             &outSize, cpr_mode, absErrBound,
                                                             relBoundRatio, compressionRatio, tolerance,
                                                             0, 0, 0, 0, send_subtree_cnt);
                // printf("root finished compression with mask %d the compressed data size is %d, the dst is %d\n", mask, outSize, dst);
                int original_size = send_subtree_cnt * sizeof(sendtype);
                MPI_Send(&original_size,
                         1,
                         MPI_INT, dst,
                         MPIR_SCATTER_TAG, comm);

                // printf("root finished send size with mask %d\n", mask);
                mpi_errno = MPI_Send(bytes,
                                     outSize,
                                     MPI_BYTE, dst,
                                     MPIR_SCATTER_TAG, comm);
                free(bytes);

                // printf("root sent to dst %d with mask %d\n", dst, mask);
            }
            else
            {
                /* non-zero root and others */
                send_subtree_cnt = curr_cnt - nbytes * mask;
                /* mask is also the size of this process's subtree */
                unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, ((char *)tmp_buf + nbytes * mask),
                                                             &outSize, cpr_mode, absErrBound,
                                                             relBoundRatio, compressionRatio, tolerance,
                                                             0, 0, 0, 0, send_subtree_cnt / sizeof(sendtype));

                MPI_Send(&send_subtree_cnt,
                         1,
                         MPI_INT, dst,
                         MPIR_SCATTER_TAG, comm);

                mpi_errno = MPI_Send(bytes,
                                     outSize,
                                     MPI_BYTE, dst, MPIR_SCATTER_TAG, comm);
                free(bytes);
            }
            if (mpi_errno)
            {
                exit(-1);
            }
            curr_cnt -= send_subtree_cnt;
        }
        mask >>= 1;
    }

    if ((rank == root) && (root == 0) && (recvbuf != MPI_IN_PLACE))
    {
        /* for root=0, put root's data in recvbuf if not MPI_IN_PLACE */
        memcpy(recvbuf, sendbuf, sendcount * sizeof(sendtype));
    }
    else if (!(relative_rank % 2) && (recvbuf != MPI_IN_PLACE))
    {
        /* for non-zero root and non-leaf nodes, copy from tmp_buf
         * into recvbuf */
        memcpy(recvbuf, tmp_buf, nbytes);
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
    free(tmp_buf);
    if (mpi_errno_ret)
        mpi_errno = mpi_errno_ret;
    return mpi_errno;
fn_fail:
    goto fn_exit;
}

int MPIR_Scatter_SZx_FXR_RI2_oa_record(const void *sendbuf,
                                       float compressionRatio,
                                       float tolerance,
                                       MPI_Aint sendcount, MPI_Datatype sendtype,
                                       void *recvbuf, MPI_Aint recvcount, MPI_Datatype recvtype, int root,
                                       MPI_Comm comm)
{
    double MPI_total_timer = 0.0;
    MPI_total_timer -= MPI_Wtime();
    MPI_Status status;
    MPI_Aint extent = 0;
    int rank, comm_size;
    int relative_rank;
    MPI_Aint curr_cnt, send_subtree_cnt;
    int mask, src, dst;
    MPI_Aint tmp_buf_size = 0;
    void *tmp_buf = NULL;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    if (rank == root)
        extent = sizeof(sendtype);

    int tmp_len = 0;
    double absErrBound = compressionRatio;
    if (cpr_mode == FXR)
    {
        tmp_len = (sendcount * comm_size - 1) / (compressionRatio * (1 - tolerance)) + 1;
    }
    if (cpr_mode == ABS)
    {
        tmp_len = sendcount * comm_size;
    }
    void *tmpbuf;
    MPI_Request reqs[2]; /* one send and one recv per transfer */
    MPI_Status stas[2];
    // Compression variables
    size_t outSize;
    size_t byteLength;
    //  The compressed size buffer
    int *compressed_sizes = (int *)malloc(comm_size * sizeof(int));
    int *compressed_offsets = (int *)malloc(comm_size * sizeof(int));
    int *compressed_lefts = (int *)malloc(comm_size * sizeof(int));

    // For comparess buffer
    unsigned char *outputBytes = (unsigned char *)malloc(tmp_len * sizeof(sendtype));

    relative_rank = (rank >= root) ? rank - root : rank - root + comm_size;

    MPI_Aint nbytes;
    if (rank == root)
    {
        /* We separate the two cases (root and non-root) because
         * in the event of recvbuf=MPI_IN_PLACE on the root,
         * recvcount and recvtype are not valid */
        MPI_Aint sendtype_size;
        sendtype_size = sizeof(sendtype);
        nbytes = sendtype_size * sendcount;
    }
    else
    {
        MPI_Aint recvtype_size;
        recvtype_size = sizeof(recvtype);
        nbytes = recvtype_size * recvcount;
    }

    curr_cnt = 0;

    /* all even nodes other than root need a temporary buffer to
     * receive data of max size (nbytes*comm_size)/2 */
    if (relative_rank && !(relative_rank % 2))
    {
        tmp_buf_size = (nbytes * comm_size) / 2;
        tmp_buf = (void *)malloc(tmp_buf_size);
    }

    /* if the root is not rank 0, we reorder the sendbuf in order of
     * relative ranks and copy it into a temporary buffer, so that
     * all the sends from the root are contiguous and in the right
     * order. */
    if (rank == root)
    {
        if (root != 0)
        {
            tmp_buf_size = nbytes * comm_size;
            tmp_buf = (void *)malloc(tmp_buf_size);

            if (recvbuf != MPI_IN_PLACE)
                memcpy(tmp_buf, ((char *)sendbuf + extent * sendcount * rank),
                       sendcount * (comm_size - rank) * sizeof(sendtype));
            else
                memcpy((char *)tmp_buf + nbytes, ((char *)sendbuf + extent * sendcount * (rank + 1)),
                       sendcount * (comm_size - rank - 1) * sizeof(sendtype));

            memcpy(((char *)tmp_buf + nbytes * (comm_size - rank)), sendbuf,
                   sendcount * rank * sizeof(sendtype));

            curr_cnt = nbytes * comm_size;
        }
        else
            curr_cnt = sendcount * comm_size;
    }

    /* root has all the data; others have zero so far */

    if (rank == root && root == 0)
    {
        int offset = 0;
        int i = 0, j = 0;
        for (i = 0; i < curr_cnt; i += sendcount)
        {
            SZ_fast_compress_args2(fast_mode, SZ_FLOAT, (char *)sendbuf + i * sizeof(recvtype),
                                   &outSize, (char *)outputBytes + offset,
                                   cpr_mode, absErrBound,
                                   relBoundRatio, compressionRatio,
                                   tolerance, 0, 0, 0, 0, sendcount);
            compressed_sizes[j] = outSize;
            compressed_offsets[j] = offset;
            offset += outSize;
            j++;
        }
        curr_cnt = offset;
        for (j = 0; j < comm_size; j++)
        {
            compressed_lefts[j] = curr_cnt - compressed_offsets[j];
        }
    }
    else if (rank == root && root != 0)
    {
        int offset = 0;
        int i = 0, j = 0;
        for (i = 0; i < curr_cnt; i += sendcount)
        {
            SZ_fast_compress_args2(fast_mode, SZ_FLOAT, (char *)tmp_buf + i * sizeof(recvtype),
                                   &outSize, (char *)outputBytes + offset,
                                   cpr_mode, absErrBound,
                                   relBoundRatio, compressionRatio,
                                   tolerance, 0, 0, 0, 0, sendcount);
            compressed_sizes[j] = outSize;
            compressed_offsets[j] = offset;
            offset += outSize;
            j++;
        }
        curr_cnt = offset;
        for (j = 0; j < comm_size; j++)
        {
            compressed_lefts[j] = curr_cnt - compressed_offsets[j];
        }
    }

    MPI_Bcast(compressed_sizes, comm_size, MPI_INT, root,
              comm);
    MPI_Bcast(compressed_offsets, comm_size, MPI_INT, root,
              comm);
    MPI_Bcast(compressed_lefts, comm_size, MPI_INT, root,
              comm);

    mask = 0x1;
    while (mask < comm_size)
    {
        if (relative_rank & mask)
        {
            src = rank - mask;
            if (src < 0)
                src += comm_size;

            /* The leaf nodes receive directly into recvbuf because
             * they don't have to forward data to anyone. Others
             * receive data into a temporary buffer. */
            if (relative_rank % 2)
            {
                mpi_errno = MPI_Recv(outputBytes, recvcount, recvtype,
                                     src, MPIR_SCATTER_TAG, comm, &status);
                if (mpi_errno)
                {
                    exit(-1);
                }
            }
            else
            {
                mpi_errno = MPI_Recv(outputBytes, tmp_buf_size, MPI_BYTE, src,
                                     MPIR_SCATTER_TAG, comm, &status);
                if (mpi_errno)
                {
                    exit(-1);
                }
                else
                    /* the recv size is larger than what may be sent in
                     * some cases. query amount of data actually received */
                    MPI_Get_count(&status, MPI_BYTE, &curr_cnt);
            }
            break;
        }
        mask <<= 1;
    }

    /* This process is responsible for all processes that have bits
     * set from the LSB up to (but not including) mask.  Because of
     * the "not including", we start by shifting mask back down
     * one. */

    mask >>= 1;
    while (mask > 0)
    {
        if (relative_rank + mask < comm_size)
        {
            dst = rank + mask;
            if (dst >= comm_size)
                dst -= comm_size;

            if ((rank == root) && (root == 0))
            {
                send_subtree_cnt = curr_cnt - compressed_offsets[mask];
                /* mask is also the size of this process's subtree */
                mpi_errno = MPI_Send(outputBytes + compressed_offsets[dst],
                                     send_subtree_cnt,
                                     MPI_BYTE, dst,
                                     MPIR_SCATTER_TAG, comm);
            }
            else
            {
                /* non-zero root and others */
                send_subtree_cnt = compressed_offsets[dst + mask - 1] - compressed_offsets[dst] + compressed_sizes[dst + mask - 1];
                /* mask is also the size of this process's subtree */
                mpi_errno = MPI_Send(((char *)outputBytes + compressed_offsets[dst] - compressed_offsets[rank]),
                                     send_subtree_cnt,
                                     MPI_BYTE, dst, MPIR_SCATTER_TAG, comm);
            }
            if (mpi_errno)
            {
                exit(-1);
            }
            curr_cnt -= send_subtree_cnt;
        }
        mask >>= 1;
    }

    // printf("compression finished at rank %d\n", rank);

    if ((relative_rank % 2))
    {
        SZ_fast_decompress_split(fast_mode, SZ_FLOAT, recvbuf,
                                 (char *)outputBytes, 0, 0, 0, 0, recvcount);
    }
    if ((rank == root) && (root == 0) && (recvbuf != MPI_IN_PLACE))
    {
        /* for root=0, put root's data in recvbuf if not MPI_IN_PLACE */
        memcpy(recvbuf, sendbuf, sendcount * sizeof(sendtype));
    }
    else if (!(relative_rank % 2) && (recvbuf != MPI_IN_PLACE))
    {
        /* for non-zero root and non-leaf nodes, copy from tmp_buf
         * into recvbuf */
        SZ_fast_decompress_split(fast_mode, SZ_FLOAT, recvbuf,
                                 (char *)outputBytes, 0, 0, 0, 0, recvcount);
        // memcpy(recvbuf, tmp_buf, nbytes);
    }

    free(outputBytes);
    free(compressed_lefts);
    free(compressed_offsets);
    free(compressed_sizes);
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
    free(tmp_buf);
    if (mpi_errno_ret)
        mpi_errno = mpi_errno_ret;
    return mpi_errno;
fn_fail:
    goto fn_exit;
}

int MPIR_Scatter_SZx_FXR_RI2_oa_mt_record(const void *sendbuf,
                                          float compressionRatio,
                                          float tolerance,
                                          MPI_Aint sendcount, MPI_Datatype sendtype,
                                          void *recvbuf, MPI_Aint recvcount, MPI_Datatype recvtype, int root,
                                          MPI_Comm comm)
{
    double MPI_total_timer = 0.0;
    MPI_total_timer -= MPI_Wtime();
    MPI_Status status;
    MPI_Aint extent = 0;
    int rank, comm_size;
    int relative_rank;
    MPI_Aint curr_cnt, send_subtree_cnt;
    int mask, src, dst;
    MPI_Aint tmp_buf_size = 0;
    void *tmp_buf = NULL;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    if (rank == root)
        extent = sizeof(sendtype);

    int tmp_len = 0;
    double absErrBound = compressionRatio;
    if (cpr_mode == FXR)
    {
        tmp_len = (sendcount * comm_size - 1) / (compressionRatio * (1 - tolerance)) + 1;
    }
    if (cpr_mode == ABS)
    {
        tmp_len = sendcount * comm_size;
    }
    void *tmpbuf;
    MPI_Request reqs[2]; /* one send and one recv per transfer */
    MPI_Status stas[2];
    // Compression variables
    size_t outSize;
    size_t byteLength;
    //  The compressed size buffer
    int *compressed_sizes = (int *)malloc(comm_size * sizeof(int));
    int *compressed_offsets = (int *)malloc(comm_size * sizeof(int));
    int *compressed_lefts = (int *)malloc(comm_size * sizeof(int));

    // For comparess buffer
    unsigned char *outputBytes = (unsigned char *)malloc(tmp_len * sizeof(sendtype));

    relative_rank = (rank >= root) ? rank - root : rank - root + comm_size;

    MPI_Aint nbytes;
    if (rank == root)
    {
        /* We separate the two cases (root and non-root) because
         * in the event of recvbuf=MPI_IN_PLACE on the root,
         * recvcount and recvtype are not valid */
        MPI_Aint sendtype_size;
        sendtype_size = sizeof(sendtype);
        nbytes = sendtype_size * sendcount;
    }
    else
    {
        MPI_Aint recvtype_size;
        recvtype_size = sizeof(recvtype);
        nbytes = recvtype_size * recvcount;
    }

    curr_cnt = 0;

    /* all even nodes other than root need a temporary buffer to
     * receive data of max size (nbytes*comm_size)/2 */
    if (relative_rank && !(relative_rank % 2))
    {
        tmp_buf_size = (nbytes * comm_size) / 2;
        tmp_buf = (void *)malloc(tmp_buf_size);
    }

    /* if the root is not rank 0, we reorder the sendbuf in order of
     * relative ranks and copy it into a temporary buffer, so that
     * all the sends from the root are contiguous and in the right
     * order. */
    if (rank == root)
    {
        if (root != 0)
        {
            tmp_buf_size = nbytes * comm_size;
            tmp_buf = (void *)malloc(tmp_buf_size);

            if (recvbuf != MPI_IN_PLACE)
                memcpy(tmp_buf, ((char *)sendbuf + extent * sendcount * rank),
                       sendcount * (comm_size - rank) * sizeof(sendtype));
            else
                memcpy((char *)tmp_buf + nbytes, ((char *)sendbuf + extent * sendcount * (rank + 1)),
                       sendcount * (comm_size - rank - 1) * sizeof(sendtype));

            memcpy(((char *)tmp_buf + nbytes * (comm_size - rank)), sendbuf,
                   sendcount * rank * sizeof(sendtype));

            curr_cnt = nbytes * comm_size;
        }
        else
            curr_cnt = sendcount * comm_size;
    }

    /* root has all the data; others have zero so far */

    if (rank == root && root == 0)
    {
        int offset = 0;
        int i = 0, j = 0;
        for (i = 0; i < curr_cnt; i += sendcount)
        {
            SZ_fast_compress_args2(4, SZ_FLOAT, (char *)sendbuf + i * sizeof(recvtype),
                                   &outSize, (char *)outputBytes + offset,
                                   cpr_mode, absErrBound,
                                   relBoundRatio, compressionRatio,
                                   tolerance, 0, 0, 0, 0, sendcount);
            compressed_sizes[j] = outSize;
            compressed_offsets[j] = offset;
            offset += outSize;
            j++;
        }
        curr_cnt = offset;
        for (j = 0; j < comm_size; j++)
        {
            compressed_lefts[j] = curr_cnt - compressed_offsets[j];
        }
    }
    else if (rank == root && root != 0)
    {
        int offset = 0;
        int i = 0, j = 0;
        for (i = 0; i < curr_cnt; i += sendcount)
        {
            SZ_fast_compress_args2(4, SZ_FLOAT, (char *)tmp_buf + i * sizeof(recvtype),
                                   &outSize, (char *)outputBytes + offset,
                                   cpr_mode, absErrBound,
                                   relBoundRatio, compressionRatio,
                                   tolerance, 0, 0, 0, 0, sendcount);
            compressed_sizes[j] = outSize;
            compressed_offsets[j] = offset;
            offset += outSize;
            j++;
        }
        curr_cnt = offset;
        for (j = 0; j < comm_size; j++)
        {
            compressed_lefts[j] = curr_cnt - compressed_offsets[j];
        }
    }

    MPI_Bcast(compressed_sizes, comm_size, MPI_INT, root,
              comm);
    MPI_Bcast(compressed_offsets, comm_size, MPI_INT, root,
              comm);
    MPI_Bcast(compressed_lefts, comm_size, MPI_INT, root,
              comm);

    mask = 0x1;
    while (mask < comm_size)
    {
        if (relative_rank & mask)
        {
            src = rank - mask;
            if (src < 0)
                src += comm_size;

            /* The leaf nodes receive directly into recvbuf because
             * they don't have to forward data to anyone. Others
             * receive data into a temporary buffer. */
            if (relative_rank % 2)
            {
                mpi_errno = MPI_Recv(outputBytes, recvcount, recvtype,
                                     src, MPIR_SCATTER_TAG, comm, &status);
                if (mpi_errno)
                {
                    exit(-1);
                }
            }
            else
            {
                mpi_errno = MPI_Recv(outputBytes, tmp_buf_size, MPI_BYTE, src,
                                     MPIR_SCATTER_TAG, comm, &status);
                if (mpi_errno)
                {
                    exit(-1);
                }
                else
                    /* the recv size is larger than what may be sent in
                     * some cases. query amount of data actually received */
                    MPI_Get_count(&status, MPI_BYTE, &curr_cnt);
            }
            break;
        }
        mask <<= 1;
    }

    /* This process is responsible for all processes that have bits
     * set from the LSB up to (but not including) mask.  Because of
     * the "not including", we start by shifting mask back down
     * one. */

    mask >>= 1;
    while (mask > 0)
    {
        if (relative_rank + mask < comm_size)
        {
            dst = rank + mask;
            if (dst >= comm_size)
                dst -= comm_size;

            if ((rank == root) && (root == 0))
            {
                send_subtree_cnt = curr_cnt - compressed_offsets[mask];
                /* mask is also the size of this process's subtree */
                mpi_errno = MPI_Send(outputBytes + compressed_offsets[dst],
                                     send_subtree_cnt,
                                     MPI_BYTE, dst,
                                     MPIR_SCATTER_TAG, comm);
            }
            else
            {
                /* non-zero root and others */
                send_subtree_cnt = compressed_offsets[dst + mask - 1] - compressed_offsets[dst] + compressed_sizes[dst + mask - 1];
                /* mask is also the size of this process's subtree */
                mpi_errno = MPI_Send(((char *)outputBytes + compressed_offsets[dst] - compressed_offsets[rank]),
                                     send_subtree_cnt,
                                     MPI_BYTE, dst, MPIR_SCATTER_TAG, comm);
            }
            if (mpi_errno)
            {
                exit(-1);
            }
            curr_cnt -= send_subtree_cnt;
        }
        mask >>= 1;
    }

    // printf("compression finished at rank %d\n", rank);

    if ((relative_rank % 2))
    {
        SZ_fast_decompress_split(4, SZ_FLOAT, recvbuf,
                                 (char *)outputBytes, 0, 0, 0, 0, recvcount);
    }
    if ((rank == root) && (root == 0) && (recvbuf != MPI_IN_PLACE))
    {
        /* for root=0, put root's data in recvbuf if not MPI_IN_PLACE */
        memcpy(recvbuf, sendbuf, sendcount * sizeof(sendtype));
    }
    else if (!(relative_rank % 2) && (recvbuf != MPI_IN_PLACE))
    {
        /* for non-zero root and non-leaf nodes, copy from tmp_buf
         * into recvbuf */
        SZ_fast_decompress_split(4, SZ_FLOAT, recvbuf,
                                 (char *)outputBytes, 0, 0, 0, 0, recvcount);
        // memcpy(recvbuf, tmp_buf, nbytes);
    }

    free(outputBytes);
    free(compressed_lefts);
    free(compressed_offsets);
    free(compressed_sizes);
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
    free(tmp_buf);
    if (mpi_errno_ret)
        mpi_errno = mpi_errno_ret;
    return mpi_errno;
fn_fail:
    goto fn_exit;
}