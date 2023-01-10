/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include <stdio.h>
#include "mpi.h"
#include "./utils.h"
#include "szx.h"
#include "szx_rw.h"

#include "./include/libs.h"
// 10MB 2 4 8 16 latencies internode

/*
 * Algorithm: Recursive Doubling
 *
 * We use this algorithm in the case of user-defined ops because in this case
 * derived datatypes are allowed, and the user could pass basic datatypes on
 * one process and derived on another as long as the type maps are the same.
 * Breaking up derived datatypes to do the reduce-scatter is tricky.
 *
 * Cost = lgp.alpha + n.lgp.beta + n.lgp.gamma
 */


int MPI_Allreduce_compre_RB(const void *sendbuf,
                         void *recvbuf,
                         MPI_Aint count,
                         MPI_Datatype datatype,
                         MPI_Op op,
                         MPI_Comm comm)
{
    // printf("This allreduce does nothing.\n");
    int comm_size, rank;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;
    int mask, dst, is_commutative, pof2, newrank, rem, newdst;
    MPI_Aint true_extent, true_lb, extent;
    // void *tmp_buf;
    void *tmp_buf = NULL;
    tmp_buf = inilize_arr_withoutset(count);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    /* copy local data into recvbuf */
    memcpy(recvbuf, sendbuf, sizeof(data_type) * count);

    /* get nearest power-of-two less than or equal to comm_size */
    pof2 = get_pof2(comm_size);

    rem = comm_size - pof2;

    /* In the non-power-of-two case, all even-numbered
     * processes of rank < 2*rem send their data to
     * (rank+1). These even-numbered processes no longer
     * participate in the algorithm until the very end. The
     * remaining processes form a nice power-of-two. */

    if (rank < 2 * rem)
    {
        if (rank % 2 == 0)
        { /* even */
            mpi_errno = MPI_Send(recvbuf, count,
                                 datatype, rank + 1, ALL_REDUCE_TAG, comm);

            /* temporarily set the rank to -1 so that this
             * process does not pariticipate in recursive
             * doubling */
            newrank = -1;
        }
        else
        { /* odd */
            mpi_errno = MPI_Recv(tmp_buf, count,
                                 datatype, rank - 1,
                                 ALL_REDUCE_TAG, comm, MPI_STATUS_IGNORE);
            /* do the reduction on received data. since the
             * ordering is right, it doesn't matter whether
             * the operation is commutative or not. */
            mpi_errno = MPI_Reduce_local(tmp_buf, recvbuf, count, datatype, op);

            /* change the rank */
            newrank = rank / 2;
        }
    }
    else /* rank >= 2*rem */
        newrank = rank - rem;

    /* If op is user-defined or count is less than pof2, use
     * recursive doubling algorithm. Otherwise do a reduce-scatter
     * followed by allgather. (If op is user-defined,
     * derived datatypes are allowed and the user could pass basic
     * datatypes on one process and derived on another as long as
     * the type maps are the same. Breaking up derived
     * datatypes to do the reduce-scatter is tricky, therefore
     * using recursive doubling in that case.) */

    if (newrank != -1)
    {
        mask = 0x1;
        while (mask < pof2)
        {
            newdst = newrank ^ mask;
            /* find real rank of dest */
            dst = (newdst < rem) ? newdst * 2 + 1 : newdst + rem;

            /* Send the most current data, which is in recvbuf. Recv
             * into tmp_buf */
            // mpi_errno = MPI_Send(recvbuf, count,
            //                       datatype, dst, ALL_REDUCE_TAG, comm);
            // mpi_errno = MPI_Recv(tmp_buf, count,
            //                       datatype, dst,
            //                       ALL_REDUCE_TAG, comm, MPI_STATUS_IGNORE);
            mpi_errno = MPI_Sendrecv(recvbuf, count, datatype,
                                     dst, ALL_REDUCE_TAG, tmp_buf,
                                     count, datatype, dst,
                                     ALL_REDUCE_TAG, comm, MPI_STATUS_IGNORE);

            /* tmp_buf contains data received in this step.
             * recvbuf contains data accumulated so far */

            if ((dst < rank))
            {
                /* op is commutative OR the order is already right */
                mpi_errno = MPI_Reduce_local(tmp_buf, recvbuf, count, datatype, op);
            }
            else
            {
                /* op is noncommutative and the order is not right */
                mpi_errno = MPI_Reduce_local(recvbuf, tmp_buf, count, datatype, op);

                /* copy result back into recvbuf */
                memcpy(recvbuf, tmp_buf, sizeof(data_type) * count);
            }
            mask <<= 1;
        }
    }
    /* In the non-power-of-two case, all odd-numbered
     * processes of rank < 2*rem send the result to
     * (rank-1), the ranks who didn't participate above. */
    if (rank < 2 * rem)
    {
        if (rank % 2) /* odd */
            mpi_errno = MPI_Send(recvbuf, count,
                                 datatype, rank - 1, ALL_REDUCE_TAG, comm);
        else /* even */
            mpi_errno = MPI_Recv(recvbuf, count,
                                 datatype, rank + 1,
                                 ALL_REDUCE_TAG, comm, MPI_STATUS_IGNORE);
    }
    free((void *)(tmp_buf));
    return mpi_errno;
}

int MPI_Allreduce_SZx_FXR_RB(const void *sendbuf,
                             void *recvbuf,
                             float compressionRatio,
                             float tolerance,
                             MPI_Aint count,
                             MPI_Datatype datatype,
                             MPI_Op op,
                             MPI_Comm comm)
{
    // printf("This allreduce does nothing.\n");
    int comm_size, rank;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;
    int mask, dst, is_commutative, pof2, newrank, rem, newdst;
    MPI_Aint true_extent, true_lb, extent;
    int tmp_len = (count - 1) / (compressionRatio * (1 - tolerance)) + 1;
    void *tmp_buf = NULL;
    tmp_buf = inilize_arr_withoutset(tmp_len);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    // Compression variables
    size_t outSize;
    size_t byteLength;
    /* copy local data into recvbuf */
    memcpy(recvbuf, sendbuf, sizeof(data_type) * count);

    /* get nearest power-of-two less than or equal to comm_size */
    pof2 = get_pof2(comm_size);

    rem = comm_size - pof2;

    /* In the non-power-of-two case, all even-numbered
     * processes of rank < 2*rem send their data to
     * (rank+1). These even-numbered processes no longer
     * participate in the algorithm until the very end. The
     * remaining processes form a nice power-of-two. */

    if (rank < 2 * rem)
    {
        if (rank % 2 == 0)
        { /* even */
            // unsigned char * bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, recvbuf, &outSize, ABS, errorBound, 0.001, 0, 0, 0, 0, 0, 0, count);
            unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, recvbuf, &outSize, FXR, 0, 0, compressionRatio, tolerance, 0, 0, 0, 0, count);
            mpi_errno = MPI_Send((void *)bytes, outSize / sizeof(data_type),
                                 datatype, rank + 1, ALL_REDUCE_TAG, comm);

            /* temporarily set the rank to -1 so that this
             * process does not pariticipate in recursive
             * doubling */
            newrank = -1;
            free(bytes);
        }
        else
        { /* odd */
            MPI_Status status;

            mpi_errno = MPI_Recv(tmp_buf, tmp_len,
                                 datatype, rank - 1,
                                 ALL_REDUCE_TAG, comm, &status);
            MPI_Get_count(&status, MPI_FLOAT, &byteLength);
            /* do the reduction on received data. since the
             * ordering is right, it doesn't matter whether
             * the operation is commutative or not. */
            // byteLength = sizeof(data_type) * outSize;
            unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, tmp_buf, byteLength, 0, 0, 0, 0, count);
            mpi_errno = MPI_Reduce_local(data, recvbuf, count, datatype, op);
            /* change the rank */
            newrank = rank / 2;
            free(data);
        }
    }
    else /* rank >= 2*rem */
        newrank = rank - rem;

    /* If op is user-defined or count is less than pof2, use
     * recursive doubling algorithm. Otherwise do a reduce-scatter
     * followed by allgather. (If op is user-defined,
     * derived datatypes are allowed and the user could pass basic
     * datatypes on one process and derived on another as long as
     * the type maps are the same. Breaking up derived
     * datatypes to do the reduce-scatter is tricky, therefore
     * using recursive doubling in that case.) */

    if (newrank != -1)
    {
        mask = 0x1;
        while (mask < pof2)
        {
            newdst = newrank ^ mask;
            /* find real rank of dest */
            dst = (newdst < rem) ? newdst * 2 + 1 : newdst + rem;

            /* Send the most current data, which is in recvbuf. Recv
             * into tmp_buf */
            unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, recvbuf, &outSize, FXR, 0, 0, compressionRatio, tolerance, 0, 0, 0, 0, count);
            // MPI_Request request;

            // printf("The outSize %d on process %d to process %d\n", outSize, rank, dst);
            MPI_Status status;

            mpi_errno = MPI_Sendrecv((void *)bytes, outSize / sizeof(data_type), datatype,
                                     dst, ALL_REDUCE_TAG, tmp_buf,
                                     tmp_len, datatype, dst,
                                     ALL_REDUCE_TAG, comm, &status);
            free(bytes);
            MPI_Get_count(&status, MPI_FLOAT, &byteLength);
            // printf("The byteLength %d on process %d from process %d\n", byteLength, rank, dst);
            unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, tmp_buf, byteLength, 0, 0, 0, 0, count);
            /* tmp_buf contains data received in this step.
             * recvbuf contains data accumulated so far */

            if ((dst < rank))
            {
                /* op is commutative OR the order is already right */
                mpi_errno = MPI_Reduce_local(data, recvbuf, count, datatype, op);
                free(data);
            }
            else
            {
                /* op is noncommutative and the order is not right */
                mpi_errno = MPI_Reduce_local(recvbuf, data, count, datatype, op);

                /* copy result back into recvbuf */
                memcpy(recvbuf, data, sizeof(data_type) * count);
                free(data);
            }
            mask <<= 1;
        }
    }
    /* In the non-power-of-two case, all odd-numbered
     * processes of rank < 2*rem send the result to
     * (rank-1), the ranks who didn't participate above. */
    if (rank < 2 * rem)
    {
        if (rank % 2) /* odd */
        {
            unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, recvbuf, &outSize, FXR, 0, 0, compressionRatio, tolerance, 0, 0, 0, 0, count);
            mpi_errno = MPI_Send((void *)bytes, outSize / sizeof(data_type),
                                 datatype, rank - 1, ALL_REDUCE_TAG, comm);
            // mpi_errno = MPI_Send(recvbuf, count,
            //      datatype, rank - 1, ALL_REDUCE_TAG, comm);
            free(bytes);
        }

        else /* even */
        {
            MPI_Status status;
            mpi_errno = MPI_Recv(tmp_buf, tmp_len,
                                 datatype, rank + 1,
                                 ALL_REDUCE_TAG, comm, &status);
            MPI_Get_count(&status, MPI_FLOAT, &byteLength);
            recvbuf = SZ_fast_decompress(fast_mode, SZ_FLOAT, tmp_buf, byteLength, 0, 0, 0, 0, count);
            // mpi_errno = MPI_Recv(recvbuf, count,
            //                       datatype, rank + 1,
            //                       ALL_REDUCE_TAG, comm, MPI_STATUS_IGNORE);
        }
    }
    // free((void *)data);
    // free((void *)bytes);
    free((void *)tmp_buf);
    return mpi_errno;
}

int MPI_Allreduce_SZx_FXR_RB_record(const void *sendbuf,
                                 void *recvbuf,
                                 float compressionRatio,
                                 float tolerance,
                                 MPI_Aint count,
                                 MPI_Datatype datatype,
                                 MPI_Op op,
                                 MPI_Comm comm)
{
    // printf("This allreduce does nothing.\n");
    int comm_size, rank;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;
    int mask, dst, is_commutative, pof2, newrank, rem, newdst;
    MPI_Aint true_extent, true_lb, extent;
    // int tmp_len = count;
    int tmp_len = (count - 1) / (compressionRatio * (1 - tolerance)) + 1;
    void *tmp_buf = NULL;
    tmp_buf = inilize_arr_withoutset(tmp_len);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);
    double MPI_timer = 0.0;
    double CPR_timer = 0.0;
    double CPT_timer = 0.0;

    // Compression variables
    size_t outSize;
    size_t byteLength;
    /* copy local data into recvbuf */
    memcpy(recvbuf, sendbuf, sizeof(data_type) * count);

    /* get nearest power-of-two less than or equal to comm_size */
    pof2 = get_pof2(comm_size);

    rem = comm_size - pof2;

    /* In the non-power-of-two case, all even-numbered
     * processes of rank < 2*rem send their data to
     * (rank+1). These even-numbered processes no longer
     * participate in the algorithm until the very end. The
     * remaining processes form a nice power-of-two. */

    if (rank < 2 * rem)
    {
        if (rank % 2 == 0)
        { /* even */
            CPR_timer -= MPI_Wtime();
            unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, recvbuf, &outSize, FXR, 0, 0, compressionRatio, tolerance, 0, 0, 0, 0, count);
            CPR_timer += MPI_Wtime();
            MPI_timer -= MPI_Wtime();
            mpi_errno = MPI_Send((void *)bytes, outSize / sizeof(data_type),
                                 datatype, rank + 1, ALL_REDUCE_TAG, comm);
            MPI_timer += MPI_Wtime();

            /* temporarily set the rank to -1 so that this
             * process does not pariticipate in recursive
             * doubling */
            newrank = -1;
            free(bytes);
        }
        else
        { /* odd */
            MPI_Status status;

            MPI_timer -= MPI_Wtime();
            mpi_errno = MPI_Recv(tmp_buf, tmp_len,
                                 datatype, rank - 1,
                                 ALL_REDUCE_TAG, comm, &status);
            MPI_timer += MPI_Wtime();

            MPI_Get_count(&status, MPI_FLOAT, &byteLength);
            /* do the reduction on received data. since the
             * ordering is right, it doesn't matter whether
             * the operation is commutative or not. */
            // byteLength = sizeof(data_type) * outSize;
            CPR_timer -= MPI_Wtime();
            unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, tmp_buf, byteLength, 0, 0, 0, 0, count);
            CPR_timer += MPI_Wtime();
            CPT_timer -= MPI_Wtime();
            mpi_errno = MPI_Reduce_local(data, recvbuf, count, datatype, op);
            CPT_timer += MPI_Wtime();
            /* change the rank */
            newrank = rank / 2;
            free(data);
        }
    }
    else /* rank >= 2*rem */
        newrank = rank - rem;

    /* If op is user-defined or count is less than pof2, use
     * recursive doubling algorithm. Otherwise do a reduce-scatter
     * followed by allgather. (If op is user-defined,
     * derived datatypes are allowed and the user could pass basic
     * datatypes on one process and derived on another as long as
     * the type maps are the same. Breaking up derived
     * datatypes to do the reduce-scatter is tricky, therefore
     * using recursive doubling in that case.) */

    if (newrank != -1)
    {
        mask = 0x1;
        while (mask < pof2)
        {
            newdst = newrank ^ mask;
            /* find real rank of dest */
            dst = (newdst < rem) ? newdst * 2 + 1 : newdst + rem;

            /* Send the most current data, which is in recvbuf. Recv
             * into tmp_buf */
            CPR_timer -= MPI_Wtime();
            unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, recvbuf, &outSize, FXR, 0, 0, compressionRatio, tolerance, 0, 0, 0, 0, count);
            // MPI_Request request;
            CPR_timer += MPI_Wtime();
            // printf("The outSize %d on process %d to process %d\n", outSize, rank, dst);
            MPI_Status status;
            MPI_timer -= MPI_Wtime();
            mpi_errno = MPI_Sendrecv((void *)bytes, outSize / sizeof(data_type), datatype,
                                     dst, ALL_REDUCE_TAG, tmp_buf,
                                     tmp_len, datatype, dst,
                                     ALL_REDUCE_TAG, comm, &status);
            MPI_timer += MPI_Wtime();
            free(bytes);
            MPI_Get_count(&status, MPI_FLOAT, &byteLength);
            // printf("The byteLength %d on process %d from process %d\n", byteLength, rank, dst);
            CPR_timer -= MPI_Wtime();
            unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, tmp_buf, byteLength, 0, 0, 0, 0, count);
            CPR_timer += MPI_Wtime();
            /* tmp_buf contains data received in this step.
             * recvbuf contains data accumulated so far */

            if ((dst < rank))
            {
                /* op is commutative OR the order is already right */
                CPT_timer -= MPI_Wtime();
                mpi_errno = MPI_Reduce_local(data, recvbuf, count, datatype, op);
                CPT_timer += MPI_Wtime();
                free(data);
            }
            else
            {
                /* op is noncommutative and the order is not right */
                CPT_timer -= MPI_Wtime();
                mpi_errno = MPI_Reduce_local(recvbuf, data, count, datatype, op);
                CPT_timer += MPI_Wtime();
                /* copy result back into recvbuf */
                memcpy(recvbuf, data, sizeof(data_type) * count);
                free(data);
            }
            mask <<= 1;
        }
    }
    /* In the non-power-of-two case, all odd-numbered
     * processes of rank < 2*rem send the result to
     * (rank-1), the ranks who didn't participate above. */
    if (rank < 2 * rem)
    {
        if (rank % 2) /* odd */
        {
            CPR_timer -= MPI_Wtime();
            unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, recvbuf, &outSize, FXR, 0, 0, compressionRatio, tolerance, 0, 0, 0, 0, count);
            CPR_timer += MPI_Wtime();
            MPI_timer -= MPI_Wtime();
            mpi_errno = MPI_Send((void *)bytes, outSize / sizeof(data_type),
                                 datatype, rank - 1, ALL_REDUCE_TAG, comm);
            MPI_timer += MPI_Wtime();
            // mpi_errno = MPI_Send(recvbuf, count,
            //      datatype, rank - 1, ALL_REDUCE_TAG, comm);
            free(bytes);
        }

        else /* even */
        {
            MPI_Status status;
            MPI_timer -= MPI_Wtime();
            mpi_errno = MPI_Recv(tmp_buf, tmp_len,
                                 datatype, rank + 1,
                                 ALL_REDUCE_TAG, comm, &status);
            MPI_timer += MPI_Wtime();
            MPI_Get_count(&status, MPI_FLOAT, &byteLength);
            CPR_timer -= MPI_Wtime();
            recvbuf = SZ_fast_decompress(fast_mode, SZ_FLOAT, tmp_buf, byteLength, 0, 0, 0, 0, count);
            CPR_timer += MPI_Wtime();
            // mpi_errno = MPI_Recv(recvbuf, count,
            //                       datatype, rank + 1,
            //                       ALL_REDUCE_TAG, comm, MPI_STATUS_IGNORE);
        }
    }
    // free((void *)data);
    // free((void *)bytes);
    free((void *)tmp_buf);

    printf("For process %d, the CPR time is %f, the MPI time is %f, the CPT time is%f\n", rank, CPR_timer, MPI_timer, CPT_timer);
    double CPR_timer_all = 0.0;
    double MPI_timer_all = 0.0;
    double CPT_timer_all = 0.0;
    MPI_Reduce(&CPR_timer, &CPR_timer_all, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&MPI_timer, &MPI_timer_all, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&CPT_timer, &CPT_timer_all, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    if (rank == 0)
    {
        double avg_total_time = CPR_timer_all / comm_size + MPI_timer_all / comm_size + CPT_timer / comm_size;
        printf("For whole allreduce, the avg CPR time is %f, the avg MPI time is %f, the avg CPT time is %f, the avg total time is %f\n", CPR_timer_all / comm_size * 1000000,
               MPI_timer_all / comm_size * 1000000, CPT_timer / comm_size * 1000000, avg_total_time * 1000000);
    }

    return mpi_errno;
}