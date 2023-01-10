/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */
#include <stdio.h>
#include "mpi.h"
// #include "./utils.h"
#include "szx.h"
#include "szx_rw.h"

#include "./include/libs.h"

// 10MB 2 4 8 16 latencies internode

int MPIR_Allgatherv_intra_ring(const void *sendbuf,
                               MPI_Aint sendcount,
                               MPI_Datatype sendtype,
                               void *recvbuf,
                               const MPI_Aint *recvcounts,
                               const MPI_Aint *displs,
                               MPI_Datatype recvtype,
                               MPI_Comm comm)
{
    int comm_size, rank, i, left, right;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;
    MPI_Status status;
    MPI_Aint recvtype_extent;
    int total_count;
    recvtype_extent = sizeof(recvtype);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    total_count = 0;
    for (i = 0; i < comm_size; i++)
        total_count += recvcounts[i];

    if (total_count == 0)
        goto fn_exit;

    if (sendbuf != MPI_IN_PLACE)
    {
        /* First, load the "local" version in the recvbuf. */
        // mpi_errno = MPIR_Localcopy(sendbuf, sendcount, sendtype,
        //                            ((char *)recvbuf + displs[rank] * recvtype_extent),
        //                            recvcounts[rank], recvtype);
        memcpy((char *)recvbuf + displs[rank] * sizeof(recvtype),
               sendbuf, sizeof(data_type) * sendcount);
    }

    left = (comm_size + rank - 1) % comm_size;
    right = (rank + 1) % comm_size;

    MPI_Aint torecv, tosend, max, chunk_count;
    torecv = total_count - recvcounts[rank];
    tosend = total_count - recvcounts[right];

    chunk_count = 0;
    max = recvcounts[0];
    for (i = 1; i < comm_size; i++)
        if (max < recvcounts[i])
            max = recvcounts[i];
    if (MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE > 0 &&
        max * recvtype_extent > MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE)
    {
        chunk_count = MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE / recvtype_extent;
        /* Handle the case where the datatype extent is larger than
         * the pipeline size. */
        if (!chunk_count)
            chunk_count = 1;
    }
    /* pipeline is disabled */
    if (!chunk_count)
        chunk_count = max;

    int soffset, roffset;
    int sidx, ridx;
    sidx = rank;
    ridx = left;
    soffset = 0;
    roffset = 0;
    while (tosend || torecv)
    { /* While we have data to send or receive */
        MPI_Aint sendnow, recvnow;
        sendnow = ((recvcounts[sidx] - soffset) >
                   chunk_count)
                      ? chunk_count
                      : (recvcounts[sidx] - soffset);
        recvnow = ((recvcounts[ridx] - roffset) >
                   chunk_count)
                      ? chunk_count
                      : (recvcounts[ridx] - roffset);

        char *sbuf, *rbuf;
        sbuf = (char *)recvbuf + ((displs[sidx] + soffset) * recvtype_extent);
        rbuf = (char *)recvbuf + ((displs[ridx] + roffset) * recvtype_extent);

        /* Protect against wrap-around of indices */
        if (!tosend)
            sendnow = 0;
        if (!torecv)
            recvnow = 0;

        /* Communicate */
        if (!sendnow && !recvnow)
        {
            /* Don't do anything. This case is possible if two
             * consecutive processes contribute 0 bytes each. */
        }
        else if (!sendnow)
        { /* If there's no data to send, just do a recv call */
            mpi_errno =
                MPI_Recv(rbuf, recvnow, recvtype, left, MPIR_ALLGATHERV_TAG, comm, &status);
            if (mpi_errno)
            {
                exit(-1);
            }
            torecv -= recvnow;
        }
        else if (!recvnow)
        { /* If there's no data to receive, just do a send call */
            mpi_errno =
                MPI_Send(sbuf, sendnow, recvtype, right, MPIR_ALLGATHERV_TAG, comm);
            if (mpi_errno)
            {
                exit(-1);
            }
            tosend -= sendnow;
        }
        else
        { /* There's data to be sent and received */
            mpi_errno = MPI_Sendrecv(sbuf, sendnow, recvtype, right, MPIR_ALLGATHERV_TAG,
                                     rbuf, recvnow, recvtype, left, MPIR_ALLGATHERV_TAG,
                                     comm, &status);
            if (mpi_errno)
            {
                exit(-1);
            }
            tosend -= sendnow;
            torecv -= recvnow;
        }

        soffset += sendnow;
        roffset += recvnow;
        if (soffset == recvcounts[sidx])
        {
            soffset = 0;
            sidx = (sidx + comm_size - 1) % comm_size;
        }
        if (roffset == recvcounts[ridx])
        {
            roffset = 0;
            ridx = (ridx + comm_size - 1) % comm_size;
        }
    }
fn_exit:
    return mpi_errno;
fn_fail:
    goto fn_exit;
}

int MPIR_Allreduce_intra_ring(const void *sendbuf,
                              void *recvbuf,
                              MPI_Aint count,
                              MPI_Datatype datatype,
                              MPI_Op op,
                              MPI_Comm comm)
{
    int mpi_errno = MPI_SUCCESS, mpi_errno_ret = MPI_SUCCESS;
    int i, src, dst;
    int nranks, is_inplace, rank;
    size_t extent;
    MPI_Aint lb, true_extent;
    MPI_Aint *cnts, *displs; /* Created for the allgatherv call */
    int send_rank, recv_rank, total_count;
    void *tmpbuf;
    int tag;
    MPI_Request reqs[2]; /* one send and one recv per transfer */

    extent = sizeof(datatype);
    is_inplace = (sendbuf == MPI_IN_PLACE);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);

    cnts = (MPI_Aint *)malloc(nranks * sizeof(MPI_Aint));
    assert(cnts != NULL);
    displs = (MPI_Aint *)malloc(nranks * sizeof(MPI_Aint));
    assert(displs != NULL);

    for (i = 0; i < nranks; i++)
        cnts[i] = 0;

    total_count = 0;
    for (i = 0; i < nranks; i++)
    {
        cnts[i] = (count + nranks - 1) / nranks;
        if (total_count + cnts[i] > count)
        {
            cnts[i] = count - total_count;
            break;
        }
        else
            total_count += cnts[i];
    }

    displs[0] = 0;
    for (i = 1; i < nranks; i++)
        displs[i] = displs[i - 1] + cnts[i - 1];

    /* Phase 1: copy to tmp buf */
    if (!is_inplace)
    {
        memcpy(recvbuf, sendbuf, sizeof(datatype) * count);
    }

    /* Phase 2: Ring based send recv reduce scatter */
    /* Need only 2 spaces for current and previous reduce_id(s) */
    tmpbuf = (void *)malloc(count * extent);

    src = (nranks + rank - 1) % nranks;
    dst = (rank + 1) % nranks;

    for (i = 0; i < nranks - 1; i++)
    {
        recv_rank = (nranks + rank - 2 - i) % nranks;
        send_rank = (nranks + rank - 1 - i) % nranks;

        /* get a new tag to prevent out of order messages */
        tag = tag_base;

        mpi_errno = MPI_Irecv(tmpbuf, cnts[recv_rank], datatype, src, tag, comm, &reqs[0]);
        if (mpi_errno)
        {
            exit(-1);
        }

        mpi_errno = MPI_Isend((char *)recvbuf + displs[send_rank] * extent, cnts[send_rank],
                              datatype, dst, tag, comm, &reqs[1]);
        if (mpi_errno)
        {
            exit(-1);
        }

        mpi_errno = MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);
        if (mpi_errno)
        {
            exit(-1);
        }

        mpi_errno =
            MPI_Reduce_local(tmpbuf, (char *)recvbuf + displs[recv_rank] * extent,
                             cnts[recv_rank], datatype, op);
        if (mpi_errno)
        {
            exit(-1);
        }
    }

    /* Phase 3: Allgatherv ring, so everyone has the reduced data */
    mpi_errno = MPIR_Allgatherv_intra_ring(MPI_IN_PLACE, -1, MPI_DATATYPE_NULL, recvbuf, cnts,
                                           displs, datatype, comm);
    if (mpi_errno)
    {
        exit(-1);
    }

    free(cnts);
    free(displs);
    free(tmpbuf);
    return mpi_errno;
}

int MPIR_Allgatherv_intra_ring_RI(const void *sendbuf,
                                  MPI_Aint sendcount,
                                  MPI_Datatype sendtype,
                                  void *recvbuf,
                                  const MPI_Aint *recvcounts,
                                  const MPI_Aint *displs,
                                  MPI_Datatype recvtype,
                                  MPI_Comm comm,
                                  float compressionRatio,
                                  float tolerance)
{
    int comm_size, rank, i, left, right;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;
    MPI_Status status;
    MPI_Aint recvtype_extent;
    int total_count;
    recvtype_extent = sizeof(recvtype);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    // Compression variables
    size_t outSize;
    size_t byteLength;

    total_count = 0;
    for (i = 0; i < comm_size; i++)
        total_count += recvcounts[i];

    if (total_count == 0)
        goto fn_exit;

    if (sendbuf != MPI_IN_PLACE)
    {
        /* First, load the "local" version in the recvbuf. */
        // mpi_errno = MPIR_Localcopy(sendbuf, sendcount, sendtype,
        //                            ((char *)recvbuf + displs[rank] * recvtype_extent),
        //                            recvcounts[rank], recvtype);
        memcpy((char *)recvbuf + displs[rank] * sizeof(recvtype),
               sendbuf, sizeof(data_type) * sendcount);
    }

    left = (comm_size + rank - 1) % comm_size;
    right = (rank + 1) % comm_size;

    MPI_Aint torecv, tosend, max, chunk_count;
    torecv = total_count - recvcounts[rank];
    tosend = total_count - recvcounts[right];

    chunk_count = 0;
    max = recvcounts[0];
    for (i = 1; i < comm_size; i++)
        if (max < recvcounts[i])
            max = recvcounts[i];
    if (MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE > 0 &&
        max * recvtype_extent > MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE)
    {
        chunk_count = MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE / recvtype_extent;
        /* Handle the case where the datatype extent is larger than
         * the pipeline size. */
        if (!chunk_count)
            chunk_count = 1;
    }
    /* pipeline is disabled */
    if (!chunk_count)
        chunk_count = max;

    int soffset, roffset;
    int sidx, ridx;
    sidx = rank;
    ridx = left;
    soffset = 0;
    roffset = 0;
    while (tosend || torecv)
    { /* While we have data to send or receive */
        MPI_Aint sendnow, recvnow;
        sendnow = ((recvcounts[sidx] - soffset) >
                   chunk_count)
                      ? chunk_count
                      : (recvcounts[sidx] - soffset);
        recvnow = ((recvcounts[ridx] - roffset) >
                   chunk_count)
                      ? chunk_count
                      : (recvcounts[ridx] - roffset);

        char *sbuf, *rbuf;
        sbuf = (char *)recvbuf + ((displs[sidx] + soffset) * recvtype_extent);
        rbuf = (char *)recvbuf + ((displs[ridx] + roffset) * recvtype_extent);

        /* Protect against wrap-around of indices */
        if (!tosend)
            sendnow = 0;
        if (!torecv)
            recvnow = 0;

        /* Communicate */
        if (!sendnow && !recvnow)
        {
            /* Don't do anything. This case is possible if two
             * consecutive processes contribute 0 bytes each. */
        }
        else if (!sendnow)
        { /* If there's no data to send, just do a recv call */
            // mpi_errno =
            //     MPI_Recv(rbuf, recvnow, recvtype, left, MPIR_ALLGATHERV_TAG, comm, &status);

            mpi_errno =
                MPI_Recv(rbuf, recvnow, recvtype, left, MPIR_ALLGATHERV_TAG, comm, &status);
            if (mpi_errno)
            {
                exit(-1);
            }
            MPI_Get_count(&status, MPI_FLOAT, &byteLength);
            unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, rbuf, byteLength, 0, 0, 0, 0, recvnow);
            memcpy(rbuf, data, sizeof(recvtype) * recvnow);
            free(data);
            torecv -= recvnow;
        }
        else if (!recvnow)
        { /* If there's no data to receive, just do a send call */
            // mpi_errno =
            //     MPI_Send(sbuf, sendnow, recvtype, right, MPIR_ALLGATHERV_TAG, comm);

            unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, sbuf, &outSize,
                                                         cpr_mode, absErrBound, relBoundRatio,
                                                         compressionRatio, tolerance,
                                                         0, 0, 0, 0, sendnow);
            mpi_errno =
                MPI_Send((void *)bytes, outSize / sizeof(recvtype), recvtype, right, MPIR_ALLGATHERV_TAG, comm);
            if (mpi_errno)
            {
                exit(-1);
            }
            free(bytes);
            tosend -= sendnow;
        }
        else
        { /* There's data to be sent and received */
            // mpi_errno = MPI_Sendrecv(sbuf, sendnow, recvtype, right, MPIR_ALLGATHERV_TAG,
            //                          rbuf, recvnow, recvtype, left, MPIR_ALLGATHERV_TAG,
            //                          comm, &status);

            unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, sbuf, &outSize, cpr_mode, absErrBound, relBoundRatio,
                                                         compressionRatio, tolerance, 0, 0, 0, 0, sendnow);
            mpi_errno = MPI_Sendrecv((void *)bytes, outSize / sizeof(recvtype), recvtype, right, MPIR_ALLGATHERV_TAG,
                                     rbuf, recvnow, recvtype, left, MPIR_ALLGATHERV_TAG,
                                     comm, &status);
            if (mpi_errno)
            {
                exit(-1);
            }
            free(bytes);
            MPI_Get_count(&status, MPI_FLOAT, &byteLength);
            unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, rbuf, byteLength, 0, 0, 0, 0, recvnow);
            memcpy(rbuf, data, sizeof(recvtype) * recvnow);
            free(data);
            tosend -= sendnow;
            torecv -= recvnow;
        }

        soffset += sendnow;
        roffset += recvnow;
        if (soffset == recvcounts[sidx])
        {
            soffset = 0;
            sidx = (sidx + comm_size - 1) % comm_size;
        }
        if (roffset == recvcounts[ridx])
        {
            roffset = 0;
            ridx = (ridx + comm_size - 1) % comm_size;
        }
    }
fn_exit:
    return mpi_errno;
fn_fail:
    goto fn_exit;
}

int MPI_Allreduce_SZx_FXR_RI(const void *sendbuf,
                             void *recvbuf,
                             float compressionRatio,
                             float tolerance,
                             MPI_Aint count,
                             MPI_Datatype datatype,
                             MPI_Op op,
                             MPI_Comm comm)
{
    int mpi_errno = MPI_SUCCESS, mpi_errno_ret = MPI_SUCCESS;
    int i, src, dst;
    int nranks, is_inplace, rank;
    size_t extent;
    MPI_Aint lb, true_extent;
    MPI_Aint *cnts, *displs; /* Created for the allgatherv call */
    int send_rank, recv_rank, total_count;
    int tmp_len = (count - 1) / (compressionRatio * (1 - tolerance)) + 1;
    void *tmpbuf;
    int tag;
    MPI_Request reqs[2]; /* one send and one recv per transfer */
    MPI_Status stas[2];
    extent = sizeof(datatype);
    is_inplace = (sendbuf == MPI_IN_PLACE);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);
    // Compression variables
    size_t outSize;
    size_t byteLength;

    cnts = (MPI_Aint *)malloc(nranks * sizeof(MPI_Aint));
    assert(cnts != NULL);
    displs = (MPI_Aint *)malloc(nranks * sizeof(MPI_Aint));
    assert(displs != NULL);

    for (i = 0; i < nranks; i++)
        cnts[i] = 0;

    total_count = 0;
    for (i = 0; i < nranks; i++)
    {
        cnts[i] = (count + nranks - 1) / nranks;
        if (total_count + cnts[i] > count)
        {
            cnts[i] = count - total_count;
            break;
        }
        else
            total_count += cnts[i];
    }

    displs[0] = 0;
    for (i = 1; i < nranks; i++)
        displs[i] = displs[i - 1] + cnts[i - 1];

    /* Phase 1: copy to tmp buf */
    if (!is_inplace)
    {
        memcpy(recvbuf, sendbuf, sizeof(datatype) * count);
    }

    /* Phase 2: Ring based send recv reduce scatter */
    /* Need only 2 spaces for current and previous reduce_id(s) */

    tmpbuf = (void *)malloc(tmp_len * extent);

    src = (nranks + rank - 1) % nranks;
    dst = (rank + 1) % nranks;

    for (i = 0; i < nranks - 1; i++)
    {
        recv_rank = (nranks + rank - 2 - i) % nranks;
        send_rank = (nranks + rank - 1 - i) % nranks;

        /* get a new tag to prevent out of order messages */
        tag = tag_base;

        // mpi_errno = MPI_Irecv(tmpbuf, cnts[recv_rank], datatype, src, tag, comm, &reqs[0]);

        mpi_errno = MPI_Irecv(tmpbuf, cnts[recv_rank], datatype, src, tag, comm, &reqs[0]);

        if (mpi_errno)
        {
            exit(-1);
        }
        unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, (char *)recvbuf + displs[send_rank] * extent,
                                                     &outSize, cpr_mode, absErrBound, relBoundRatio,
                                                     compressionRatio, tolerance, 0, 0, 0, 0, cnts[send_rank]);
        mpi_errno = MPI_Isend(bytes, outSize / sizeof(datatype),
                              datatype, dst, tag, comm, &reqs[1]);
        if (mpi_errno)
        {
            exit(-1);
        }

        mpi_errno = MPI_Waitall(2, reqs, stas);
        if (mpi_errno)
        {
            exit(-1);
        }
        MPI_Get_count(&stas[0], MPI_FLOAT, &byteLength);
        unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, tmpbuf, byteLength, 0, 0, 0, 0, cnts[recv_rank]);

        mpi_errno =
            MPI_Reduce_local(data, (char *)recvbuf + displs[recv_rank] * extent,
                             cnts[recv_rank], datatype, op);
        if (mpi_errno)
        {
            exit(-1);
        }
        free(bytes);
        free(data);
    }

    /* Phase 3: Allgatherv ring, so everyone has the reduced data */
    mpi_errno = MPIR_Allgatherv_intra_ring_RI(MPI_IN_PLACE, -1, MPI_DATATYPE_NULL, recvbuf, cnts,
                                              displs, datatype, comm, compressionRatio, tolerance);
    if (mpi_errno)
    {
        exit(-1);
    }

    free(cnts);
    free(displs);
    free(tmpbuf);
    return mpi_errno;
}

int MPIR_Allgatherv_intra_ring_RI_record(const void *sendbuf,
                                         MPI_Aint sendcount,
                                         MPI_Datatype sendtype,
                                         void *recvbuf,
                                         const MPI_Aint *recvcounts,
                                         const MPI_Aint *displs,
                                         MPI_Datatype recvtype,
                                         MPI_Comm comm,
                                         float compressionRatio,
                                         float tolerance,
                                         double *MPI_timer,
                                         double *CPR_timer)
{
    int comm_size, rank, i, left, right;
    int mpi_errno = MPI_SUCCESS;
    int mpi_errno_ret = MPI_SUCCESS;
    MPI_Status status;
    MPI_Aint recvtype_extent;
    int total_count;
    recvtype_extent = sizeof(recvtype);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    // Compression variables
    size_t outSize;
    size_t byteLength;

    total_count = 0;
    for (i = 0; i < comm_size; i++)
        total_count += recvcounts[i];

    if (total_count == 0)
        goto fn_exit;

    if (sendbuf != MPI_IN_PLACE)
    {
        /* First, load the "local" version in the recvbuf. */
        // mpi_errno = MPIR_Localcopy(sendbuf, sendcount, sendtype,
        //                            ((char *)recvbuf + displs[rank] * recvtype_extent),
        //                            recvcounts[rank], recvtype);
        memcpy((char *)recvbuf + displs[rank] * sizeof(recvtype),
               sendbuf, sizeof(data_type) * sendcount);
    }

    left = (comm_size + rank - 1) % comm_size;
    right = (rank + 1) % comm_size;

    MPI_Aint torecv, tosend, max, chunk_count;
    torecv = total_count - recvcounts[rank];
    tosend = total_count - recvcounts[right];

    chunk_count = 0;
    max = recvcounts[0];
    for (i = 1; i < comm_size; i++)
        if (max < recvcounts[i])
            max = recvcounts[i];
    if (MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE > 0 &&
        max * recvtype_extent > MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE)
    {
        chunk_count = MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE / recvtype_extent;
        /* Handle the case where the datatype extent is larger than
         * the pipeline size. */
        if (!chunk_count)
            chunk_count = 1;
    }
    /* pipeline is disabled */
    if (!chunk_count)
        chunk_count = max;

    int soffset, roffset;
    int sidx, ridx;
    sidx = rank;
    ridx = left;
    soffset = 0;
    roffset = 0;
    while (tosend || torecv)
    { /* While we have data to send or receive */
        MPI_Aint sendnow, recvnow;
        sendnow = ((recvcounts[sidx] - soffset) >
                   chunk_count)
                      ? chunk_count
                      : (recvcounts[sidx] - soffset);
        recvnow = ((recvcounts[ridx] - roffset) >
                   chunk_count)
                      ? chunk_count
                      : (recvcounts[ridx] - roffset);

        char *sbuf, *rbuf;
        sbuf = (char *)recvbuf + ((displs[sidx] + soffset) * recvtype_extent);
        rbuf = (char *)recvbuf + ((displs[ridx] + roffset) * recvtype_extent);

        /* Protect against wrap-around of indices */
        if (!tosend)
            sendnow = 0;
        if (!torecv)
            recvnow = 0;

        /* Communicate */
        if (!sendnow && !recvnow)
        {
            /* Don't do anything. This case is possible if two
             * consecutive processes contribute 0 bytes each. */
        }
        else if (!sendnow)
        { /* If there's no data to send, just do a recv call */
            // mpi_errno =
            //     MPI_Recv(rbuf, recvnow, recvtype, left, MPIR_ALLGATHERV_TAG, comm, &status);
            *MPI_timer -= MPI_Wtime();
            mpi_errno =
                MPI_Recv(rbuf, recvnow, recvtype, left, MPIR_ALLGATHERV_TAG, comm, &status);
            *MPI_timer += MPI_Wtime();
            if (mpi_errno)
            {
                exit(-1);
            }
            MPI_Get_count(&status, MPI_FLOAT, &byteLength);
            *CPR_timer -= MPI_Wtime();
            unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, rbuf, byteLength, 0, 0, 0, 0, recvnow);
            memcpy(rbuf, data, sizeof(recvtype) * recvnow);
            free(data);
            *CPR_timer += MPI_Wtime();
            torecv -= recvnow;
        }
        else if (!recvnow)
        { /* If there's no data to receive, just do a send call */
            // mpi_errno =
            //     MPI_Send(sbuf, sendnow, recvtype, right, MPIR_ALLGATHERV_TAG, comm);
            *CPR_timer -= MPI_Wtime();
            unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, sbuf, &outSize,
                                                         cpr_mode, absErrBound, relBoundRatio,
                                                         compressionRatio, tolerance,
                                                         0, 0, 0, 0, sendnow);
            *CPR_timer += MPI_Wtime();
            *MPI_timer -= MPI_Wtime();
            mpi_errno =
                MPI_Send((void *)bytes, outSize / sizeof(recvtype), recvtype, right, MPIR_ALLGATHERV_TAG, comm);
            *MPI_timer += MPI_Wtime();
            if (mpi_errno)
            {
                exit(-1);
            }
            free(bytes);
            tosend -= sendnow;
        }
        else
        { /* There's data to be sent and received */
            // mpi_errno = MPI_Sendrecv(sbuf, sendnow, recvtype, right, MPIR_ALLGATHERV_TAG,
            //                          rbuf, recvnow, recvtype, left, MPIR_ALLGATHERV_TAG,
            //                          comm, &status);
            *CPR_timer -= MPI_Wtime();
            unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, sbuf, &outSize, cpr_mode, absErrBound, relBoundRatio,
                                                         compressionRatio, tolerance, 0, 0, 0, 0, sendnow);
            *CPR_timer += MPI_Wtime();
            *MPI_timer -= MPI_Wtime();
            mpi_errno = MPI_Sendrecv((void *)bytes, outSize / sizeof(recvtype), recvtype, right, MPIR_ALLGATHERV_TAG,
                                     rbuf, recvnow, recvtype, left, MPIR_ALLGATHERV_TAG,
                                     comm, &status);
            *MPI_timer += MPI_Wtime();
            if (mpi_errno)
            {
                exit(-1);
            }
            free(bytes);
            MPI_Get_count(&status, MPI_FLOAT, &byteLength);
            *CPR_timer -= MPI_Wtime();
            unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, rbuf, byteLength, 0, 0, 0, 0, recvnow);
            memcpy(rbuf, data, sizeof(recvtype) * recvnow);
            free(data);
            *CPR_timer += MPI_Wtime();
            tosend -= sendnow;
            torecv -= recvnow;
        }

        soffset += sendnow;
        roffset += recvnow;
        if (soffset == recvcounts[sidx])
        {
            soffset = 0;
            sidx = (sidx + comm_size - 1) % comm_size;
        }
        if (roffset == recvcounts[ridx])
        {
            roffset = 0;
            ridx = (ridx + comm_size - 1) % comm_size;
        }
    }
fn_exit:
    return mpi_errno;
fn_fail:
    goto fn_exit;
}

int MPI_Allreduce_SZx_FXR_RI_record(const void *sendbuf,
                                    void *recvbuf,
                                    float compressionRatio,
                                    float tolerance,
                                    MPI_Aint count,
                                    MPI_Datatype datatype,
                                    MPI_Op op,
                                    MPI_Comm comm)
{
    int mpi_errno = MPI_SUCCESS, mpi_errno_ret = MPI_SUCCESS;
    int i, src, dst;
    int nranks, is_inplace, rank;
    size_t extent;
    MPI_Aint lb, true_extent;
    MPI_Aint *cnts, *displs; /* Created for the allgatherv call */
    int send_rank, recv_rank, total_count;
    int tmp_len = 0;
    if (cpr_mode == FXR)
    {
        tmp_len = (count - 1) / (compressionRatio * (1 - tolerance)) + 1;
    }
    if (cpr_mode == ABS)
    {
        tmp_len = count;
    }
    void *tmpbuf;
    int tag;
    MPI_Request reqs[2]; /* one send and one recv per transfer */
    MPI_Status stas[2];
    extent = sizeof(datatype);
    is_inplace = (sendbuf == MPI_IN_PLACE);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);
    double MPI_timer = 0.0;
    double CPR_timer = 0.0;
    double CPT_timer = 0.0;

    // Compression variables
    size_t outSize;
    size_t byteLength;

    cnts = (MPI_Aint *)malloc(nranks * sizeof(MPI_Aint));
    assert(cnts != NULL);
    displs = (MPI_Aint *)malloc(nranks * sizeof(MPI_Aint));
    assert(displs != NULL);

    for (i = 0; i < nranks; i++)
        cnts[i] = 0;

    total_count = 0;
    for (i = 0; i < nranks; i++)
    {
        cnts[i] = (count + nranks - 1) / nranks;
        if (total_count + cnts[i] > count)
        {
            cnts[i] = count - total_count;
            break;
        }
        else
            total_count += cnts[i];
    }

    displs[0] = 0;
    for (i = 1; i < nranks; i++)
        displs[i] = displs[i - 1] + cnts[i - 1];

    /* Phase 1: copy to tmp buf */
    if (!is_inplace)
    {
        memcpy(recvbuf, sendbuf, sizeof(datatype) * count);
    }

    /* Phase 2: Ring based send recv reduce scatter */
    /* Need only 2 spaces for current and previous reduce_id(s) */

    tmpbuf = (void *)malloc(tmp_len * extent);

    src = (nranks + rank - 1) % nranks;
    dst = (rank + 1) % nranks;

    for (i = 0; i < nranks - 1; i++)
    {
        recv_rank = (nranks + rank - 2 - i) % nranks;
        send_rank = (nranks + rank - 1 - i) % nranks;

        /* get a new tag to prevent out of order messages */
        tag = tag_base;

        // mpi_errno = MPI_Irecv(tmpbuf, cnts[recv_rank], datatype, src, tag, comm, &reqs[0]);
        MPI_timer -= MPI_Wtime();
        mpi_errno = MPI_Irecv(tmpbuf, cnts[recv_rank], datatype, src, tag, comm, &reqs[0]);
        MPI_timer += MPI_Wtime();
        if (mpi_errno)
        {
            exit(-1);
        }
        CPR_timer -= MPI_Wtime();
        unsigned char *bytes = SZ_fast_compress_args(fast_mode, SZ_FLOAT, (char *)recvbuf + displs[send_rank] * extent,
                                                     &outSize, cpr_mode, absErrBound, relBoundRatio, compressionRatio, tolerance, 0, 0, 0, 0, cnts[send_rank]);
        CPR_timer += MPI_Wtime();
        MPI_timer -= MPI_Wtime();
        mpi_errno = MPI_Isend(bytes, outSize / sizeof(datatype),
                              datatype, dst, tag, comm, &reqs[1]);
        MPI_timer += MPI_Wtime();
        if (mpi_errno)
        {
            exit(-1);
        }
        MPI_timer -= MPI_Wtime();
        mpi_errno = MPI_Waitall(2, reqs, stas);
        MPI_timer += MPI_Wtime();
        if (mpi_errno)
        {
            exit(-1);
        }
        MPI_Get_count(&stas[0], MPI_FLOAT, &byteLength);
        CPR_timer -= MPI_Wtime();
        unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, tmpbuf, byteLength, 0, 0, 0, 0, cnts[recv_rank]);
        CPR_timer += MPI_Wtime();
        CPT_timer -= MPI_Wtime();
        mpi_errno =
            MPI_Reduce_local(data, (char *)recvbuf + displs[recv_rank] * extent,
                             cnts[recv_rank], datatype, op);
        CPT_timer += MPI_Wtime();
        if (mpi_errno)
        {
            exit(-1);
        }
        free(bytes);
        free(data);
    }

    /* Phase 3: Allgatherv ring, so everyone has the reduced data */
    mpi_errno = MPIR_Allgatherv_intra_ring_RI_record(MPI_IN_PLACE, -1, MPI_DATATYPE_NULL, recvbuf, cnts,
                                                     displs, datatype, comm, compressionRatio, tolerance, &MPI_timer, &CPR_timer);
    if (mpi_errno)
    {
        exit(-1);
    }

    free(cnts);
    free(displs);
    free(tmpbuf);
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
        double avg_total_time = CPR_timer_all / nranks + MPI_timer_all / nranks + CPT_timer / nranks;
        printf("For whole allreduce, the avg CPR time is %f, the avg MPI time is %f, the avg CPT time is %f, the avg total time is %f\n", CPR_timer_all / nranks * 1000000,
               MPI_timer_all / nranks * 1000000, CPT_timer / nranks * 1000000, avg_total_time * 1000000);
    }
    return mpi_errno;
}