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

int MPIR_Allgatherv_intra_ring_RI2_op_record(const void *sendbuf,
                                             MPI_Aint sendcount,
                                             MPI_Datatype sendtype,
                                             void *recvbuf,
                                             const MPI_Aint *recvcounts,
                                             const MPI_Aint *displs,
                                             MPI_Datatype recvtype,
                                             MPI_Comm comm,
                                             unsigned char *outputBytes,
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
    if (MPIR_CVAR_ALLGATHERV_CPR_PIPELINE_MSG_SIZE > 0 &&
        max * recvtype_extent > MPIR_CVAR_ALLGATHERV_CPR_PIPELINE_MSG_SIZE)
    {
        chunk_count = MPIR_CVAR_ALLGATHERV_CPR_PIPELINE_MSG_SIZE / recvtype_extent;
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

    //prepare decompress buffer
    float *newData = (float *)malloc(chunk_count * recvtype_extent);

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
            // MPI_Get_count(&status, MPI_FLOAT, &byteLength);
            *CPR_timer -= MPI_Wtime();
            SZ_fast_decompress_split(fast_mode, SZ_FLOAT, newData, rbuf, 0, 0, 0, 0, recvnow);
            // unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, rbuf, byteLength, 0, 0, 0, 0, recvnow);
            memcpy(rbuf, newData, sizeof(recvtype) * recvnow);
            // free(data);
            *CPR_timer += MPI_Wtime();
            torecv -= recvnow;
        }
        else if (!recvnow)
        { /* If there's no data to receive, just do a send call */
            // mpi_errno =
            //     MPI_Send(sbuf, sendnow, recvtype, right, MPIR_ALLGATHERV_TAG, comm);
            *CPR_timer -= MPI_Wtime();
            SZ_fast_compress_args2(fast_mode, SZ_FLOAT, sbuf, &outSize,
                                                         outputBytes,
                                                         cpr_mode, absErrBound, relBoundRatio,
                                                         compressionRatio, tolerance,
                                                         0, 0, 0, 0, sendnow);
            unsigned char *bytes =  outputBytes;
            *CPR_timer += MPI_Wtime();
            *MPI_timer -= MPI_Wtime();
            mpi_errno =
                MPI_Send((void *)bytes, outSize / sizeof(recvtype), recvtype, right, MPIR_ALLGATHERV_TAG, comm);
            *MPI_timer += MPI_Wtime();
            if (mpi_errno)
            {
                exit(-1);
            }
            // free(bytes);
            tosend -= sendnow;
        }
        else
        { /* There's data to be sent and received */
            // mpi_errno = MPI_Sendrecv(sbuf, sendnow, recvtype, right, MPIR_ALLGATHERV_TAG,
            //                          rbuf, recvnow, recvtype, left, MPIR_ALLGATHERV_TAG,
            //                          comm, &status);
            *CPR_timer -= MPI_Wtime();
            SZ_fast_compress_args2(fast_mode, SZ_FLOAT, sbuf, &outSize, outputBytes, cpr_mode, absErrBound, relBoundRatio,
                                                         compressionRatio, tolerance, 0, 0, 0, 0, sendnow);
            unsigned char *bytes = outputBytes;
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
            // free(bytes);
            // MPI_Get_count(&status, MPI_FLOAT, &byteLength);
            *CPR_timer -= MPI_Wtime();
            SZ_fast_decompress_split(fast_mode, SZ_FLOAT, newData, rbuf, 0, 0, 0, 0, recvnow);
            // unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, rbuf, byteLength, 0, 0, 0, 0, recvnow);
            memcpy(rbuf, newData, sizeof(recvtype) * recvnow);
            // free(data);
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
    free(newData);
fn_exit:
    return mpi_errno;
fn_fail:
    goto fn_exit;
}

int MPI_Allreduce_SZx_FXR_RI2_op_record(const void *sendbuf,
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
    double ReSca_timer = 0.0;
    double CPR_timer = 0.0;
    double CPT_timer = 0.0;
    double WAIT_timer = 0.0;

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
    // For comparess buffer
    int chunk_size;
    if (fast_mode == 4)
    {
        chunk_size = MULTITHREAD_CHUNK_SIZE;
    }
    else
    {
        chunk_size = SINGLETHREAD_CHUNK_SIZE; 
    }
    int chunk_num = cnts[0] / chunk_size;
    unsigned char *outputBytes = (unsigned char *)malloc(tmp_len * extent + sizeof(size_t) * (chunk_num + 1));
    float *newData = (float *)malloc(cnts[0] * extent);

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

        // int flag;
        // MPI_Test(&reqs[0], &flag, &stas[0]);
        // SZ_fast_compress_args2(fast_mode, SZ_FLOAT, (char *)recvbuf + displs[send_rank] * extent,
        //                                              &outSize, outputBytes, cpr_mode, absErrBound, relBoundRatio, compressionRatio, tolerance, 0, 0, 0, 0, cnts[send_rank]);

        outSize = 0;
        chunk_num = cnts[send_rank] / chunk_size;
        // printf("chunk_num = %d\n", chunk_num);
        int chunk_remainder_size = cnts[send_rank] % chunk_size;
        // printf("chunk_reminder_size = %d\n", chunk_remainder_size);
        int flag, iter;
        size_t chunk_out_size = 0;

        for (iter = 0; iter < chunk_num; iter++)
        {
            MPI_Test(&reqs[0], &flag, &stas[0]);
            // SZ_fast_compress_args2(fast_mode, SZ_FLOAT, (char *)recvbuf + displs[send_rank] * extent + iter * chunk_size * extent,
            //                        &chunk_out_size, outputBytes + outSize,
            //                        cpr_mode, absErrBound, relBoundRatio, compressionRatio, tolerance, 0, 0, 0, 0, chunk_size);
            // outSize += chunk_out_size;
            SZ_fast_compress_args2_split(fast_mode, SZ_FLOAT, (char *)recvbuf + displs[send_rank] * extent + iter * chunk_size * extent,
                                         &chunk_out_size, outputBytes + outSize + sizeof(size_t) * (chunk_num + 1), outputBytes, iter,
                                         cpr_mode, absErrBound, relBoundRatio, compressionRatio, tolerance, 0, 0, 0, 0, chunk_size);
            outSize += chunk_out_size;
            // printf("outSize = %d\n", outSize);
        }
        // printf("Main iterations finished\n");
        // printf("iterations = %d\n", iter);
        // Handling the remainder size
        if (chunk_remainder_size != 0)
        {
            MPI_Test(&reqs[0], &flag, &stas[0]);
            // SZ_fast_compress_args2(fast_mode, SZ_FLOAT, (char *)recvbuf + displs[send_rank] * extent + (cnts[send_rank] - chunk_remainder_size) * extent,
            //                        &chunk_out_size, outputBytes + outSize,
            //                        cpr_mode, absErrBound, relBoundRatio, compressionRatio, tolerance, 0, 0, 0, 0, chunk_remainder_size);
            SZ_fast_compress_args2_split(fast_mode, SZ_FLOAT, (char *)recvbuf + displs[send_rank] * extent + iter * chunk_size * extent,
                                         &chunk_out_size, outputBytes + outSize + sizeof(size_t) * (chunk_num + 1), outputBytes, iter,
                                         cpr_mode, absErrBound, relBoundRatio, compressionRatio, tolerance, 0, 0, 0, 0, chunk_remainder_size);
            outSize += chunk_out_size;
            // printf("outSize = %d\n", outSize);
        }
        // printf("Remainder iterations finished\n");

        // unsigned char *outputBytes_copy = (unsigned char *)malloc(tmp_len * extent);
        // int outSize_copy;
        // SZ_fast_compress_args2(fast_mode, SZ_FLOAT, (char *)recvbuf + displs[send_rank] * extent,
        //                                              &outSize_copy, outputBytes, cpr_mode, absErrBound, relBoundRatio, compressionRatio, tolerance, 0, 0, 0, 0, cnts[send_rank]);
        // printf("compression finished\n");
        unsigned char *bytes = outputBytes;
        CPR_timer += MPI_Wtime();
        MPI_timer -= MPI_Wtime();
        mpi_errno = MPI_Isend(bytes, (outSize + sizeof(size_t) * (chunk_num + 1)) / sizeof(datatype),
                              datatype, dst, tag, comm, &reqs[1]);

        MPI_timer += MPI_Wtime();
        if (mpi_errno)
        {
            exit(-1);
        }
        // MPI_timer -= MPI_Wtime();
        // mpi_errno = MPI_Waitall(2, reqs, stas);
        // MPI_timer += MPI_Wtime();
        if (mpi_errno)
        {
            exit(-1);
        }
        MPI_timer -= MPI_Wtime();
        WAIT_timer -= MPI_Wtime();
        MPI_Wait(&reqs[0], &stas[0]);
        WAIT_timer += MPI_Wtime();
        MPI_timer += MPI_Wtime();

        // MPI_Get_count(&stas[0], MPI_FLOAT, &byteLength);
        CPR_timer -= MPI_Wtime();

        outSize = 0;
        chunk_num = cnts[recv_rank] / chunk_size;
        chunk_remainder_size = cnts[recv_rank] % chunk_size;
        chunk_out_size = 0;
        // printf("data_size = %d chunk_num = %d chunk_reminder_size = %d\n", cnts[send_rank], chunk_num, chunk_remainder_size);
        int decom_offset = 0;
        int decom_out_offset = 0;
        // printf("Decompressing started\n");
        // printf("Chunk_num = %d\n", chunk_num);
        size_t *chunk_arr = (size_t *)tmpbuf;
        unsigned char *cmpBytes;
        for (iter = 0; iter < chunk_num; iter++)
        {
            MPI_Test(&reqs[1], &flag, &stas[1]);
            // SZ_fast_decompress(fast_mode, SZ_FLOAT, tmpbuf, byteLength, 0, 0, 0, 0, cnts[recv_rank]);
            cmpBytes = (unsigned char *)tmpbuf + sizeof(size_t) * (chunk_num + 1) + decom_offset;
            // if (cmpBytes[4] != 128 || iter > 1870)
            // {
            //     printf("iter = %d blocksize %d blocksize address %d\n", iter, cmpBytes[4], &cmpBytes[4]);
            // }

            SZ_fast_decompress_split(fast_mode, SZ_FLOAT, (unsigned char *)newData + decom_out_offset, cmpBytes, 0, 0, 0, 0, chunk_size);
            decom_offset += chunk_arr[iter];
            decom_out_offset += chunk_size * extent;
        }
        // Handling the reminder size
        if (chunk_remainder_size != 0)
        {
            MPI_Test(&reqs[1], &flag, &stas[1]);
            // SZ_fast_decompress(fast_mode, SZ_FLOAT, tmpbuf, byteLength, 0, 0, 0, 0, cnts[recv_rank]);
            cmpBytes = (unsigned char *)tmpbuf + sizeof(size_t) * (chunk_num + 1) + decom_offset;
            // if (cmpBytes[4] != 128)
            // {
            //     printf("iter = %d blocksize %d blocksize address %d\n", iter, cmpBytes[4], &cmpBytes[4]);
            // }
            SZ_fast_decompress_split(fast_mode, SZ_FLOAT, (unsigned char *)newData + decom_out_offset, cmpBytes, 0, 0, 0, 0, chunk_remainder_size);
        }
        // printf("Decompressing finished\n");
        // MPI_Test(&reqs[1], &flag, &stas[1]);
        // unsigned char *data = SZ_fast_decompress(fast_mode, SZ_FLOAT, tmpbuf, byteLength, 0, 0, 0, 0, cnts[recv_rank]);
        CPR_timer += MPI_Wtime();

        CPT_timer -= MPI_Wtime();
        mpi_errno =
            MPI_Reduce_local(newData, (char *)recvbuf + displs[recv_rank] * extent,
                             cnts[recv_rank], datatype, op);
        CPT_timer += MPI_Wtime();

        MPI_timer -= MPI_Wtime();
        WAIT_timer -= MPI_Wtime();
        MPI_Wait(&reqs[1], &stas[1]);
        WAIT_timer += MPI_Wtime();
        MPI_timer += MPI_Wtime();
        if (mpi_errno)
        {
            exit(-1);
        }
        // free(bytes);
        // free(data);
    }

    ReSca_timer = MPI_timer;
    /* Phase 3: Allgatherv ring, so everyone has the reduced data */
    mpi_errno = MPIR_Allgatherv_intra_ring_RI2_op_record(MPI_IN_PLACE, -1, MPI_DATATYPE_NULL, recvbuf, cnts,
                                                         displs, datatype, comm, outputBytes, compressionRatio, tolerance, &MPI_timer, &CPR_timer);
    if (mpi_errno)
    {
        exit(-1);
    }

    free(outputBytes);
    free(newData);
    free(cnts);
    free(displs);
    free(tmpbuf);
    printf("For process %d, the CPR time is %f, the MPI time is %f, the WAIT time is %f, the CPT time is%f\n", rank, CPR_timer, MPI_timer, WAIT_timer, CPT_timer);
    double CPR_timer_all = 0.0;
    double MPI_timer_all = 0.0;
    double CPT_timer_all = 0.0;
    double WAIT_timer_all = 0.0;
    double ReSca_timer_all = 0.0;
    MPI_Reduce(&CPR_timer, &CPR_timer_all, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&MPI_timer, &MPI_timer_all, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&ReSca_timer, &ReSca_timer_all, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&WAIT_timer, &WAIT_timer_all, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&CPT_timer, &CPT_timer_all, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    if (rank == 0)
    {
        double avg_total_time = CPR_timer_all / nranks + MPI_timer_all / nranks + CPT_timer / nranks;
        printf("For whole allreduce, the avg CPR time is %f, the avg MPI time is %f, the avg Reduce_Scatter time is %f, the avg WAIT time is %f, the avg CPT time is %f, the avg total time is %f\n", CPR_timer_all / nranks * 1000000,
               MPI_timer_all / nranks * 1000000, ReSca_timer_all / nranks * 1000000, WAIT_timer_all / nranks * 1000000, CPT_timer / nranks * 1000000, avg_total_time * 1000000);
    }
    return mpi_errno;
}