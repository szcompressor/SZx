/**
 *  @file szx_defines.h
 *  @author Sheng Di
 *  @date Jan, 2022
 *  @brief Header file for the dataCompression.c.
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _SZX_DEFINES_H
#define _SZX_DEFINES_H

#define SZx_VERNUM 0x0200
#define SZx_VER_MAJOR 1
#define SZx_VER_MINOR 1
#define SZx_VER_BUILD 0
#define SZx_VER_REVISION 0

#define ABS 0
#define REL 1
#define FXR 2
#define VR_REL 1  //alternative name to REL
#define ABS_AND_REL 2
#define ABS_OR_REL 3
#define PSNR 4
#define NORM 5
#define FIX_RATE 6

#define PW_REL 10
#define ABS_AND_PW_REL 11
#define ABS_OR_PW_REL 12
#define REL_AND_PW_REL 13
#define REL_OR_PW_REL 14


#define SZx_FLOAT 0
#define SZx_DOUBLE 1
#define SZx_UINT8 2
#define SZx_INT8 3
#define SZx_UINT16 4
#define SZx_INT16 5
#define SZx_UINT32 6
#define SZx_INT32 7
#define SZx_UINT64 8
#define SZx_INT64 9

#define LITTLE_ENDIAN_DATA 0 //refers to the endian type of the data read from the disk
#define BIG_ENDIAN_DATA 1 //big_endian (ppc, max, etc.) ; little_endian (x86, x64, etc.)

#define LITTLE_ENDIAN_SYSTEM 0 //refers to the endian type of the system
#define BIG_ENDIAN_SYSTEM 1


#define SZx_NO_BLOCK_FAST_CMPR 1
#define SZx_WITH_BLOCK_FAST_CMPR 2
#define SZx_RANDOMACCESS_FAST_CMPR 3
#define SZx_OPENMP_FAST_CMPR 4

//SUCCESS returning status
#define SZx_SCES 0  //successful
#define SZx_NSCS -1 //Not successful
#define SZx_FERR -2 //Failed to open input file
#define SZx_TERR -3 //wrong data type (should be only float or double)
#define SZx_DERR -4 //dimension error
#define SZx_MERR -5 //sz_mode error
#define SZx_BERR -6 //bound-mode error (should be only ABS, REL, ABS_AND_REL, ABS_OR_REL, or PW_REL)

#endif /* _SZX_DEFINES_H */
