/**
 *  @file szx_rw.h
 *  @author Sheng Di
 *  @date Jan, 2022
 *  @brief Header file for the whole io interface.
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _SZX_RW_H
#define _SZX_RW_H

#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

int SZx_checkFileExistance(char* filePath);

size_t SZx_checkFileSize(char *srcFilePath, int *status);

unsigned char *SZx_readByteData(char *srcFilePath, size_t *byteLength, int *status);
double *SZx_readDoubleData(char *srcFilePath, size_t *nbEle, int *status);
float *SZx_readFloatData(char *srcFilePath, size_t *nbEle, int *status);

double *SZx_readDoubleData_systemEndian(char *srcFilePath, size_t *nbEle, int *status);
float *SZx_readFloatData_systemEndian(char *srcFilePath, size_t *nbEle, int *status);

void SZx_writeByteData(unsigned char *bytes, size_t byteLength, char *tgtFilePath, int *status);
void SZx_writeDoubleData(double *data, size_t nbEle, char *tgtFilePath, int *status);
void SZx_writeFloatData(float *data, size_t nbEle, char *tgtFilePath, int *status);
void SZx_writeIntData(int *data, size_t nbEle, char *tgtFilePath, int *status);
void SZx_writeData(void *data, int dataType, size_t nbEle, char *tgtFilePath, int *status);
void SZx_writeFloatData_inBytes(float *data, size_t nbEle, char* tgtFilePath, int *status);
void SZx_writeDoubleData_inBytes(double *data, size_t nbEle, char* tgtFilePath, int *status);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _SZX_RW_H  ----- */
