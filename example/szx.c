#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "szx.h"
#include "szx_rw.h"
#ifdef _OPENMP
#include "omp.h"
#endif

#define bool _Bool
#define true 1
#define false 0
#define MAXRANGERADIUS 65536

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;



// typedef struct LossyCompressionElement
// {
// 	int leadingZeroBytes; //0,1,2,or 3
// 	unsigned char integerMidBytes[8];
// 	int integerMidBytes_Length; //they are mid_bits actually
// 	//char curBytes[8];
// 	//int curBytes_Length; //4 for single_precision or 8 for double_precision
// 	int resMidBitsLength;
// 	int residualMidBits;
// } LossyCompressionElement;

// typedef struct DoubleValueCompressElement
// {
// 	double data;
// 	long curValue;
// 	unsigned char curBytes[8]; //big_endian
// 	int reqBytesLength;
// 	int resiBitsLength;
// } DoubleValueCompressElement;

// typedef struct DataCompressLogTransform
// {
// 	double data;
// 	long curValue;
// 	unsigned char curBytes[8]; //big_endian
// 	int reqBytesLength;
// 	int resiBitsLength;
// } DataCompressLogTransform;

// typedef struct TightDataPointStorageD
// {
// 	size_t dataSeriesLength;
// 	int allSameData;
// 	double realPrecision;
// 	double medianValue;
// 	char reqLength;	
// 	char radExpo; //used to compute reqLength based on segmented precisions in "pw_rel_compression"

// 	double minLogValue;

// 	int stateNum;
// 	int allNodes;

// 	size_t exactDataNum;
// 	double reservedValue;
	
// 	unsigned char* rtypeArray;
// 	size_t rtypeArray_size;
	
// 	unsigned char* typeArray; //its size is dataSeriesLength/4 (or xxx/4+1) 
// 	size_t typeArray_size;
	
// 	unsigned char* leadNumArray; //its size is exactDataNum/4 (or exactDataNum/4+1)
// 	size_t leadNumArray_size;
	
// 	unsigned char* exactMidBytes;
// 	size_t exactMidBytes_size;
	
// 	unsigned char* residualMidBits;
// 	size_t residualMidBits_size;
	
// 	unsigned int intervals;
	
// 	unsigned char isLossless; //a mark to denote whether it's lossless compression (1 is yes, 0 is no)
	
// 	size_t segment_size;
	
// 	unsigned char* pwrErrBoundBytes;
// 	int pwrErrBoundBytes_size;
		
// 	unsigned char* raBytes;
// 	size_t raBytes_size;
	
// 	unsigned char plus_bits;
// 	unsigned char max_bits;
	
// } TightDataPointStorageD;

// size_t convertIntArray2ByteArray_fast_dynamic(unsigned char* timeStepType, unsigned char resiBitLength, size_t nbEle, unsigned char **bytes)
// {
// 	size_t i = 0, j = 0, k = 0; 
// 	int value;
// 	DynamicByteArray* dba;
// 	new_DBA(&dba, 1024);
// 	int tmp = 0, leftMovSteps = 0;
// 	for(j = 0;j<nbEle;j++)
// 	{
// 		if(resiBitLength==0)
// 			continue;
// 		value = timeStepType[i];
// 		leftMovSteps = getLeftMovingSteps(k, resiBitLength);
// 		if(leftMovSteps < 0)
// 		{
// 			tmp = tmp | (value >> (-leftMovSteps));
// 			addDBA_Data(dba, (unsigned char)tmp);
// 			tmp = 0 | (value << (8+leftMovSteps));
// 		}
// 		else if(leftMovSteps > 0)
// 		{
// 			tmp = tmp | (value << leftMovSteps);
// 		}
// 		else //==0
// 		{
// 			tmp = tmp | value;
// 			addDBA_Data(dba, (unsigned char)tmp);
// 			tmp = 0;
// 		}
// 		i++;
// 		k += resiBitLength;
// 	}
// 	if(leftMovSteps != 0)
// 		addDBA_Data(dba, (unsigned char)tmp);
// 	convertDBAtoBytes(dba, bytes);
// 	size_t size = dba->size;
// 	free_DBA(dba);
// 	return size;
// }

// void new_TightDataPointStorageD(TightDataPointStorageD **this, 
// 		size_t dataSeriesLength, size_t exactDataNum, 
// 		int* type, unsigned char* exactMidBytes, size_t exactMidBytes_size,
// 		unsigned char* leadNumIntArray,  //leadNumIntArray contains readable numbers....
// 		unsigned char* resiMidBits, size_t resiMidBits_size,
// 		unsigned char resiBitLength, 
// 		double realPrecision, double medianValue, char reqLength, unsigned int intervals,
// 		unsigned char* pwrErrBoundBytes, size_t pwrErrBoundBytes_size, unsigned char radExpo) {
// 	//int i = 0;
// 	*this = (TightDataPointStorageD *)malloc(sizeof(TightDataPointStorageD));
// 	(*this)->allSameData = 0;
// 	(*this)->realPrecision = realPrecision;
// 	(*this)->medianValue = medianValue;
// 	(*this)->reqLength = reqLength;

// 	(*this)->dataSeriesLength = dataSeriesLength;
// 	(*this)->exactDataNum = exactDataNum;

// 	(*this)->rtypeArray = NULL;
// 	(*this)->rtypeArray_size = 0;

// 	int stateNum = 2*intervals;

// 	&(*this)->typeArray = SZ_fast_compress_args(2, SZ_DOUBLE, (void *)type, &(*this)->typeArray_size, ABS, 0, 0, 0, 0, 0, 0, dataSeriesLength);

		
// 	(*this)->exactMidBytes = exactMidBytes;
// 	(*this)->exactMidBytes_size = exactMidBytes_size;

// 	(*this)->leadNumArray_size = convertIntArray2ByteArray_fast_2b(leadNumIntArray, exactDataNum, &((*this)->leadNumArray));

// 	(*this)->residualMidBits_size = convertIntArray2ByteArray_fast_dynamic(resiMidBits, resiBitLength, exactDataNum, &((*this)->residualMidBits));
	
// 	(*this)->intervals = intervals;
	
// 	(*this)->isLossless = 0;
	
// 	if(confparams_cpr->errorBoundMode>=PW_REL)
// 		(*this)->pwrErrBoundBytes = pwrErrBoundBytes;
// 	else
// 		(*this)->pwrErrBoundBytes = NULL;
		
// 	(*this)->radExpo = radExpo;
	
// 	(*this)->pwrErrBoundBytes_size = pwrErrBoundBytes_size;
// }

inline void longToBytes_bigEndian(unsigned char *b, unsigned long num) 
{
	b[0] = (unsigned char)(num>>56);
	b[1] = (unsigned char)(num>>48);
	b[2] = (unsigned char)(num>>40);
	b[3] = (unsigned char)(num>>32);
	b[4] = (unsigned char)(num>>24);
	b[5] = (unsigned char)(num>>16);
	b[6] = (unsigned char)(num>>8);
	b[7] = (unsigned char)(num);
//	if(dataEndianType==LITTLE_ENDIAN_DATA)
//		symTransform_8bytes(*b);
}

inline void listAdd_double(double last3CmprsData[3], double value)
{
	last3CmprsData[2] = last3CmprsData[1];
	last3CmprsData[1] = last3CmprsData[0];
	last3CmprsData[0] = value;
}

int compIdenticalLeadingBytesCount_double(unsigned char* preBytes, unsigned char* curBytes)
{
	int i, n = 0;
	for(i=0;i<8;i++)
		if(preBytes[i]==curBytes[i])
			n++;
		else
			break;
	if(n>3) n = 3;
	return n;
}

// size_t dataLength
				// size_t exactDataNum
				// int* type
				// unsigned char* exactMidByteArray->array
				// size_t exactMidByteArray->size
				// unsigned char* exactLeadNumArray->array
				// unsigned char* resiBitArray->array
				// size_t resiBitArray->size,
				// int resiBitsLength,
				// double realPrecision,
				// double medianValue_1
				// (char)reqLength
				// unsigned int quantization_intervals

// unsigned char* get_bytes_for_compress( size_t dataLength,
// 				size_t exactDataNum,
// 				int* type,
// 				unsigned char* exactMidByteArrayarray,
// 				size_t exactMidByteArraysize,
// 				unsigned char* exactLeadNumArrayarray,
// 				unsigned char* resiBitArrayarray,
// 				size_t resiBitArraysize,
// 				int resiBitsLength,
// 				double realPrecision,
// 				double medianValue_1,
// 				(char)reqLength,
// 				unsigned int quantization_intervals){
// 	unsigned char* retBytes = (unsigned char*)malloc(sizeof(double)*(10)+sizeof(unsigned char)*(exactMidByteArraysize+resiBitArraysize));
// }

// inline void addExactData(DynamicByteArray *exactMidByteArray, DynamicIntArray *exactLeadNumArray,
// 		DynamicIntArray *resiBitArray, LossyCompressionElement *lce)
// {
// 	int i;
// 	int leadByteLength = lce->leadingZeroBytes;
// 	addDIA_Data(exactLeadNumArray, leadByteLength);
// 	unsigned char* intMidBytes = lce->integerMidBytes;
// 	int integerMidBytesLength = lce->integerMidBytes_Length;
// 	int resMidBitsLength = lce->resMidBitsLength;
// 	if(intMidBytes!=NULL||resMidBitsLength!=0)
// 	{
// 		if(intMidBytes!=NULL)
// 			for(i = 0;i<integerMidBytesLength;i++)
// 				addDBA_Data(exactMidByteArray, intMidBytes[i]);
// 		if(resMidBitsLength!=0)
// 			addDIA_Data(resiBitArray, lce->residualMidBits);
// 	}
// }

// void updateLossyCompElement_Double(unsigned char* curBytes, unsigned char* preBytes, 
// 		int reqBytesLength, int resiBitsLength,  LossyCompressionElement *lce)
// {
// 	int resiIndex, intMidBytes_Length = 0;
// 	int leadingNum = compIdenticalLeadingBytesCount_double(preBytes, curBytes); //in fact, float is enough for both single-precision and double-precisiond ata.
// 	int fromByteIndex = leadingNum;
// 	int toByteIndex = reqBytesLength; //later on: should use "< toByteIndex" to tarverse....
// 	if(fromByteIndex < toByteIndex)
// 	{
// 		intMidBytes_Length = reqBytesLength - leadingNum;
// 		memcpy(lce->integerMidBytes, &(curBytes[fromByteIndex]), intMidBytes_Length);
// 	}
// 	int resiBits = 0;
// 	if(resiBitsLength!=0)
// 	{
// 		resiIndex = reqBytesLength;
// 		if(resiIndex < 8)
// 			resiBits = (curBytes[resiIndex] & 0xFF) >> (8-resiBitsLength);
// 	}
// 	lce->leadingZeroBytes = leadingNum;
// 	lce->integerMidBytes_Length = intMidBytes_Length;
// 	lce->resMidBitsLength = resiBitsLength;
// 	lce->residualMidBits = resiBits;
// }

// void compressSingleDoubleValue_MSST19(DoubleValueCompressElement *vce, double tgtValue, double precision, int reqLength, int reqBytesLength, int resiBitsLength)
// {
//     ldouble lfBuf;
//     lfBuf.value = tgtValue;

//     int ignBytesLength = 64 - reqLength;
//     if(ignBytesLength<0)
//         ignBytesLength = 0;

//     long tmp_long = lfBuf.lvalue;
//     longToBytes_bigEndian(vce->curBytes, tmp_long);

//     lfBuf.lvalue = (lfBuf.lvalue >> ignBytesLength) << ignBytesLength;

//     //float tmpValue = lfBuf.value;

//     vce->data = lfBuf.value;
//     vce->curValue = tmp_long;
//     vce->reqBytesLength = reqBytesLength;
//     vce->resiBitsLength = resiBitsLength;
// }

inline short getPrecisionReqLength_double(double precision)
{
	ldouble lbuf;
	lbuf.value = precision;
	long lvalue = lbuf.lvalue;
	
	int expValue = (int)((lvalue & 0x7FF0000000000000) >> 52);
	expValue -= 1023;
//	unsigned char the1stManBit = (unsigned char)((lvalue & 0x0008000000000000) >> 51);
//	if(the1stManBit==1)
//		expValue--;
	return (short)expValue;
}

inline void intToBytes_bigEndian(unsigned char *b, unsigned int num)
{
	b[0] = (unsigned char)(num >> 24);	
	b[1] = (unsigned char)(num >> 16);	
	b[2] = (unsigned char)(num >> 8);	
	b[3] = (unsigned char)(num);	
	
	//note: num >> xxx already considered endian_type...
//if(dataEndianType==LITTLE_ENDIAN_DATA)
//		symTransform_4bytes(*b); //change to BIG_ENDIAN_DATA
}

inline short computeReqLength_double_MSST19(double realPrecision)
{
	short reqExpo = getPrecisionReqLength_double(realPrecision);
	return 12-reqExpo;
}

unsigned int roundUpToPowerOf2(unsigned int base)
{
  base -= 1;

  base = base | (base >> 1);
  base = base | (base >> 2);
  base = base | (base >> 4);
  base = base | (base >> 8);
  base = base | (base >> 16);

  return base + 1;
}

unsigned int optimize_intervals_double_1D_opt_MSST19(double *oriData, size_t dataLength, double realPrecision)
{
	size_t i = 0, radiusIndex;
	double pred_value = 0;
	double pred_err;
	size_t *intervals = (size_t*)malloc(MAXRANGERADIUS*sizeof(size_t));
	memset(intervals, 0, MAXRANGERADIUS*sizeof(size_t));
	size_t totalSampleSize = 0;//dataLength/confparams_cpr->sampleDistance;

	double * data_pos = oriData + 2;
	double divider = log2(1+realPrecision)*2;
	int tempIndex = 0;
	while(data_pos - oriData < dataLength){
		if(*data_pos == 0){
        		data_pos += confparams_cpr->sampleDistance;
        		continue;
		}
		tempIndex++;
		totalSampleSize++;
		pred_value = data_pos[-1];
		pred_err = fabs((double)*data_pos / pred_value);
		radiusIndex = (unsigned long)fabs(log2(pred_err)/divider+0.5);
		if(radiusIndex>=MAXRANGERADIUS)
			radiusIndex = MAXRANGERADIUS - 1;
		intervals[radiusIndex]++;

		data_pos += confparams_cpr->sampleDistance;
	}
	//compute the appropriate number
	size_t targetCount = totalSampleSize*confparams_cpr->predThreshold;
	size_t sum = 0;
	for(i=0;i<MAXRANGERADIUS;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=MAXRANGERADIUS)
		i = MAXRANGERADIUS-1;

	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<64)
		powerOf2 = 64;

	free(intervals);
	return powerOf2;
}

double computeRangeSize_double_MSST19(double* oriData, size_t size, double* valueRangeSize, double* medianValue, unsigned char * signs, bool* positive, double* nearZero)
{
    size_t i = 0;
    double min = oriData[0];
    double max = min;
    *nearZero = min;

    for(i=1;i<size;i++)
    {
        double data = oriData[i];
        if(data <0){
            signs[i] = 1;
            *positive = false;
        }
        if(oriData[i] != 0 && fabs(oriData[i]) < fabs(*nearZero)){
            *nearZero = oriData[i];
        }
        if(min>data)
            min = data;
        else if(max<data)
            max = data;
    }

    *valueRangeSize = max - min;
    *medianValue = min + *valueRangeSize/2;
    return min;
}

void cost_start()
{
	totalCost = 0;
	gettimeofday(&costStart, NULL);
}

void cost_end()
{
	double elapsed;
	struct timeval costEnd;
	gettimeofday(&costEnd, NULL);
	elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
	totalCost += elapsed;
}


void usage()
{
	printf("Usage: szx <options>\n");
	printf("Options:\n");
	printf("* operation type:\n");
	printf("	-z <compressed file>: the compression operation with an optionally specified output file.\n");
	printf("                          (the compressed file will be named as <input_file>.sz if not specified)\n");
	printf("	-x <decompressed file>: the decompression operation with an optionally specified output file\n");
	printf("                      (the decompressed file will be named as <cmpred_file>.out if not specified)\n");
	printf("	-h: print the help information\n");
	printf("	-v: print the version number\n");	
	printf("* data type:\n");
	printf("	-f: single precision (float type)\n");
	printf("	-d: double precision (double type)\n");
	printf("* execution mode:\n");
	printf("	-m <mode>:\n");
	printf("		-m 1: nonblock+serial\n");
	printf("		-m 2: blocked+serial\n");
	printf("		-m 3: blocked+openmp\n");
	printf("		-m 4: blocked+randomaccess+serial\n");
	printf("* error control: (the error control parameters here will overwrite the setting in sz.config)\n");
	printf("	-M <error bound mode> : 2 options as follows. \n");
	printf("		ABS (absolute error bound)\n");
	printf("		REL (value range based error bound, so a.k.a., VR_REL)\n");
	printf("	-A <absolute error bound>: specifying absolute error bound\n");
	printf("	-R <value_range based relative error bound>: specifying relative error bound\n");
	printf("* input data file:\n");
	printf("	-i <original data file> : original data file\n");
	printf("	-s <compressed data file> : compressed data file in decompression\n");
	printf("* output type of decompressed file: \n");
	printf("	-b (by default) : decompressed file stored in binary format\n");
	printf("	-t : decompreadded file stored in text format\n");
	printf("* dimensions: \n");
	printf("	-1 <nx> : dimension for 1D data such as data[nx]\n");
	printf("	-2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]\n");
	printf("	-3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] \n");
	printf("	-4 <nx> <ny> <nz> <np>: dimensions for 4D data such as data[np][nz][ny][nx] \n");
	printf("* print compression results: \n");
	printf("	-a : print compression results such as distortions\n");
	printf("* examples: \n");
	printf("	sz -z -f -i testdata/x86/testfloat_8_8_128.dat -3 8 8 128 -M ABS -A 1E-3\n");
	printf("	sz -x -f -s testdata/x86/testfloat_8_8_128.dat.sz -3 8 8 128 -a\n");
	exit(0);
}


int main(int argc, char* argv[])
{
	int binaryOutput = 1;
	int printCmpResults = 0;
	int isCompression = -1000; //1 : compression ; 0: decompression
	int dataType = 0; //0: single precision ; 1: double precision
	char* inPath = NULL;
	char* cmpPath = NULL;
	char* decPath = NULL;
	
	char* errBoundMode = NULL;
	char* absErrBound = NULL;
	char* relErrBound = NULL;
	float absErrorBound = 0, relBoundRatio = 0;

	int doPredQuant = 0;

	int fastMode = SZx_WITH_BLOCK_FAST_CMPR; //1: non-blocked+serial, 2: blocked+serial, 3: blocked+openmp, 4: blocked+randomaccess+serial
	size_t r5 = 0;
	size_t r4 = 0;
	size_t r3 = 0;
	size_t r2 = 0; 
	size_t r1 = 0;
	
	size_t i = 0;
	int status;
	size_t nbEle;
	if(argc==1)
		usage();
	
	for(i=1;i<argc;i++)
	{
		if (argv[i][0] != '-' || argv[i][2])
			usage();
		switch (argv[i][1])
		{
		case 'Q':
			doPredQuant = 1;
			break;
		case 'h':
			usage();
			exit(0);
		case 'v':
			printf("version: %d.%d.%d\n", SZx_VER_MAJOR, SZx_VER_MINOR, SZx_VER_BUILD);
			exit(0);
		case 'b': 
			binaryOutput = 1;
			break;
		case 't': 
			binaryOutput = 0;
			break;
		case 'a':
			printCmpResults = 1;
			break;
		case 'z':
			isCompression = 1;
			if (i+1 < argc)
			{
				cmpPath = argv[i+1];
				if(cmpPath[0]!='-')
					i++;
				else
					cmpPath = NULL;
			}
			break;
		case 'x': 
			isCompression = 0;
			if (i+1 < argc)
			{
				decPath = argv[i+1];
				if(decPath[0]!='-')
					i++;
				else
					decPath = NULL;
			}			
			break;
		case 'f': 
			dataType = 0;
			break;
		case 'd':
			dataType = 1;
			break;
		case 'i':
			if (++i == argc)
				usage();
			inPath = argv[i];		
			break;
		case 's':
			if (++i == argc)
				usage();
			cmpPath = argv[i];
			break;
		case '1': 
			if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1)
				usage();
			break;
		case '2':
			if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
				++i == argc || sscanf(argv[i], "%zu", &r2) != 1)
				usage();
			break;
		case '3':
			if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
				++i == argc || sscanf(argv[i], "%zu", &r2) != 1 ||
				++i == argc || sscanf(argv[i], "%zu", &r3) != 1)
				usage();		
			break;
		case '4':
			if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
				++i == argc || sscanf(argv[i], "%zu", &r2) != 1 ||
				++i == argc || sscanf(argv[i], "%zu", &r3) != 1 ||
				++i == argc || sscanf(argv[i], "%zu", &r4) != 1)
				usage();		
			break;
		case 'm':
			if (++i == argc)
				usage();
			fastMode = atoi(argv[i]);
			break;
		case 'M':
			if (++i == argc)
				usage();
			errBoundMode = argv[i];
			break;
		case 'A':
			if (++i == argc)
				usage();
			absErrBound = argv[i];
			break;
		case 'R':
			if (++i == argc)
				usage();
			relErrBound = argv[i];
			break;
		default: 
			usage();
			break;
		}
	}

	if((inPath==NULL) & (cmpPath == NULL))
	{
		printf("Error: you need to specify either a raw binary data file or a compressed data file as input\n");
		usage();
		exit(0);
	}

	int errorBoundMode = 0;
	if(isCompression == 1 && errBoundMode != NULL)
	{
		if(strcmp(errBoundMode, "ABS")==0)
			errorBoundMode = ABS;
		else if(strcmp(errBoundMode, "REL")==0||strcmp(errBoundMode, "VR_REL")==0)
			errorBoundMode = REL;
		else
		{
			printf("Error: wrong error bound mode setting by using the option '-M'\n");
			usage();
			exit(0);
		}
	}
	
	char outputFilePath[256];	
	unsigned char *bytes = NULL; //the binary data read from "compressed data file"
	size_t byteLength = 0; 
	if(isCompression == 1)
	{
		if(absErrBound != NULL)
			absErrorBound = atof(absErrBound);
		
		if(relErrBound != NULL)
			relBoundRatio = atof(relErrBound);
		
		size_t outSize;	
		if(dataType == SZ_FLOAT) //single precision
		{
			float *data = readFloatData(inPath, &nbEle, &status);
			
			//exclude the value range calculation time from the compression to get more accurate evaluation for performance
			if(errorBoundMode == REL)
			{
				float max = data[0], min = data[0];
				for (size_t i = 0; i < nbEle; i++) {
					if (data[i] > max) {
						max = data[i];
					}
					if (data[i] < min) {
						min = data[i];
					}
				}
				absErrorBound = (max - min) * relBoundRatio;
				errorBoundMode = ABS;
			}
			
			if(status!=SZ_SCES)
			{
				printf("Error: cannot read the input file: %s\n", inPath);
				exit(0);
			}
			cost_start();
			bytes = SZ_fast_compress_args(fastMode, SZ_FLOAT, data, &outSize, errorBoundMode, absErrorBound, relBoundRatio, r5, r4, r3, r2, r1);
			cost_end();
			if(cmpPath == NULL)
				sprintf(outputFilePath, "%s.szx", inPath);
			else
				strcpy(outputFilePath, cmpPath);
			writeByteData(bytes, outSize, outputFilePath, &status);		
			free(data);
			if(status != SZ_SCES)
			{
				printf("Error: data file %s cannot be written!\n", outputFilePath);
				exit(0);
			}
			printf("compression time = %f\n", totalCost);
			printf("compressed data file: %s\n", outputFilePath);			
		}
		else //dataType == 1: double precision
		{
			double *data = readDoubleData(inPath, &nbEle, &status);	
			
			//exclude the value range calculation time from the compression to get more accurate evaluation for performance
			if(errorBoundMode == REL)
			{
				double max = data[0], min = data[0];
				for (size_t i = 0; i < nbEle; i++) {
					if (data[i] > max) {
						max = data[i];
					}
					if (data[i] < min) {
						min = data[i];
					}
				}
				absErrorBound = (max - min) * relBoundRatio;
				errorBoundMode = ABS;
			}			
			
			if(status!=SZ_SCES)
			{
				printf("Error: cannot read the input file: %s\n", inPath);
				exit(0);
			}
			cost_start();


			if(doPredQuant){

				unsigned char *bytesInt = NULL;
				unsigned char *bytesUnpred = NULL;
				size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);

				double *oriData = (double *)data;
				double prediction = oriData[0];

				

				int plus_bits = 3;
				// double* oriData = (double *) data;
				double valueRangeSize = 0, medianValue = 0;

				double pwrErrRatio = (double) absErrorBound;

				unsigned char * signs = NULL;
				bool positive = true;
				double nearZero = 0.0;
				double min = 0;

				signs = (unsigned char *) malloc(dataLength);
				memset(signs, 0, dataLength);
				min = computeRangeSize_double_MSST19(oriData, dataLength, &valueRangeSize, &medianValue, signs, &positive, &nearZero);
			
				double max = min+valueRangeSize;
				// confparams_cpr->dmin = min;
				// confparams_cpr->dmax = max;


				double realPrecision = 0;
				// realPrecision = getRealPrecision_double(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio, &status);
				size_t outSize = 0;
				unsigned char* newByteData;

				double multiplier = pow((1+pwrErrRatio), -3.0001);
				for(int i=0; i<dataLength; i++){
					if(oriData[i] == 0){
						oriData[i] = nearZero * multiplier;
					}
				}

				double median_log = sqrt(fabs(nearZero * max));

				unsigned int quantization_intervals;
				quantization_intervals = optimize_intervals_double_1D_opt_MSST19(oriData, dataLength, realPrecision);

				int intvRadius = quantization_intervals/2;

				double* precisionTable = (double*)malloc(sizeof(double) * quantization_intervals);
				double inv = 2.0-pow(2, -(plus_bits));
				for(int i=0; i<quantization_intervals; i++){
					double test = pow((1+realPrecision), inv*(i - intvRadius));
					precisionTable[i] = test;
				}
				struct TopLevelTableWideInterval levelTable;
   				MultiLevelCacheTableWideIntervalBuild(&levelTable, precisionTable, quantization_intervals, realPrecision, plus_bits);

				int* type = (int*) malloc(dataLength*sizeof(int));
				int state;
				//double checkRadius;
				double curData;
				double pred = prediction;

				double* spaceFillingValue = oriData;

				double predRelErrRatio;

				const uint64_t top = levelTable.topIndex, base = levelTable.baseIndex;
				const uint64_t range = top - base;
				const int bits = levelTable.bits;
				uint64_t* const buffer = (uint64_t*)&predRelErrRatio;
				const int shift = 52-bits;
				uint64_t expoIndex, mantiIndex;
				uint16_t* tables[range+1];
				for(int i=0; i<=range; i++){
					tables[i] = levelTable.subTables[i].table;
				}

				size_t unpred_count = 1;
				type[0] = 0;

				for(i=1;i<dataLength;i++)
				{
					curData = spaceFillingValue[i];
					predRelErrRatio = curData / pred;

					expoIndex = ((*buffer & 0x7fffffffffffffff) >> 52) - base;
					if(expoIndex <= range){
						mantiIndex = (*buffer & 0x000fffffffffffff) >> shift;
						state = tables[expoIndex][mantiIndex];
					}else{
						state = 0;
					}

					if(state)
					{
						type[i] = state;
						pred *= precisionTable[state];
						//hit++;
						continue;
					}

					// //unpredictable data processing
					type[i] = 0;
					// compressSingleDoubleValue_MSST19(vce, curData, realPrecision, reqLength, reqBytesLength, resiBitsLength);
					// updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					// memcpy(preDataBytes,vce->curBytes,8);
					// addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					// pred =  vce->data;
					unpred_count++;

				}

				double *unpredData = (double *)malloc(sizeof(double)*unpred_count);

				int unpred_ind = 0;
				for (size_t i = 0; i < dataLength; i++)
				{
					if (type[i] == 0)
					{
						unpredData[unpred_ind] = spaceFillingValue[i];
						unpred_ind++;
					}
					
				}
				
				size_t sizeUnpred, sizeType;

				bytesInt = SZ_fast_compress_args(fastMode, SZ_DOUBLE, (void *)unpredData, &sizeUnpred, ABS, 0.0, 0.0, 0, 0, 0, 0, unpred_count);
				bytesUnpred = SZ_fast_compress_args(fastMode, SZ_FLOAT, (void *)type, &sizeType, ABS, 0.0, 0.0, 0, 0, 0, 0, dataLength);

				cost_end();
				if(cmpPath == NULL)
					sprintf(outputFilePath, "%s-unpred.szx", inPath);
				else
					strcpy(outputFilePath, cmpPath);
				writeByteData(bytesUnpred, sizeType, outputFilePath, &status);		
				free(data);
				if(status != SZ_SCES)
				{
					printf("Error: data file %s cannot be written!\n", outputFilePath);
					exit(0);
				}

				if(cmpPath == NULL)
					sprintf(outputFilePath, "%s-quants.szx", inPath);
				else
					strcpy(outputFilePath, cmpPath);
				writeByteData(bytesUnpred, sizeUnpred, outputFilePath, &status);		
				free(data);
				if(status != SZ_SCES)
				{
					printf("Error: data file %s cannot be written!\n", outputFilePath);
					exit(0);
				}

				printf("compression time = %f\n", totalCost);
				printf("compressed data file: %s\n", outputFilePath);

				freeTopLevelTableWideInterval(&levelTable);
			}else{

				bytes = SZ_fast_compress_args(fastMode, SZ_DOUBLE, data, &outSize, errorBoundMode, absErrorBound, relBoundRatio, r5, r4, r3, r2, r1);
				cost_end();
				if(cmpPath == NULL)
					sprintf(outputFilePath, "%s.szx", inPath);
				else
					strcpy(outputFilePath, cmpPath);
				writeByteData(bytes, outSize, outputFilePath, &status);		
				free(data);
				if(status != SZ_SCES)
				{
					printf("Error: data file %s cannot be written!\n", outputFilePath);
					exit(0);
				}		
				printf("compression time = %f\n", totalCost);
				printf("compressed data file: %s\n", outputFilePath);
			}
		}

		if (printCmpResults == 1)
		{
			printf ("Error: -a can be only used in decompression.\n");
		}
	}
	else if(isCompression == 0) //decompression
	{
		if(printCmpResults)
		{
			if(inPath==NULL)
			{
				printf("Error: Since you add -a option (analysis), please specify the original data path by -i <path>.\n");
				exit(0);
			}
		}		
		
		char outputFilePath[256];
		
		if(r2==0)
			nbEle = r1;
		else if(r3==0)
			nbEle = r1*r2;
		else if(r4==0)
			nbEle = r1*r2*r3;
		else if(r5==0)
			nbEle = r1*r2*r3*r4;
		else
			nbEle = r1*r2*r3*r4*r5;

		if(checkFileExistance(cmpPath)==0)
		{
			printf("Error: compression file (%s) is not readable.\n", cmpPath);
			exit(0);
		}

		if(dataType == SZ_FLOAT)
		{
			bytes = readByteData(cmpPath, &byteLength, &status);
			if(status!=SZ_SCES)
			{
				printf("Error: %s cannot be read!\n", cmpPath);
				exit(0);
			}
			cost_start();
			float *data = SZ_fast_decompress(fastMode, SZ_FLOAT, bytes, byteLength, r5, r4, r3, r2, r1);
			cost_end();
			if(decPath == NULL)
				sprintf(outputFilePath, "%s.out", cmpPath);	
			else
				strcpy(outputFilePath, decPath);
			if(binaryOutput==1)		
				writeFloatData_inBytes(data, nbEle, outputFilePath, &status);
			else //txt output
				writeFloatData(data, nbEle, outputFilePath, &status);

			if(status!=SZ_SCES)
			{
				printf("Error: %s cannot be written!\n", outputFilePath);
				exit(0);
			}
			
			if(printCmpResults)
			{
				if(inPath==NULL)
				{
					printf("Error: Since you add -a option (analysis), please specify the original data path by -i <path>.\n");
					exit(0);
				}
				//compute the distortion / compression errors...
				size_t totalNbEle;
				float *ori_data = readFloatData(inPath, &totalNbEle, &status);
				if(status!=SZ_SCES)
				{
					printf("Error: %s cannot be read!\n", inPath);
					exit(0);
				}

				size_t i = 0;
				float Max = 0, Min = 0, diffMax = 0;
				Max = ori_data[0];
				Min = ori_data[0];
				diffMax = fabs(data[0] - ori_data[0]);
				double sum1 = 0, sum2 = 0, sum22 = 0;
				for (i = 0; i < nbEle; i++)
				{
					sum1 += ori_data[i];
					sum2 += data[i];
					sum22 += data[i]*data[i];
				}
				double mean1 = sum1/nbEle;
				double mean2 = sum2/nbEle;

				double sum3 = 0, sum4 = 0;
				double sum = 0, prodSum = 0, relerr = 0;

				double maxpw_relerr = 0; 
				for (i = 0; i < nbEle; i++)
				{
					if (Max < ori_data[i]) Max = ori_data[i];
					if (Min > ori_data[i]) Min = ori_data[i];
					
					float err = fabs(data[i] - ori_data[i]);

					if(ori_data[i]!=0)
					{
						relerr = err/fabs(ori_data[i]);
						if(maxpw_relerr<relerr)
							maxpw_relerr = relerr;
					}

					if (diffMax < err)
						diffMax = err;
					prodSum += (ori_data[i]-mean1)*(data[i]-mean2);
					sum3 += (ori_data[i] - mean1)*(ori_data[i]-mean1);
					sum4 += (data[i] - mean2)*(data[i]-mean2);
					sum += err*err;	
				}
				double std1 = sqrt(sum3/nbEle);
				double std2 = sqrt(sum4/nbEle);
				double ee = prodSum/nbEle;
				double acEff = ee/std1/std2;

				double mse = sum/nbEle;
				double range = Max - Min;
				double psnr = 20*log10(range)-10*log10(mse);
				double nrmse = sqrt(mse)/range;
				double compressionRatio = 1.0*nbEle*sizeof(float)/byteLength;
				double normErr = sqrt(sum);
				double normErr_norm = normErr/sqrt(sum22);		
		
				printf ("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
				printf ("Max absolute error = %.10f\n", diffMax);
				printf ("Max relative error = %f\n", diffMax/(Max-Min));
				printf ("Max pw relative error = %f\n", maxpw_relerr);
				printf ("PSNR = %f, NRMSE= %.20G\n", psnr,nrmse);
				printf ("normError = %f, normErr_norm = %f\n", normErr, normErr_norm);
				printf ("acEff=%f\n", acEff);	
				printf ("compressionRatio=%f\n", compressionRatio);
				
				free(ori_data);
			}
			free(data);	
			
			printf("decompression time = %f seconds.\n", totalCost);
			printf("decompressed data file: %s\n", outputFilePath);							
		}
		else //double-data
		{
			bytes = readByteData(cmpPath, &byteLength, &status);
			if(status!=SZ_SCES)
			{
				printf("Error: %s cannot be read!\n", cmpPath);
				exit(0);
			}
			cost_start();
			double* data = SZ_fast_decompress(fastMode, SZ_DOUBLE, bytes, byteLength, r5, r4, r3, r2, r1);
			cost_end();
			if(decPath == NULL)
				sprintf(outputFilePath, "%s.out", cmpPath);	
			else
				strcpy(outputFilePath, decPath);
			if(binaryOutput==1)		
				writeDoubleData_inBytes(data, nbEle, outputFilePath, &status);
			else //txt output
				writeDoubleData(data, nbEle, outputFilePath, &status);			
			if(status!=SZ_SCES)
			{
				printf("Error: %s cannot be written!\n", outputFilePath);
				exit(0);
			}
						
			printf("decompression time = %f seconds.\n", totalCost);
			printf("decompressed data file: %s\n", outputFilePath);
			
			if(printCmpResults)
			{
				if(inPath==NULL)
				{
					printf("Error: Since you add -a option (analysis), please specify the original data path by -i <path>.\n");
					exit(0);
				}
				size_t totalNbEle;

				//compute the distortion / compression errors...
				double *ori_data = readDoubleData(inPath, &totalNbEle, &status);
				if(status!=SZ_SCES)
				{
					printf("Error: %s cannot be read!\n", inPath);
					exit(0);
				}

				size_t i = 0;
				double Max = 0, Min = 0, diffMax = 0;
				Max = ori_data[0];
				Min = ori_data[0];
				diffMax = data[0]>ori_data[0]?data[0]-ori_data[0]:ori_data[0]-data[0];

				//diffMax = fabs(data[0] - ori_data[0]);
				double sum1 = 0, sum2 = 0, sum22 = 0;

				for (i = 0; i < nbEle; i++)
				{
					sum1 += ori_data[i];
					sum2 += data[i];
					sum22 += data[i]*data[i];
				}
				double mean1 = sum1/nbEle;
				double mean2 = sum2/nbEle;

				double sum3 = 0, sum4 = 0;
				double sum = 0, prodSum = 0, relerr = 0;

				double maxpw_relerr = 0; 
				for (i = 0; i < nbEle; i++)
				{
					if (Max < ori_data[i]) Max = ori_data[i];
					if (Min > ori_data[i]) Min = ori_data[i];

					float err = fabs(data[i] - ori_data[i]);
					if(ori_data[i]!=0)
					{
						relerr = err/fabs(ori_data[i]);
						if(maxpw_relerr<relerr)
						  maxpw_relerr = relerr;
					}

					if (diffMax < err)
					  diffMax = err;
					prodSum += (ori_data[i]-mean1)*(data[i]-mean2);
					sum3 += (ori_data[i] - mean1)*(ori_data[i]-mean1);
					sum4 += (data[i] - mean2)*(data[i]-mean2);
					sum += err*err;	
				}
				double std1 = sqrt(sum3/nbEle);
				double std2 = sqrt(sum4/nbEle);
				double ee = prodSum/nbEle;
				double acEff = ee/std1/std2;

				double mse = sum/nbEle;
				double range = Max - Min;
				double psnr = 20*log10(range)-10*log10(mse);
				double normErr = sqrt(sum);
				double normErr_norm = normErr/sqrt(sum22);
				double nrmse = sqrt(mse)/range;

				double compressionRatio = 1.0*nbEle*sizeof(double)/byteLength;

				printf ("Min = %.20G, Max = %.20G, range = %.20G\n", Min, Max, range);
				printf ("Max absolute error = %.10f\n", diffMax);
				printf ("Max relative error = %f\n", diffMax/(Max-Min));
				printf ("Max pw relative error = %f\n", maxpw_relerr);
				printf ("PSNR = %f, NRMSE = %.20G\n", psnr,nrmse);
				printf ("normErr = %f, normErr_norm = %f\n", normErr, normErr_norm);
				printf ("acEff = %f\n", acEff);
				printf ("compressionRatio = %f\n", compressionRatio);
				
				free(ori_data);
			}			
			free(data);								
		}	
	}
	
	free(bytes);
	
}
