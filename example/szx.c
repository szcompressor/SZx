#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "szx.h"
#include "szx_rw.h"
#ifdef _OPENMP
#include "omp.h"
#endif

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;


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
	printf("	-M <error bound mode> : 3 options as follows. \n");
	printf("		ABS (absolute error bound)\n");
	printf("		REL (value range based error bound, so a.k.a., VR_REL)\n");
	printf("		FXR (fix ratio)\n");
	printf("	-A <absolute error bound>: specifying absolute error bound\n");
	printf("	-R <value_range based relative error bound>: specifying relative error bound\n");
	printf("	-C <compression ratio>: specifying the compression ratio.\n");
	printf("	-T <tolerance of fixed compression ratio>: specifying the acceptable relative tolerance\n");
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
	printf("	sz -z -f -i ~/Data/Hurricane-ISA/CLOUDf48.bin.dat -3 500 500 100 -M FXR -C 10 -T 0.1\n");
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
	char* compRatio = NULL;
	float absErrorBound = 0, relBoundRatio = 0;
	float compressRatio = 0;
	float tolerance = 0;

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
		case 'C':
			if (++i == argc)
				usage();
			compRatio = argv[i];
			break;
		case 'T':
			if (++i == argc)
				usage();
			tolerance = atof(argv[i]);
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
		else if(strcmp(errBoundMode, "FXR")==0)
			errorBoundMode = FXR;
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
	
		if(compRatio != NULL)
			compressRatio = atof(compRatio);

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
			bytes = SZ_fast_compress_args_test(fastMode, SZ_FLOAT, data, &outSize, errorBoundMode, absErrorBound, relBoundRatio, compressRatio, tolerance, r5, r4, r3, r2, r1);
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
			bytes = SZ_fast_compress_args_test(fastMode, SZ_DOUBLE, data, &outSize, errorBoundMode, absErrorBound, relBoundRatio, compressRatio, tolerance, r5, r4, r3, r2, r1);
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
			float *data = SZ_fast_decompress_test(fastMode, SZ_FLOAT, bytes, byteLength, r5, r4, r3, r2, r1);
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
			double* data = SZ_fast_decompress_test(fastMode, SZ_DOUBLE, bytes, byteLength, r5, r4, r3, r2, r1);
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
