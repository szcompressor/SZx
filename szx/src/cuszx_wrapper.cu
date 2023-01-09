#include "cuszx_entry.h"
#include "szx_defines.h"
#include "szx_BytesToolkit.h"
#include "szx_TypeManager.h"
#include "timingGPU.h"

extern "C"{
    void cuSZx_integrated_compress(unsigned char *bytes, float *data, float r2r_threshold, float r2r_err, size_t nbEle, int blockSize, size_t *outSize){
        float max,min;
        max = data[0];
        min = data[0];
        for (size_t i = 0; i < nbEle; i++)
        {
            if(data[i] > max) max = data[i];
            if(data[i] < min) min = data[i];
        }
        
        float threshold = r2r_threshold*(max-min);
        float errBound = r2r_err*(max-min);
        bytes = cuSZx_fast_compress_args_unpredictable_blocked_float(data, outSize, errBound, nbEle, blockSize, threshold);
   	printf("outSize %p\n", bytes);
    }

    void cuSZx_integrated_decompress(float *data, unsigned char *bytes, size_t nbEle){
        printf("test\n");
	//cuSZx_fast_decompress_args_unpredictable_blocked_float(&data, nbEle, bytes);
    }

    unsigned char* cuSZx_device_compress(float *oriData, size_t *outSize, float absErrBound, size_t nbEle, int blockSize, float threshold){
        return device_ptr_cuSZx_compress_float(oriData, outSize, absErrBound, nbEle, blockSize, threshold);
    }
    
}
