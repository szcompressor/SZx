#include "cuszx_entry.h"
#include "szx_defines.h"
#include "szx_BytesToolkit.h"
#include "szx_TypeManager.h"
#include "timingGPU.h"
#include <stdio.h>

#define DO_ARTIFACT_ANALYSIS

TimingGPU timer_GPU;
void bin(unsigned n)
{
    unsigned i;
    for (i = 1 << 31; i > 0; i = i / 2)
        (n & i) ? printf("1") : printf("0");
}

int _post_proc(float *oriData, unsigned char *meta, short *offsets, unsigned char *midBytes, unsigned char *outBytes, size_t nbEle, int blockSize)
{
    int out_size = 0;
    size_t nbConstantBlocks = 0;
    size_t nbBlocks = nbEle/blockSize;
    size_t ncBytes = blockSize/4;
    size_t mSize = sizeof(float)+1+ncBytes; //Number of bytes for each data block's metadata.
    out_size += 5+sizeof(size_t)+sizeof(float)*nbBlocks;
    if (nbBlocks%8==0)
        out_size += nbBlocks/8;
    else
        out_size += nbBlocks/8+1;
    for (int i=0; i<nbBlocks; i++){
        if (meta[i]==0) nbConstantBlocks++;
        else out_size += 1+(blockSize/4)+offsets[i];
    }
    out_size += (nbBlocks-nbConstantBlocks)*sizeof(short)+(nbEle%blockSize)*sizeof(float);

    //outBytes = (unsigned char*)malloc(out_size);
	unsigned char* r = outBytes; 
	r[0] = SZx_VER_MAJOR;
	r[1] = SZx_VER_MINOR;
	r[2] = 1;
	r[3] = 0; // indicates this is not a random access version
	r[4] = (unsigned char)blockSize;
	r=r+5; //1 byte
	sizeToBytes(r, nbConstantBlocks);
	r += sizeof(size_t); 
	r += convertIntArray2ByteArray_fast_1b_args(meta, nbBlocks, r);
    memcpy(r, oriData+nbBlocks*blockSize, (nbEle%blockSize)*sizeof(float));
    r += (nbEle%blockSize)*sizeof(float);
    unsigned char* c = r;
    unsigned char* o = c+nbConstantBlocks*sizeof(float);
    unsigned char* nc = o+(nbBlocks-nbConstantBlocks)*sizeof(short);
    for (int i=0; i<nbBlocks; i++){
        
        if (meta[i]==0){
            memcpy(c, meta+(nbBlocks+i*mSize), sizeof(float));
            c += sizeof(float);
        }else{
            shortToBytes(o, offsets[i]);
            o += sizeof(short);
            memcpy(nc, meta+(nbBlocks+i*mSize), mSize);
            nc += mSize; 
            memcpy(nc, midBytes+(i*blockSize*sizeof(float)), offsets[i]);
            nc += offsets[i];
        } 
    }

    return out_size;
}

__device__ void fits_linear(float y1, float y2, float* block_data, int blockSize, float absErrBound, float median){
    int fits_bound = 1;
    float mse_l = 0.0;
    float mse_m = 0.0;
    for (float x = 0.0, int i=0; i < blockSize; i++, x+=1.0)
    {
        float y = y1+((y2-y1)/((float)blockSize))*x;
        float diff_l = y - block_data[i];
        float diff_m = median - block_data[i];
        if (fabs(diff_l) >absErrBound)
        {
            fits_bound = 0;
        }
        mse_l+= diff_l*diff_l;
        mse_m+= diff_m*diff_m;
    }

    mse_l = mse_l/(float)blockSize;
    mse_m = mse_m/(float)blockSize;
}

__global__ void artifact_analysis(float *oriData, int blockSize, float absErrBound, size_t nbEle, unsigned char *meta, int xdim, int ydim, int zdim, unsigned char *artifact_type){
    size_t ncBytes = blockSize/4;
    size_t mSize = sizeof(float)+1+ncBytes;
    size_t nbBlocks = nbEle/blockSize;

    int blocks_per_z = xdim*ydim/blockSize;

    for (int d = blockIdx.x; d < zdim; d+=gridDim.x) // Each block handles 2D grid
    {
        int startblock = d*blocks_per_z;
        int blocks_per_y = xdim/blockSize;
        for (int r = threadIdx.x; r < ydim; r+=blockDim.x) // Each thread handles one row of 2D grid
        {
            int rowstartblock = startblock + r*blocks_per_y;
            float prevblock_interpolate, nextblock_interpolate, median;
            memcpy(&prevblock_interpolate, meta+(nbBlocks+rowstartblock*mSize), sizeof(float));
            for (int c = rowstartblock; c < rowstartblock+blocks_per_y; c++) // For each chunk in row...
            {
                if (c > rowstartblock && c < (rowstartblock+blocks_per_y)-1)
                {
                    
                    if (meta[c-1] == 0 && meta[c] == 0 && meta[c+1] == 0) //three constant blocks
                    {
                        memcpy(&nextblock_interpolate, meta+(nbBlocks+(c+1)*mSize), sizeof(float));
                        memcpy(&median, meta+(nbBlocks+(c)*mSize), sizeof(float));
                        int fits_bound = 1;
                        float mse_l = 0.0;
                        float mse_m = 0.0;
                        float y1 = prevblock_interpolate;
                        float y2 = nextblock_interpolate;
                        float last_point = 0.0;
                        for (float x = 0.0, int i=0; i < blockSize; i++, x+=1.0)
                        {
                            float y = y1+((y2-y1)/((float)blockSize))*x;
                            float diff_l = y - block_data[i];
                            float diff_m = median - block_data[i];
                            if (fabs(diff_l) >absErrBound)
                            {
                                fits_bound = 0;
                            }
                            if (i == blockSize-1)
                            {
                                last_point = y;
                            }
                            
                            mse_l+= diff_l*diff_l;
                            mse_m+= diff_m*diff_m;
                        }

                        mse_l = mse_l/(float)blockSize;
                        mse_m = mse_m/(float)blockSize;
                        if (mse_l < mse_m && fits_bound==1)
                        {
                            artifact_type[c] = 1;
                        }else{
                            artifact_type[c] = 0;
                        }
                        
                    }
                    
                }
                
            }
            
        }
        
    }
    
}

// Note that median values must be of length nbBlocks not number of nonconstant blocks
__global__ void artifact_reconstruct(float *newData, float *medianvalues, unsigned char* states, int blockSize, float absErrBound, size_t nbEle, unsigned char *meta, int xdim, int ydim, int zdim, unsigned char *artifact_type){
    size_t ncBytes = blockSize/4;
    size_t mSize = sizeof(float)+1+ncBytes;
    size_t nbBlocks = nbEle/blockSize;

    int blocks_per_z = xdim*ydim/blockSize;

    for (int d = blockIdx.x; d < zdim; d+=gridDim.x) // Each block handles 2D grid
    {
        int startblock = d*blocks_per_z;
        int blocks_per_y = xdim/blockSize;
        for (int r = threadIdx.x; r < ydim; r+=blockDim.x) // Each thread handles one row of 2D grid
        {
            int rowstartblock = startblock + r*blocks_per_y;
            float prevblock_interpolate, nextblock_interpolate, median;
            prevblock_interpolate = medianvalues[rowstartblock];
            for (int c = rowstartblock; c < rowstartblock+blocks_per_y; c++) // For each chunk in row...
            {
                if (c > rowstartblock && c < (rowstartblock+blocks_per_y)-1)
                {
                    
                    if (artifact_type[c]==1) //three constant blocks
                    {

                        nextblock_interpolate = medianvalues[c];
                        median = medianvalues[c];
                        int fits_bound = 1;
                        float mse_l = 0.0;
                        float mse_m = 0.0;
                        float y1 = prevblock_interpolate;
                        float y2 = nextblock_interpolate;
                        float last_point = 0.0;
                        for (float x = 0.0, int i=0; i < blockSize; i++, x+=1.0)
                        {
                            float y = y1+((y2-y1)/((float)blockSize))*x;
                            float diff_l = y - block_data[i];
                            float diff_m = median - block_data[i];
                            if (fabs(diff_l) >absErrBound)
                            {
                                fits_bound = 0;
                            }
                            if (i == blockSize-1)
                            {
                                last_point = y;
                            }
                            
                            mse_l+= diff_l*diff_l;
                            mse_m+= diff_m*diff_m;
                        }

                        mse_l = mse_l/(float)blockSize;
                        mse_m = mse_m/(float)blockSize;
                        if (mse_l < mse_m && fits_bound==1)
                        {
                            artifact_type[c] = 1;
                        }else{
                            artifact_type[c] = 0;
                        }
                        
                    }
                    
                }
                
            }
            
        }
        
    }
    
}


unsigned char* cuSZx_fast_compress_args_unpredictable_blocked_float(float *oriData, size_t *outSize, float absErrBound, size_t nbEle, int blockSize)
{

	float* d_oriData;
    cudaMalloc((void**)&d_oriData, sizeof(float)*nbEle); 
    cudaMemcpy(d_oriData, oriData, sizeof(float)*nbEle, cudaMemcpyHostToDevice); 

	size_t nbBlocks = nbEle/blockSize;
	size_t remainCount = nbEle%blockSize;
	size_t actualNBBlocks = remainCount==0 ? nbBlocks : nbBlocks+1;

    size_t ncBytes = blockSize/4;
    //ncBytes = (blockSize+1)%4==0 ? ncBytes : ncBytes+1; //Bytes to store one non-constant block data.
    size_t mSize = sizeof(float)+1+ncBytes; //Number of bytes for each data block's metadata.
    size_t msz = (1+mSize) * nbBlocks * sizeof(unsigned char);
    size_t mbsz = sizeof(float) * nbEle * sizeof(unsigned char);

    unsigned char *meta = (unsigned char*)malloc(msz);
    short *offsets = (short*)malloc(nbBlocks*sizeof(short));
    unsigned char *midBytes = (unsigned char*)malloc(mbsz);

	unsigned char* d_meta;
	unsigned char* d_midBytes;
	short* d_offsets;
    checkCudaErrors(cudaMalloc((void**)&d_meta, msz)); 
    //checkCudaErrors(cudaMemcpy(d_meta, meta, msz, cudaMemcpyHostToDevice)); 
    checkCudaErrors(cudaMemset(d_meta, 0, msz));
    checkCudaErrors(cudaMalloc((void**)&d_offsets, nbBlocks*sizeof(short))); 
    checkCudaErrors(cudaMemset(d_offsets, 0, nbBlocks*sizeof(short)));
    checkCudaErrors(cudaMalloc((void**)&d_midBytes, mbsz)); 
    checkCudaErrors(cudaMemset(d_midBytes, 0, mbsz));

    timer_GPU.StartCounter();
    dim3 dimBlock(32, blockSize/32);
    dim3 dimGrid(65536, 1);
    const int sMemsize = blockSize * sizeof(float) + dimBlock.y * sizeof(int);
    compress_float<<<dimGrid, dimBlock, sMemsize>>>(d_oriData, d_meta, d_offsets, d_midBytes, absErrBound, blockSize, nbBlocks, mSize);
    cudaError_t err = cudaGetLastError();        // Get error code
    printf("CUDA Error: %s\n", cudaGetErrorString(err));
    printf("GPU compression timing: %f ms\n", timer_GPU.GetCounter());
    checkCudaErrors(cudaMemcpy(meta, d_meta, msz, cudaMemcpyDeviceToHost)); 
    checkCudaErrors(cudaMemcpy(offsets, d_offsets, nbBlocks*sizeof(short), cudaMemcpyDeviceToHost)); 
    checkCudaErrors(cudaMemcpy(midBytes, d_midBytes, mbsz, cudaMemcpyDeviceToHost)); 

    size_t maxPreservedBufferSize = sizeof(float)*nbEle;
    unsigned char* outBytes = (unsigned char*)malloc(maxPreservedBufferSize);
    memset(outBytes, 0, maxPreservedBufferSize);

    *outSize = _post_proc(oriData, meta, offsets, midBytes, outBytes, nbEle, blockSize);

    free(meta);
    free(offsets);
    free(midBytes);
    checkCudaErrors(cudaFree(d_meta));
    checkCudaErrors(cudaFree(d_offsets));
    checkCudaErrors(cudaFree(d_midBytes));
    return outBytes;
}

unsigned char* cuSZx_fast_compress_args_artifact_aware(float *oriData, size_t *outSize, float absErrBound, size_t nbEle, int blockSize, int xdim, int ydim, int zdim)
{

	float* d_oriData;
    unsigned char *artifact_type, *d_artifact_type;
    cudaMalloc((void**)&d_oriData, sizeof(float)*nbEle); 
    cudaMemcpy(d_oriData, oriData, sizeof(float)*nbEle, cudaMemcpyHostToDevice); 

	size_t nbBlocks = nbEle/blockSize;
	size_t remainCount = nbEle%blockSize;
	size_t actualNBBlocks = remainCount==0 ? nbBlocks : nbBlocks+1;

    size_t ncBytes = blockSize/4;
    //ncBytes = (blockSize+1)%4==0 ? ncBytes : ncBytes+1; //Bytes to store one non-constant block data.
    size_t mSize = sizeof(float)+1+ncBytes; //Number of bytes for each data block's metadata.
    size_t msz = (1+mSize) * nbBlocks * sizeof(unsigned char);
    size_t mbsz = sizeof(float) * nbEle * sizeof(unsigned char);

    unsigned char *meta = (unsigned char*)malloc(msz);
    short *offsets = (short*)malloc(nbBlocks*sizeof(short));
    unsigned char *midBytes = (unsigned char*)malloc(mbsz);

	unsigned char* d_meta;
	unsigned char* d_midBytes;
	short* d_offsets;
    checkCudaErrors(cudaMalloc((void**)&d_meta, msz)); 
    //checkCudaErrors(cudaMemcpy(d_meta, meta, msz, cudaMemcpyHostToDevice)); 
    checkCudaErrors(cudaMemset(d_meta, 0, msz));
    checkCudaErrors(cudaMalloc((void**)&d_offsets, nbBlocks*sizeof(short))); 
    checkCudaErrors(cudaMemset(d_offsets, 0, nbBlocks*sizeof(short)));
    checkCudaErrors(cudaMalloc((void**)&d_midBytes, mbsz)); 
    checkCudaErrors(cudaMemset(d_midBytes, 0, mbsz));
    checkCudaErrors(cudaMalloc(&d_artifact_type, sizeof(char)*nbBlocks));

    timer_GPU.StartCounter();
    dim3 dimBlock(32, blockSize/32);
    dim3 dimGrid(65536, 1);
    const int sMemsize = blockSize * sizeof(float) + dimBlock.y * sizeof(int);
    compress_float<<<dimGrid, dimBlock, sMemsize>>>(d_oriData, d_meta, d_offsets, d_midBytes, absErrBound, blockSize, nbBlocks, mSize);
    cudaError_t err = cudaGetLastError();        // Get error code
    printf("CUDA Error: %s\n", cudaGetErrorString(err));
    printf("GPU compression timing: %f ms\n", timer_GPU.GetCounter());
    checkCudaErrors(cudaMemcpy(meta, d_meta, msz, cudaMemcpyDeviceToHost)); 
    checkCudaErrors(cudaMemcpy(offsets, d_offsets, nbBlocks*sizeof(short), cudaMemcpyDeviceToHost)); 
    checkCudaErrors(cudaMemcpy(midBytes, d_midBytes, mbsz, cudaMemcpyDeviceToHost)); 
    
    artifact_analysis<<<40,256>>>(d_oriData, blockSize, absErrBound, nbEle, d_meta, xdim, ydim, zdim, d_artifact_type);

    cudaDeviceSynchronize();
    checkCudaErrors(cudaMemcpy(artifact_type, d_artifact_type,sizeof(char)*nbBlocks, cudaMemcpyDeviceToHost));
    size_t maxPreservedBufferSize = sizeof(float)*nbEle;
    unsigned char* outBytes = (unsigned char*)malloc(maxPreservedBufferSize);
    memset(outBytes, 0, maxPreservedBufferSize);

    *outSize = _post_proc(oriData, meta, offsets, midBytes, outBytes, nbEle, blockSize);

    FILE *file = fopen("artifacts.bin", "wb");
    fwrite(artifact_type, 1, nbBlocks, file);
    fclose(file);

    free(artifact_type);
    free(meta);
    free(offsets);
    free(midBytes);
    cudaFree(d_artifact_type);
    checkCudaErrors(cudaFree(d_meta));
    checkCudaErrors(cudaFree(d_offsets));
    checkCudaErrors(cudaFree(d_midBytes));
    return outBytes;
}

void cuSZx_fast_decompress_args_unpredictable_blocked_float(float** newData, size_t nbEle, unsigned char* cmpBytes)
{
	*newData = (float*)malloc(sizeof(float)*nbEle);
    memset(*newData, 0, sizeof(float)*nbEle);
	
	unsigned char* r = cmpBytes;
	r += 4;
	int blockSize = r[0];  //get block size
	r++;
	size_t nbConstantBlocks = bytesToLong_bigEndian(r); //get number of constant blocks
	r += sizeof(size_t);
		
	size_t nbBlocks = nbEle/blockSize;
	size_t ncBlocks = nbBlocks - nbConstantBlocks; //get number of constant blocks
	size_t stateNBBytes = nbBlocks%8==0 ? nbBlocks/8 : nbBlocks/8+1;
    size_t ncLeading = blockSize/4;
    size_t mSize = sizeof(float)+1+ncLeading; //Number of bytes for each data block's metadata.
	unsigned char* stateArray = (unsigned char*)malloc(nbBlocks);
	float* constantMedianArray = (float*)malloc(nbConstantBlocks*sizeof(float));			
	unsigned char* data = (unsigned char*)malloc(ncBlocks*blockSize*sizeof(float));
    memset(data, 0, ncBlocks*blockSize*sizeof(float));
		
	convertByteArray2IntArray_fast_1b_args(nbBlocks, r, stateNBBytes, stateArray); //get the stateArray
	
	r += stateNBBytes;
	size_t i = 0, j = 0, k = 0; //k is used to keep track of constant block index
    memcpy((*newData)+nbBlocks*blockSize, r, (nbEle%blockSize)*sizeof(float));
    r += (nbEle%blockSize)*sizeof(float);
	float* fr = (float*)r; //fr is the starting address of constant median values.
	for(i = 0;i < nbConstantBlocks;i++, j+=4) //get the median values for constant-value blocks
		constantMedianArray[i] = fr[i];
    r += nbConstantBlocks*sizeof(float);
    unsigned char* p = r + ncBlocks * sizeof(short);
    for(i = 0;i < ncBlocks;i++){
        int leng = (int)bytesToShort(r)+mSize;
        r += sizeof(short);
        if (leng > blockSize*sizeof(float))
        {
            printf("Warning: compressed block is larger than the original block!\n");
            exit(0);
        }
        memcpy(data+i*blockSize*sizeof(float), p, leng);
        p += leng;
    } 

    unsigned char* d_data;
    checkCudaErrors(cudaMalloc((void**)&d_data, ncBlocks*blockSize*sizeof(float))); 
    checkCudaErrors(cudaMemcpy(d_data, data, ncBlocks*blockSize*sizeof(float), cudaMemcpyHostToDevice)); 

    timer_GPU.StartCounter();
    dim3 dimBlock(32, blockSize/32);
    dim3 dimGrid(65536, 1);
    const int sMemsize = blockSize * sizeof(float) + dimBlock.y * sizeof(int);
    decompress_float<<<dimGrid, dimBlock, sMemsize>>>(d_data, blockSize, ncBlocks, mSize);
    cudaError_t err = cudaGetLastError();        // Get error code
    printf("CUDA Error: %s\n", cudaGetErrorString(err));
    printf("GPU decompression timing: %f ms\n", timer_GPU.GetCounter());
    checkCudaErrors(cudaMemcpy(data, d_data, ncBlocks*blockSize*sizeof(float), cudaMemcpyDeviceToHost)); 
    float* fdata = (float*)data;

    int nb=0, nc=0;
    for (i=0;i<nbBlocks;i++){
        if (stateArray[i]==0){
            float Median = constantMedianArray[nb];
            if (Median>1) printf("data%i:%f\n",i, Median);
            for (j=0;j<blockSize;j++)
                *((*newData)+i*blockSize+j) = Median;
            nb++;
        }else{
            for (j=0;j<blockSize;j++)
                *((*newData)+i*blockSize+j) = fdata[nc*blockSize+j];
            nc++;
        }
    }

	free(stateArray);
	free(constantMedianArray);
	free(data);
    checkCudaErrors(cudaFree(d_data));

}

void cuSZx_fast_decompress_args_unpredictable_blocked_float_artifacts(float** newData, size_t nbEle, unsigned char* cmpBytes)
{
	*newData = (float*)malloc(sizeof(float)*nbEle);
    memset(*newData, 0, sizeof(float)*nbEle);
	
	unsigned char* r = cmpBytes;
	r += 4;
	int blockSize = r[0];  //get block size
	r++;
	size_t nbConstantBlocks = bytesToLong_bigEndian(r); //get number of constant blocks
	r += sizeof(size_t);
		
	size_t nbBlocks = nbEle/blockSize;
	size_t ncBlocks = nbBlocks - nbConstantBlocks; //get number of constant blocks
	size_t stateNBBytes = nbBlocks%8==0 ? nbBlocks/8 : nbBlocks/8+1;
    size_t ncLeading = blockSize/4;
    size_t mSize = sizeof(float)+1+ncLeading; //Number of bytes for each data block's metadata.
	unsigned char* stateArray = (unsigned char*)malloc(nbBlocks);
	float* constantMedianArray = (float*)malloc(nbConstantBlocks*sizeof(float));			
	unsigned char* data = (unsigned char*)malloc(ncBlocks*blockSize*sizeof(float));
    memset(data, 0, ncBlocks*blockSize*sizeof(float));
	
    unsigned char *artifact_type = (unsigned char*)malloc(nbBlocks*sizeof(char));
    FILE *file = fopen("artifacts.bin", "wb");
    fread(artifact_type, 1, nbBlocks, file);
    fclose(file);

	convertByteArray2IntArray_fast_1b_args(nbBlocks, r, stateNBBytes, stateArray); //get the stateArray
	
	r += stateNBBytes;
	size_t i = 0, j = 0, k = 0; //k is used to keep track of constant block index
    memcpy((*newData)+nbBlocks*blockSize, r, (nbEle%blockSize)*sizeof(float));
    r += (nbEle%blockSize)*sizeof(float);
	float* fr = (float*)r; //fr is the starting address of constant median values.
	for(i = 0;i < nbConstantBlocks;i++, j+=4) //get the median values for constant-value blocks
		constantMedianArray[i] = fr[i];
    r += nbConstantBlocks*sizeof(float);
    unsigned char* p = r + ncBlocks * sizeof(short);
    for(i = 0;i < ncBlocks;i++){
        int leng = (int)bytesToShort(r)+mSize;
        r += sizeof(short);
        if (leng > blockSize*sizeof(float))
        {
            printf("Warning: compressed block is larger than the original block!\n");
            exit(0);
        }
        memcpy(data+i*blockSize*sizeof(float), p, leng);
        p += leng;
    } 

    unsigned char* d_data;
    checkCudaErrors(cudaMalloc((void**)&d_data, ncBlocks*blockSize*sizeof(float))); 
    checkCudaErrors(cudaMemcpy(d_data, data, ncBlocks*blockSize*sizeof(float), cudaMemcpyHostToDevice)); 

    timer_GPU.StartCounter();
    dim3 dimBlock(32, blockSize/32);
    dim3 dimGrid(65536, 1);
    const int sMemsize = blockSize * sizeof(float) + dimBlock.y * sizeof(int);
    decompress_float<<<dimGrid, dimBlock, sMemsize>>>(d_data, blockSize, ncBlocks, mSize);
    cudaError_t err = cudaGetLastError();        // Get error code
    printf("CUDA Error: %s\n", cudaGetErrorString(err));
    printf("GPU decompression timing: %f ms\n", timer_GPU.GetCounter());
    checkCudaErrors(cudaMemcpy(data, d_data, ncBlocks*blockSize*sizeof(float), cudaMemcpyDeviceToHost)); 
    float* fdata = (float*)data;

    int nb=0, nc=0;
    for (i=0;i<nbBlocks;i++){
        if (stateArray[i]==0){
            float Median = constantMedianArray[nb];
            if (Median>1) printf("data%i:%f\n",i, Median);
            for (j=0;j<blockSize;j++)
                *((*newData)+i*blockSize+j) = Median;
            nb++;
        }else{
            for (j=0;j<blockSize;j++)
                *((*newData)+i*blockSize+j) = fdata[nc*blockSize+j];
            nc++;
        }
    }

	free(stateArray);
	free(constantMedianArray);
	free(data);
    checkCudaErrors(cudaFree(d_data));

}