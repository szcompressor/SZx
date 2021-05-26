#include "cuSZx_entry.h"
#include "timingGPU.h"

TimingGPU timer_GPU;

unsigned char* cuSZx_fast_compress_args_unpredictable_blocked_float(float *oriData, size_t *outSize, float absErrBound, size_t nbEle, int blockSize, unsigned char *test)
{

	float* d_oriData;
    cudaMalloc((void**)&d_oriData, sizeof(float)*nbEle); 
    cudaMemcpy(d_oriData, oriData, sizeof(float)*nbEle, cudaMemcpyHostToDevice); 

	size_t nbBlocks = nbEle/blockSize;
	size_t remainCount = nbEle%blockSize;
	size_t actualNBBlocks = remainCount==0 ? nbBlocks : nbBlocks+1;

    size_t ncBytes = blockSize/4;
    //ncBytes = (blockSize+1)%4==0 ? ncBytes : ncBytes+1; //Bytes to store one non-constant block data.
    size_t mSize = 1+sizeof(float)+1+ncBytes+sizeof(unsigned int); //Number of bytes for each data block's metadata.
    size_t msz = mSize * nbBlocks * sizeof(unsigned char);
    size_t mbsz = sizeof(float) * nbEle * sizeof(unsigned char);

    unsigned char *meta = (unsigned char*)malloc(msz);
    unsigned char *midBytes = (unsigned char*)malloc(mbsz);

	unsigned char* d_meta;
	unsigned char* d_midBytes;
    checkCudaErrors(cudaMalloc((void**)&d_meta, msz)); 
    //checkCudaErrors(cudaMemcpy(d_meta, meta, msz, cudaMemcpyHostToDevice)); 
    checkCudaErrors(cudaMemset(d_meta, 0, msz));
    checkCudaErrors(cudaMalloc((void**)&d_midBytes, mbsz)); 
    //checkCudaErrors(cudaMemset(d_midBytes, 0, mbsz));
    //cudaMemcpy(dresults, results, sizeof(unsigned char)*reSize*nbBlocks, cudaMemcpyHostToDevice); 
    for (int i=0; i<nbBlocks; i++) 
        printf ("test:%d\n", i);

    //timer_GPU.StartCounter();
    //dim3 dimBlock(32, blockSize/32);
    //dim3 dimGrid(512, 1);
    //const int sMemsize = 10 * dimBlock.y * sizeof(double);
    ////compress_float<<<dimGrid, dimBlock, sMemsize>>>(d_oriData, d_meta, d_midBytes, absErrBound, blockSize, nbBlocks, mSize);
    //cudaError_t err = cudaGetLastError();        // Get error code
    //printf("CUDA Error: %s\n", cudaGetErrorString(err));
    //checkCudaErrors(cudaMemcpy(meta, d_meta, msz, cudaMemcpyDeviceToHost)); 

    //for (int i=0; i<nbBlocks; i++){ 
    //    //if (meta[i*mSize]!=test[i]) 
    //    //printf("state %d : %u\n", i, test[i]);
    //    printf ("test:%d\n", nbBlocks);
    //}
}
