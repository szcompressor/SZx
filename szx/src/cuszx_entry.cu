#include "cuszx_entry.h"
#include "szx_defines.h"
#include "szx_BytesToolkit.h"
#include "szx_TypeManager.h"
#include "timingGPU.h"

#define SPARSITY_LEVEL 0.25

TimingGPU timer_GPU;
void bin(unsigned n)
{
    unsigned i;
    for (i = 1 << 31; i > 0; i = i / 2)
        (n & i) ? printf("1") : printf("0");
}

size_t convert_state_to_out(unsigned char* meta, size_t length, unsigned char *result){
    size_t out_length;

    if(length%4==0)
		out_length = length/4;
	else
		out_length = length/4+1;

    for (size_t i = 0; i < out_length; i++)
    {
        uint8_t tmp = 0;

        for (size_t j = 0; j < 4; j++)
        {
            if (i*4 + j < length)
            {
                tmp |= (0x03 & meta[i*4+j]) << 2*j;
            }
            
        }
        result[i] = tmp;
    }
    return out_length;
}

size_t convert_block2_to_out(unsigned char *result, uint32_t numBlocks, uint64_t num_sig, uint32_t *blk_idx, float *blk_vals, uint8_t *blk_subidx){
    size_t out_length = 0;
    memcpy(result, blk_idx, numBlocks*sizeof(uint32_t));
    out_length += numBlocks*4;
    memcpy(result+out_length, blk_vals, num_sig*sizeof(float));
    out_length += num_sig*sizeof(float);
    memcpy(result+out_length, blk_subidx, num_sig*sizeof(uint8_t));
    out_length += num_sig*sizeof(uint8_t);
    
    return out_length;
}

int _post_proc(float *oriData, unsigned char *meta, short *offsets, unsigned char *midBytes, unsigned char *outBytes, size_t nbEle, int blockSize, uint64_t num_sig, uint32_t *blk_idx, float *blk_vals, uint8_t *blk_subidx)
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
    unsigned char* r_old = outBytes; 
	r[0] = SZx_VER_MAJOR;
	r[1] = SZx_VER_MINOR;
	r[2] = 1;
	r[3] = 0; // indicates this is not a random access version
	r[4] = (unsigned char)blockSize;
	r=r+5; //1 byte
	sizeToBytes(r, nbConstantBlocks);
	r += sizeof(size_t); 
	r += convert_state_to_out(meta, nbBlocks, r);
    r += convert_block2_to_out(r, nbBlocks,num_sig, blk_idx, blk_vals, blk_subidx);
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

    // return out_size;
    return (uint32_t) (r-r_old);
}

unsigned char* cuSZx_fast_compress_args_unpredictable_blocked_float(float *oriData, size_t *outSize, float absErrBound, size_t nbEle, int blockSize, float threshold)
{
    float sparsity_level = SPARSITY_LEVEL;
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

    uint32_t *blk_idx, *d_blk_idx;
    uint8_t *blk_subidx, *d_blk_subidx;
    float *blk_vals, *d_blk_vals;
    uint64_t *num_sig, *d_num_sig;

    checkCudaErrors(cudaMalloc((void **)&d_num_sig, sizeof(uint64_t)));
    num_sig = (uint64_t *)malloc(sizeof(uint64_t));
    checkCudaErrors(cudaMalloc((void **)&d_blk_idx, nbBlocks*sizeof(uint32_t)));
    // blk_idx = malloc()
    checkCudaErrors(cudaMalloc((void **)&d_blk_subidx, nbEle*sizeof(uint8_t)));

    checkCudaErrors(cudaMalloc((void **)&d_blk_vals, nbEle*sizeof(float)));

    checkCudaErrors(cudaMalloc((void**)&d_meta, msz)); 
    //checkCudaErrors(cudaMemcpy(d_meta, meta, msz, cudaMemcpyHostToDevice)); 
    checkCudaErrors(cudaMemset(d_meta, 0, msz));
    checkCudaErrors(cudaMalloc((void**)&d_offsets, nbBlocks*sizeof(short))); 
    checkCudaErrors(cudaMemset(d_offsets, 0, nbBlocks*sizeof(short)));
    checkCudaErrors(cudaMalloc((void**)&d_midBytes, mbsz)); 
    checkCudaErrors(cudaMemset(d_midBytes, 0, mbsz));

    timer_GPU.StartCounter();
    apply_threshold<<<80,256>>>(d_oriData, threshold, nbEle);
    cudaDeviceSynchronize();
    dim3 dimBlock(32, blockSize/32);
    dim3 dimGrid(65536, 1);
    const int sMemsize = blockSize * sizeof(float) + dimBlock.y * sizeof(int);
    compress_float<<<dimGrid, dimBlock, sMemsize>>>(d_oriData, d_meta, d_offsets, d_midBytes, absErrBound, blockSize, nbBlocks, mSize, sparsity_level, d_blk_idx, d_blk_subidx,d_blk_vals);
    cudaError_t err = cudaGetLastError();        // Get error code
    printf("CUDA Error: %s\n", cudaGetErrorString(err));
    printf("GPU compression timing: %f ms\n", timer_GPU.GetCounter());
    cudaDeviceSynchronize();
    get_numsig<<<1,1>>>(d_num_sig);
    cudaDeviceSynchronize();

    checkCudaErrors(cudaMemcpy(num_sig, d_num_sig, sizeof(uint64_t), cudaMemcpyDeviceToHost));

    blk_idx = (uint32_t *)malloc(nbBlocks*sizeof(uint32_t));
    blk_vals= (float *)malloc((*num_sig)*sizeof(float));
    blk_subidx = (uint8_t *)malloc((*num_sig)*sizeof(uint8_t));

    checkCudaErrors(cudaMemcpy(meta, d_meta, msz, cudaMemcpyDeviceToHost)); 
    checkCudaErrors(cudaMemcpy(offsets, d_offsets, nbBlocks*sizeof(short), cudaMemcpyDeviceToHost)); 
    checkCudaErrors(cudaMemcpy(midBytes, d_midBytes, mbsz, cudaMemcpyDeviceToHost)); 
    
    
    checkCudaErrors(cudaMemcpy(blk_idx, d_blk_idx, nbBlocks*sizeof(uint32_t), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(blk_vals,d_blk_vals, (*num_sig)*sizeof(float), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(blk_subidx,d_blk_subidx, (*num_sig)*sizeof(uint8_t), cudaMemcpyDeviceToHost));

    size_t maxPreservedBufferSize = sizeof(float)*nbEle;
    unsigned char* outBytes = (unsigned char*)malloc(maxPreservedBufferSize);
    memset(outBytes, 0, maxPreservedBufferSize);

    *outSize = _post_proc(oriData, meta, offsets, midBytes, outBytes, nbEle, blockSize, *num_sig, blk_idx, blk_vals, blk_subidx);

    free(blk_idx);
    free(blk_subidx);
    free(blk_vals);
    free(meta);
    free(offsets);
    free(midBytes);
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
