#include "cuszx_entry.h"
#include "szx_defines.h"
#include "szx_BytesToolkit.h"
#include "szx_TypeManager.h"
#include "timingGPU.h"
#include "szx.h"

#define SPARSITY_LEVEL 0.25

TimingGPU timer_GPU;
void bin(unsigned n)
{
    unsigned i;
    for (i = 1 << 31; i > 0; i = i / 2)
        (n & i) ? printf("1") : printf("0");
}

__host__ __device__ size_t convert_state_to_out(unsigned char* meta, size_t length, unsigned char *result){
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

// nbBlocks, r, stateNBBytes, stateArray
__host__ __device__ size_t convert_out_to_state(size_t nbBlocks, unsigned char* cmp, unsigned char* out_state){
    size_t state_length;
    if(nbBlocks%4==0)
		state_length = nbBlocks/4;
	else
		state_length = nbBlocks/4+1;

    for (size_t i = 0; i < state_length; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            if (4*i + j < nbBlocks)
            {
                out_state[4*i + j]= (cmp[i] >> 2*j) & 0x03;
            }
            
        }
    }
    return nbBlocks;
}

__host__ __device__ size_t convert_block2_to_out(unsigned char *result, uint32_t numBlocks, uint64_t num_sig, uint32_t *blk_idx, float *blk_vals, uint8_t *blk_subidx, uint8_t *blk_sig){
    size_t out_length = 0;
    memcpy(result, blk_idx, numBlocks*sizeof(uint32_t));
    out_length += numBlocks*4;
    memcpy(result+out_length, blk_vals, num_sig*sizeof(float));
    out_length += num_sig*sizeof(float);
    memcpy(result+out_length, blk_subidx, num_sig*sizeof(uint8_t));
    out_length += num_sig*sizeof(uint8_t);
    memcpy(result+out_length, blk_sig, numBlocks*sizeof(uint8_t));
    out_length+= numBlocks*sizeof(uint8_t);

    return out_length;
}

__host__ __device__ size_t convert_out_to_block2(unsigned char *in_cmp, uint32_t numBlocks, uint64_t num_sig, uint32_t *blk_idx, float *blk_vals, uint8_t *blk_subidx, uint8_t *blk_sig){
    size_t out_length = 0;
    memcpy(blk_idx, in_cmp, numBlocks*sizeof(uint32_t));
    out_length += numBlocks*4;
    memcpy(blk_vals, in_cmp+out_length,num_sig*sizeof(float));
    out_length += num_sig*sizeof(float);
    memcpy(blk_subidx, in_cmp+out_length, num_sig*sizeof(uint8_t));
    out_length += num_sig*sizeof(uint8_t);
    memcpy(blk_sig, in_cmp+out_length, numBlocks*sizeof(uint8_t));
    out_length += numBlocks*sizeof(uint8_t);

    return out_length;
}

int _post_proc(float *oriData, unsigned char *meta, short *offsets, unsigned char *midBytes, unsigned char *outBytes, size_t nbEle, int blockSize, uint64_t num_sig, uint32_t *blk_idx, float *blk_vals, uint8_t *blk_subidx, uint8_t *blk_sig)
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
    int s0 = 0;
    int s1 = 0;
    int s2 = 0;
    int s3 = 0;
    for (int i=0; i<nbBlocks; i++){
        if (meta[i]==0 || meta[i]==1 || meta[i] == 2) nbConstantBlocks++;
        else out_size += 1+(blockSize/4)+offsets[i];
    
    	if(meta[i]==0) s0++;
    	if(meta[i]==1) s1++;
    	if(meta[i]==2) s2++;
    	if(meta[i]==3) s3++;
    }
    printf("%d %d %d %d\n", s0, s1, s2, s3);
    out_size += (nbBlocks-nbConstantBlocks)*sizeof(short)+(nbEle%blockSize)*sizeof(float);

    //outBytes = (unsigned char*)malloc(out_size);
  //  printf("accessing outbytes now...\n");
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
    sizeToBytes(r, (size_t) num_sig);
    r += sizeof(size_t); 
	r += convert_state_to_out(meta, nbBlocks, r);
    r += convert_block2_to_out(r, nbBlocks,num_sig, blk_idx, blk_vals, blk_subidx, blk_sig);
    memcpy(r, oriData+nbBlocks*blockSize, (nbEle%blockSize)*sizeof(float));
    r += (nbEle%blockSize)*sizeof(float);
    unsigned char* c = r;
    unsigned char* o = c+nbConstantBlocks*sizeof(float);
    unsigned char* nc = o+(nbBlocks-nbConstantBlocks)*sizeof(short);
    for (int i=0; i<nbBlocks; i++){
        
        if (meta[i]==0 || meta[i] == 1){
	    memcpy(c, meta+(nbBlocks+i*mSize), sizeof(float));
            c += sizeof(float);
        }else if(meta[i] == 3){
            shortToBytes(o, offsets[i]);
	   
            o += sizeof(short);
            memcpy(nc, meta+(nbBlocks+i*mSize), mSize);
            
	    nc += mSize; 
            memcpy(nc, midBytes+(i*blockSize*sizeof(float)), offsets[i]);
            
	    nc += offsets[i];
	   
        } 
    }

    // return out_size;
    return (uint32_t) (nc-r_old);
}

unsigned char* cuSZx_fast_compress_args_unpredictable_blocked_float(float *oriData, size_t *outSize, float absErrBound, size_t nbEle, int blockSize, float threshold)
{
//    printf("tr thresh abs %f %f\n", threshold, absErrBound);
  //  printf("first: %f %f %f\n", oriData[0], oriData[1], oriData[2]);
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
    uint8_t *blk_sig, *d_blk_sig;
    uint8_t *blk_subidx, *d_blk_subidx;
    float *blk_vals, *d_blk_vals;
    uint64_t *num_sig, *d_num_sig;

    checkCudaErrors(cudaMalloc((void **)&d_num_sig, sizeof(uint64_t)));
    num_sig = (uint64_t *)malloc(sizeof(uint64_t));
    checkCudaErrors(cudaMalloc((void **)&d_blk_idx, nbBlocks*sizeof(uint32_t)));
    // blk_idx = malloc()
    checkCudaErrors(cudaMalloc((void **)&d_blk_subidx, nbEle*sizeof(uint8_t)));

    checkCudaErrors(cudaMalloc((void **)&d_blk_vals, nbEle*sizeof(float)));

    checkCudaErrors(cudaMalloc((void **)&d_blk_sig, nbBlocks*sizeof(uint8_t)));

    checkCudaErrors(cudaMalloc((void**)&d_meta, msz)); 
    //checkCudaErrors(cudaMemcpy(d_meta, meta, msz, cudaMemcpyHostToDevice)); 
    checkCudaErrors(cudaMemset(d_meta, 0, msz));
    checkCudaErrors(cudaMalloc((void**)&d_offsets, nbBlocks*sizeof(short))); 
    checkCudaErrors(cudaMemset(d_offsets, 0, nbBlocks*sizeof(short)));
    checkCudaErrors(cudaMalloc((void**)&d_midBytes, mbsz)); 
    checkCudaErrors(cudaMemset(d_midBytes, 0, mbsz));

    timer_GPU.StartCounter();
    // apply_threshold<<<80,256>>>(d_oriData, threshold, nbEle);
    // cudaDeviceSynchronize();
    dim3 dimBlock(32, blockSize/32);
    dim3 dimGrid(65536, 1);
    const int sMemsize = blockSize * sizeof(float) + dimBlock.y * sizeof(int);
    compress_float<<<dimGrid, dimBlock, sMemsize>>>(d_oriData, d_meta, d_offsets, d_midBytes, absErrBound, blockSize, nbBlocks, mSize, sparsity_level, d_blk_idx, d_blk_subidx,d_blk_vals, threshold, d_blk_sig);
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
    blk_sig = (uint8_t *)malloc(nbBlocks*sizeof(uint8_t));

    checkCudaErrors(cudaMemcpy(meta, d_meta, msz, cudaMemcpyDeviceToHost)); 
    checkCudaErrors(cudaMemcpy(offsets, d_offsets, nbBlocks*sizeof(short), cudaMemcpyDeviceToHost)); 
    checkCudaErrors(cudaMemcpy(midBytes, d_midBytes, mbsz, cudaMemcpyDeviceToHost)); 
    
    
    checkCudaErrors(cudaMemcpy(blk_idx, d_blk_idx, nbBlocks*sizeof(uint32_t), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(blk_vals,d_blk_vals, (*num_sig)*sizeof(float), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(blk_subidx,d_blk_subidx, (*num_sig)*sizeof(uint8_t), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(blk_sig,d_blk_sig, (nbBlocks)*sizeof(uint8_t), cudaMemcpyDeviceToHost));

    size_t maxPreservedBufferSize = sizeof(float)*nbEle;
    unsigned char* outBytes = (unsigned char*)malloc(maxPreservedBufferSize);
    memset(outBytes, 0, maxPreservedBufferSize);

    outSize = (size_t *)malloc(sizeof(size_t));
    //outSize[0] = _post_proc(oriData, meta, offsets, midBytes, outBytes, nbEle, blockSize, *num_sig, blk_idx, blk_vals, blk_subidx, blk_sig);

    *outSize = _post_proc(oriData, meta, offsets, midBytes, outBytes, nbEle, blockSize, *num_sig, blk_idx, blk_vals, blk_subidx, blk_sig);
//    printf("Beginning free\n");
    printf("outsize %p \n", outBytes);
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
    uint32_t *blk_idx, *d_blk_idx;
    uint8_t *blk_subidx, *d_blk_subidx;
    uint8_t *blk_sig, *d_blk_sig;
    float *blk_vals, *d_blk_vals;
    size_t num_sig, *d_num_sig;

	*newData = (float*)malloc(sizeof(float)*nbEle);
    memset(*newData, 0, sizeof(float)*nbEle);
	
	unsigned char* r = cmpBytes;
	r += 4;
	int blockSize = r[0];  //get block size
	if(blockSize == 0)blockSize = 256;
	r++;
	size_t nbConstantBlocks = bytesToLong_bigEndian(r); //get number of constant blocks
	r += sizeof(size_t);
	num_sig = bytesToSize(r);
    r += sizeof(size_t);
	size_t nbBlocks = nbEle/blockSize;
    size_t ncBlocks = 0;
    size_t num_state2_blks = 0;
	// size_t ncBlocks = nbBlocks - nbConstantBlocks; //get number of constant blocks
	size_t stateNBBytes = nbBlocks%4==0 ? nbBlocks/4 : nbBlocks/4+1;
    size_t ncLeading = blockSize/4;
    size_t mSize = sizeof(float)+1+ncLeading; //Number of bytes for each data block's metadata.
	unsigned char* stateArray = (unsigned char*)malloc(nbBlocks);
    unsigned char* d_stateArray;
    cudaMalloc(&d_stateArray, nbBlocks);
	float* constantMedianArray = (float*)malloc(nbConstantBlocks*sizeof(float));			
	
    

    blk_idx = (uint32_t *)malloc(nbBlocks*sizeof(uint32_t));
    blk_vals= (float *)malloc((num_sig)*sizeof(float));
    blk_subidx = (uint8_t *)malloc((num_sig)*sizeof(uint8_t));
    blk_sig = (uint8_t *)malloc(nbBlocks*sizeof(uint8_t));

	printf("Converting state array\n");
    convert_out_to_state(nbBlocks, r, stateArray);
	// convertByteArray2IntArray_fast_1b_args(nbBlocks, r, stateNBBytes, stateArray); //get the stateArray
	for (size_t i = 0; i < nbBlocks; i++)
    {
        if (stateArray[i] == 2)
        {
            num_state2_blks++;
        }else if(stateArray[i] == 3){
            ncBlocks++;
        }
    }
    
	r += stateNBBytes;
    unsigned char* data = (unsigned char*)malloc(ncBlocks*blockSize*sizeof(float));
    memset(data, 0, ncBlocks*blockSize*sizeof(float));
    printf("converting block vals\n");
    size_t to_add = convert_out_to_block2(r, nbBlocks, (uint64_t)num_sig, blk_idx, blk_vals, blk_subidx, blk_sig);
    r+= to_add;
    // checkCudaErrors(cudaMalloc((void **)&d_num_sig, sizeof(uint64_t)));
    // num_sig = (uint64_t *)malloc(sizeof(uint64_t));
    checkCudaErrors(cudaMalloc((void **)&d_blk_idx, nbBlocks*sizeof(uint32_t)));
    // blk_idx = malloc()
    checkCudaErrors(cudaMalloc((void **)&d_blk_subidx, num_sig*sizeof(uint8_t)));

    checkCudaErrors(cudaMalloc((void **)&d_blk_vals, num_sig*sizeof(float)));

    checkCudaErrors(cudaMalloc((void **)&d_blk_sig, nbBlocks*sizeof(uint8_t)));

    checkCudaErrors(cudaMemcpy(d_blk_idx, blk_idx, nbBlocks*sizeof(uint32_t), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_blk_vals, blk_vals, (num_sig)*sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_blk_subidx, blk_subidx, (num_sig)*sizeof(uint8_t), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_stateArray, stateArray, nbBlocks, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_blk_sig, blk_sig, nbBlocks*sizeof(uint8_t), cudaMemcpyHostToDevice));


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
    float *d_newdata;
    checkCudaErrors(cudaMalloc((void**)&d_data, ncBlocks*blockSize*sizeof(float))); 
    checkCudaErrors(cudaMemcpy(d_data, data, ncBlocks*blockSize*sizeof(float), cudaMemcpyHostToDevice)); 
    checkCudaErrors(cudaMalloc(&d_newdata, nbBlocks*blockSize*sizeof(float)));

    timer_GPU.StartCounter();
    dim3 dimBlock(32, blockSize/32);
    dim3 dimGrid(65536, 1);
    const int sMemsize = blockSize * sizeof(float) + dimBlock.y * sizeof(int);
    decompress_state2<<<nbBlocks, 64>>>(d_newdata, d_stateArray,d_blk_idx, d_blk_vals, d_blk_subidx,blockSize, d_blk_sig);
    decompress_float<<<dimGrid, dimBlock, sMemsize>>>(d_data, blockSize, ncBlocks, mSize);
    cudaError_t err = cudaGetLastError();        // Get error code
    printf("CUDA Error: %s\n", cudaGetErrorString(err));
    printf("GPU decompression timing: %f ms\n", timer_GPU.GetCounter());
    cudaDeviceSynchronize();
    checkCudaErrors(cudaMemcpy(data, d_data, ncBlocks*blockSize*sizeof(float), cudaMemcpyDeviceToHost)); 
    checkCudaErrors(cudaMemcpy(*newData, d_newdata, nbBlocks*blockSize*sizeof(float), cudaMemcpyDeviceToHost));
    float* fdata = (float*)data;

    int nb=0, nc=0;
    for (i=0;i<nbBlocks;i++){
        if (stateArray[i]==0 || stateArray[i]==1){
            float Median = constantMedianArray[nb];
            if (Median>1) printf("data%i:%f\n",i, Median);
            for (j=0;j<blockSize;j++)
                *((*newData)+i*blockSize+j) = Median;
            nb++;
        }else if(stateArray[i]==3){
            for (j=0;j<blockSize;j++)
                *((*newData)+i*blockSize+j) = fdata[nc*blockSize+j];
            nc++;
        }
    }

	free(stateArray);
	free(constantMedianArray);
	free(data);
    cudaFree(d_newdata);
    cudaFree(d_stateArray);
    checkCudaErrors(cudaFree(d_data));

}

__device__ inline void longToBytes_bigEndian_d(unsigned char *b, unsigned long num) 
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

__device__ inline void shortToBytes_d(unsigned char* b, short value)
{
	lint16 buf;
	buf.svalue = value;
	memcpy(b, buf.byte, 2);
}

__global__ void device_post_proc(size_t *outSize, float *oriData, unsigned char *meta, short *offsets, unsigned char *midBytes, unsigned char *outBytes, size_t nbEle, int blockSize, uint64_t num_sig, uint32_t *blk_idx, float *blk_vals, uint8_t *blk_subidx, uint8_t *blk_sig)
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
    int s0 = 0;
    int s1 = 0;
    int s2 = 0;
    int s3 = 0;
    for (int i=0; i<nbBlocks; i++){
        if (meta[i]==0 || meta[i]==1 || meta[i] == 2) nbConstantBlocks++;
        else out_size += 1+(blockSize/4)+offsets[i];
    
    	if(meta[i]==0) s0++;
    	if(meta[i]==1) s1++;
    	if(meta[i]==2) s2++;
    	if(meta[i]==3) s3++;
    }
    printf("%d %d %d %d\n", s0, s1, s2, s3);
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
	//sizeToBytes(r, nbConstantBlocks);
    longToBytes_bigEndian_d(r, nbConstantBlocks);
	r += sizeof(size_t);
    //sizeToBytes(r, (size_t) num_sig);
    longToBytes_bigEndian_d(r, (unsigned long)num_sig);
    r += sizeof(size_t); 
	r += convert_state_to_out(meta, nbBlocks, r);
    r += convert_block2_to_out(r, nbBlocks,num_sig, blk_idx, blk_vals, blk_subidx, blk_sig);
    memcpy(r, oriData+nbBlocks*blockSize, (nbEle%blockSize)*sizeof(float));
    r += (nbEle%blockSize)*sizeof(float);
    unsigned char* c = r;
    unsigned char* o = c+nbConstantBlocks*sizeof(float);
    unsigned char* nc = o+(nbBlocks-nbConstantBlocks)*sizeof(short);
    for (int i=0; i<nbBlocks; i++){
        
        if (meta[i]==0 || meta[i] == 1){
            memcpy(c, meta+(nbBlocks+i*mSize), sizeof(float));
            c += sizeof(float);
        }else if(meta[i] == 3){
           shortToBytes_d(o, offsets[i]);
            o += sizeof(short);
            memcpy(nc, meta+(nbBlocks+i*mSize), mSize);
            nc += mSize; 
            memcpy(nc, midBytes+(i*blockSize*sizeof(float)), offsets[i]);
            nc += offsets[i];
        } 
    }

    // return out_size;
    *outSize = (uint32_t) (nc-r_old);
    // return (uint32_t) (nc-r_old);
}

unsigned char* device_ptr_cuSZx_compress_float(float *oriData, size_t *outSize, float absErrBound, size_t nbEle, int blockSize, float threshold)
{
    /**
     * Assuming the following are device pointers:
     *  float *oriData
     *  size_t *outSize
     *  unsigned char* outBytes
     * 
     */

    float sparsity_level = SPARSITY_LEVEL;

    // Set the input data as the function parameter, this should be a device pointer

	float* d_oriData = oriData;
    // cudaMalloc((void**)&d_oriData, sizeof(float)*nbEle); 
    // cudaMemcpy(d_oriData, oriData, sizeof(float)*nbEle, cudaMemcpyHostToDevice); 

	size_t nbBlocks = nbEle/blockSize;
	size_t remainCount = nbEle%blockSize;
	size_t actualNBBlocks = remainCount==0 ? nbBlocks : nbBlocks+1;

    size_t ncBytes = blockSize/4;
    //ncBytes = (blockSize+1)%4==0 ? ncBytes : ncBytes+1; //Bytes to store one non-constant block data.
    size_t mSize = sizeof(float)+1+ncBytes; //Number of bytes for each data block's metadata.
    size_t msz = (1+mSize) * nbBlocks * sizeof(unsigned char);
    size_t mbsz = sizeof(float) * nbEle * sizeof(unsigned char);

    // These are host pointers and do not need to be allocated

    // unsigned char *meta = (unsigned char*)malloc(msz);
    // short *offsets = (short*)malloc(nbBlocks*sizeof(short));
    // unsigned char *midBytes = (unsigned char*)malloc(mbsz);

	unsigned char* d_meta;
	unsigned char* d_midBytes;
	short* d_offsets;

    uint32_t *blk_idx, *d_blk_idx;
    uint8_t *blk_sig, *d_blk_sig;
    uint8_t *blk_subidx, *d_blk_subidx;
    float *blk_vals, *d_blk_vals;
    uint64_t *num_sig, *d_num_sig;

    checkCudaErrors(cudaMalloc((void **)&d_num_sig, sizeof(uint64_t)));
    num_sig = (uint64_t *)malloc(sizeof(uint64_t));
    checkCudaErrors(cudaMalloc((void **)&d_blk_idx, nbBlocks*sizeof(uint32_t)));
    // blk_idx = malloc()
    checkCudaErrors(cudaMalloc((void **)&d_blk_subidx, nbEle*sizeof(uint8_t)));

    checkCudaErrors(cudaMalloc((void **)&d_blk_vals, nbEle*sizeof(float)));

    checkCudaErrors(cudaMalloc((void **)&d_blk_sig, nbBlocks*sizeof(uint8_t)));

    checkCudaErrors(cudaMalloc((void**)&d_meta, msz)); 
    //checkCudaErrors(cudaMemcpy(d_meta, meta, msz, cudaMemcpyHostToDevice)); 
    checkCudaErrors(cudaMemset(d_meta, 0, msz));
    checkCudaErrors(cudaMalloc((void**)&d_offsets, nbBlocks*sizeof(short))); 
    checkCudaErrors(cudaMemset(d_offsets, 0, nbBlocks*sizeof(short)));
    checkCudaErrors(cudaMalloc((void**)&d_midBytes, mbsz)); 
    checkCudaErrors(cudaMemset(d_midBytes, 0, mbsz));

    timer_GPU.StartCounter();
    // apply_threshold<<<80,256>>>(d_oriData, threshold, nbEle);
    // cudaDeviceSynchronize();
    dim3 dimBlock(32, blockSize/32);
    dim3 dimGrid(65536, 1);
    const int sMemsize = blockSize * sizeof(float) + dimBlock.y * sizeof(int);
    compress_float<<<dimGrid, dimBlock, sMemsize>>>(d_oriData, d_meta, d_offsets, d_midBytes, absErrBound, blockSize, nbBlocks, mSize, sparsity_level, d_blk_idx, d_blk_subidx,d_blk_vals, threshold, d_blk_sig);
    cudaError_t err = cudaGetLastError();        // Get error code
    printf("CUDA Error: %s\n", cudaGetErrorString(err));
    printf("GPU compression timing: %f ms\n", timer_GPU.GetCounter());
    cudaDeviceSynchronize();
    get_numsig<<<1,1>>>(d_num_sig);
    cudaDeviceSynchronize();

    checkCudaErrors(cudaMemcpy(num_sig, d_num_sig, sizeof(uint64_t), cudaMemcpyDeviceToHost));

    // These are allocations and memcpys to host pointers, do not need them

    // blk_idx = (uint32_t *)malloc(nbBlocks*sizeof(uint32_t));
    // blk_vals= (float *)malloc((*num_sig)*sizeof(float));
    // blk_subidx = (uint8_t *)malloc((*num_sig)*sizeof(uint8_t));
    // blk_sig = (uint8_t *)malloc(nbBlocks*sizeof(uint8_t));

    // checkCudaErrors(cudaMemcpy(meta, d_meta, msz, cudaMemcpyDeviceToHost)); 
    // checkCudaErrors(cudaMemcpy(offsets, d_offsets, nbBlocks*sizeof(short), cudaMemcpyDeviceToHost)); 
    // checkCudaErrors(cudaMemcpy(midBytes, d_midBytes, mbsz, cudaMemcpyDeviceToHost)); 
    
    
    // checkCudaErrors(cudaMemcpy(blk_idx, d_blk_idx, nbBlocks*sizeof(uint32_t), cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(blk_vals,d_blk_vals, (*num_sig)*sizeof(float), cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(blk_subidx,d_blk_subidx, (*num_sig)*sizeof(uint8_t), cudaMemcpyDeviceToHost));
    // checkCudaErrors(cudaMemcpy(blk_sig,d_blk_sig, (nbBlocks)*sizeof(uint8_t), cudaMemcpyDeviceToHost));


    size_t maxPreservedBufferSize = sizeof(float)*nbEle;
    unsigned char *d_outBytes;
    // unsigned char* outBytes = (unsigned char*)malloc(maxPreservedBufferSize);
    // memset(outBytes, 0, maxPreservedBufferSize);
    checkCudaErrors(cudaMalloc(&d_outBytes, maxPreservedBufferSize));

    size_t *d_outSize;

    checkCudaErrors(cudaMalloc(&d_outSize, sizeof(size_t)));

    device_post_proc<<<1,1>>>(d_outSize, d_oriData, d_meta, d_offsets, d_midBytes, d_outBytes, nbEle, blockSize, *num_sig, d_blk_idx, d_blk_vals, d_blk_subidx, d_blk_sig);

    cudaDeviceSynchronize();
    checkCudaErrors(cudaMemcpy(outSize, d_outSize, sizeof(size_t), cudaMemcpyDeviceToHost));

    printf("completed compression\n");
    //free(blk_idx);
    //free(blk_subidx);
    //free(blk_vals);
    // free(meta);
    // free(offsets);
    // free(midBytes);
    checkCudaErrors(cudaFree(d_meta));
    checkCudaErrors(cudaFree(d_offsets));
    checkCudaErrors(cudaFree(d_midBytes));
//    printf("completed compression\n");
    return d_outBytes;
}

__device__ inline long bytesToLong_bigEndian(unsigned char* b) {
	long temp = 0;
	long res = 0;

	res <<= 8;
	temp = b[0] & 0xff;
	res |= temp;

	res <<= 8;
	temp = b[1] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[2] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[3] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[4] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[5] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[6] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[7] & 0xff;
	res |= temp;						
	
	return res;
}

__device__ inline size_t bytesToSize(unsigned char* bytes)
{
	size_t result = bytesToLong_bigEndian(bytes);//8	
	return result;
}

__global__ void decompress_startup(float **newData, size_t nbEle, unsigned char* cmpBytes, 
    uint32_t *blk_idx, uint8_t *blk_subidx, uint8_t *blk_sig,
    float *blk_vals, size_t *numSigValues, int *bs,
    size_t *numConstantBlks, size_t *numBlks, size_t *ncBlks,
    unsigned char *stateArray, float* constantMedianArray, unsigned char *data,
    size_t *mSizeptr
){
    /**
     * Structures to return:
     * blk_idx, blk_subidx, blk_sig, blk_vals, numSigValues (pointer)
     * bs (pointer to blockSize), numConstantBlks (pointer), numBlks (pointer)
     * ncBlks (pointer), stateArray, constantMedianArray
     */
	
	unsigned char* r = cmpBytes;
    size_t num_sig;
	r += 4;
	int blockSize = r[0];  //get block size
	if(blockSize == 0)blockSize = 256;
	r++;
	size_t nbConstantBlocks = bytesToLong_bigEndian(r); //get number of constant blocks
	r += sizeof(size_t);
	num_sig = bytesToSize(r);
    r += sizeof(size_t);
	size_t nbBlocks = nbEle/blockSize;
    size_t ncBlocks = 0;
    size_t num_state2_blks = 0;
	// size_t ncBlocks = nbBlocks - nbConstantBlocks; //get number of constant blocks
	size_t stateNBBytes = nbBlocks%4==0 ? nbBlocks/4 : nbBlocks/4+1;
    size_t ncLeading = blockSize/4;
    size_t mSize = sizeof(float)+1+ncLeading; //Number of bytes for each data block's metadata.
	*mSizeptr = mSize;
    unsigned char* stateArray = (unsigned char*)malloc(nbBlocks);
    // unsigned char* d_stateArray;
    // cudaMalloc(&d_stateArray, nbBlocks);
	float* constantMedianArray = (float*)malloc(nbConstantBlocks*sizeof(float));			

    blk_idx = (uint32_t *)malloc(nbBlocks*sizeof(uint32_t));
    blk_vals= (float *)malloc((num_sig)*sizeof(float));
    blk_subidx = (uint8_t *)malloc((num_sig)*sizeof(uint8_t));
    blk_sig = (uint8_t *)malloc(nbBlocks*sizeof(uint8_t));

	printf("Converting state array\n");
    convert_out_to_state(nbBlocks, r, stateArray);
	// convertByteArray2IntArray_fast_1b_args(nbBlocks, r, stateNBBytes, stateArray); //get the stateArray
	for (size_t i = 0; i < nbBlocks; i++)
    {
        if (stateArray[i] == 2)
        {
            num_state2_blks++;
        }else if(stateArray[i] == 3){
            ncBlocks++;
        }
    }
    
	r += stateNBBytes;
    data = (unsigned char*)malloc(ncBlocks*blockSize*sizeof(float));
    memset(data, 0, ncBlocks*blockSize*sizeof(float));
    printf("converting block vals\n");
    size_t to_add = convert_out_to_block2(r, nbBlocks, (uint64_t)num_sig, blk_idx, blk_vals, blk_subidx, blk_sig);
    r+= to_add;

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

    *numConstantBlks = nbConstantBlocks;
    *numBlks = nbBlocks;
    *ncBlks = ncBlocks;
    *numSigValues = num_sig;
    *bs = blockSize;
}

__global__ void decompress_post_proc(unsigned char *data, float **newData, int blockSize, 
    size_t nbBlocks, size_t ncBlocks, unsigned char *stateArray,
    float *constantMedianArray
){
    // checkCudaErrors(cudaMemcpy(data, d_data, ncBlocks*blockSize*sizeof(float), cudaMemcpyDeviceToHost)); 
    // checkCudaErrors(cudaMemcpy(*newData, d_newdata, nbBlocks*blockSize*sizeof(float), cudaMemcpyDeviceToHost));
    float* fdata = (float*)data;
    int i,j;
    int nb=0, nc=0;
    for (i=0;i<nbBlocks;i++){
        if (stateArray[i]==0 || stateArray[i]==1){
            float Median = constantMedianArray[nb];
            if (Median>1) printf("data%i:%f\n",i, Median);
            for (j=0;j<blockSize;j++)
                *((*newData)+i*blockSize+j) = Median;
            nb++;
        }else if(stateArray[i]==3){
            for (j=0;j<blockSize;j++)
                *((*newData)+i*blockSize+j) = fdata[nc*blockSize+j];
            nc++;
        }
    }
}

void device_ptr_cuSZx_decompress_float(float** newData, size_t nbEle, unsigned char* cmpBytes)
{
    /**
     * Assume the following are device pointers
     * 
     * unsigned char* cmpBytes
     * float** newData (must be already allocated on device with size = sizeof(float)*nbEle)
     * 
     */
    
    uint32_t *blk_idx;
    uint8_t *blk_subidx;
    uint8_t *blk_sig;
    float *blk_vals, *constantMedianArray;
    size_t *num_sig, *mSize, mSize_h;
    int *blockSize, bs;
    size_t *nbConstantBlocks, *nbBlocks, *ncBlocks, nbBlocks_h, ncBlocks_h;
    unsigned char *stateArray, *data;

    checkCudaErrors(cudaMalloc((void**)&num_sig, sizeof(size_t)));
    checkCudaErrors(cudaMalloc((void**)&blockSize, sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&nbConstantBlocks, sizeof(size_t)));
    checkCudaErrors(cudaMalloc((void**)&nbBlocks, sizeof(size_t)));
    checkCudaErrors(cudaMalloc((void**)&ncBlocks, sizeof(size_t)));
    

    decompress_startup<<<1,1>>>(newData, nbEle, cmpBytes, 
    blk_idx, blk_subidx, blk_sig,
    blk_vals, num_sig, blockSize,
    nbConstantBlocks, nbBlocks, ncBlocks,
    stateArray, constantMedianArray, data, mSize);
    cudaDeviceSynchronize();

    checkCudaErrors(cudaMemcpy(&nbBlocks_h, nbBlocks, sizeof(size_t), cudaMemcpyDeviceToHost)); 
    checkCudaErrors(cudaMemcpy(&ncBlocks_h, ncBlocks, sizeof(size_t), cudaMemcpyDeviceToHost)); 
    checkCudaErrors(cudaMemcpy(&bs, blockSize, sizeof(int), cudaMemcpyDeviceToHost)); 
    checkCudaErrors(cudaMemcpy(&mSize_h, mSize, sizeof(size_t), cudaMemcpyDeviceToHost)); 

    // unsigned char* d_data;
    float *d_newdata;
    // checkCudaErrors(cudaMalloc((void**)&d_data, ncBlocks*blockSize*sizeof(float))); 
    // checkCudaErrors(cudaMemcpy(d_data, data, ncBlocks*blockSize*sizeof(float), cudaMemcpyHostToDevice)); 
    checkCudaErrors(cudaMalloc(&d_newdata, nbBlocks*blockSize*sizeof(float)));

    timer_GPU.StartCounter();
    dim3 dimBlock(32, bs/32);
    dim3 dimGrid(65536, 1);
    const int sMemsize = bs * sizeof(float) + dimBlock.y * sizeof(int);
    decompress_state2<<<nbBlocks_h, 64>>>(d_newdata, stateArray,blk_idx, blk_vals, blk_subidx, bs, blk_sig);
    decompress_float<<<dimGrid, dimBlock, sMemsize>>>(data, bs, ncBlocks_h, mSize_h);
    cudaError_t err = cudaGetLastError();        // Get error code
    printf("CUDA Error: %s\n", cudaGetErrorString(err));
    printf("GPU decompression timing: %f ms\n", timer_GPU.GetCounter());
    cudaDeviceSynchronize();
    // checkCudaErrors(cudaMemcpy(data, d_data, ncBlocks*blockSize*sizeof(float), cudaMemcpyDeviceToHost)); 
    checkCudaErrors(cudaMemcpy(*newData, d_newdata, nbBlocks*blockSize*sizeof(float), cudaMemcpyDeviceToDevice));
    cudaFree(d_newdata);

    decompress_post_proc<<<1,1>>>(data, newData, bs, 
    nbBlocks_h, ncBlocks_h, stateArray,
    constantMedianArray);
    cudaDeviceSynchronize();

	cudaFree(stateArray);
	cudaFree(constantMedianArray);
	cudaFree(data);
    cudaFree(blk_idx);
    cudaFree(blk_subidx);
    cudaFree(blk_vals);
    cudaFree(blk_sig);

}

