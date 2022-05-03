SZx: An Ultrafast Error-Bounded Lossy Compressor for Scientific Datasets
=====
 (C) 2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
       See COPYRIGHT in top-level directory.

* Major Authors: Sheng Di, Xiaodong Yu, Kai Zhao 
* Supervisor: Franck Cappello

## Citations
* SZx: Xiaodong Yu, Sheng Di, Kai Zhao, Jiannan Tian, Dingwen Tao, Xin Liang, and Franck Cappello. "[Ultrafast Error-Bounded Lossy Compression for Scientific Datasets](https:)", In proceedings of the 31st International ACM Symposium on High-Performance Parallel and Distributed Computing (HPDC '22), Minneapolis, Minnesota, USA, 2022.

SZx is a novel ultrafast error-bounded lossy compressor that can obtain fairly high compression performance on both CPUs and GPUs and with reasonably high compression ratios. It composes only lightweight operations such as bitwise operations, additions, and subtractions. It also can support strict control of the compression errors within user-specified error bounds. For more details please refer to the HPDC'22 paper. This document simply introduces how to install and use the SZx compressor.  

## Installation

### Prerequisites
cmake >= 3.19 <br />
CUDA toolkit >= 10.0 with samples/common/inc

## How to install
Under the SZx directory:
```bash
mkdir build
cd build
cmake ../
make -j12
```

After the installation, you'll find all the executables in build/bin.

## Testing Examples
--------------------------------------

Please refer to ```testfloat_compress_fastmode2 --help``` and ```testfloat_decompress_fastmode2 --help``` for more details.

## Compression
--------------
```bash
testfloat_compress_fastmode2 testfloat_8_8_128.dat 64 1E-3 --cuda
```

`Description: `

testfloat_8_8_128.dat is the binary testing file (small-endian format), which contains a 3d array (128X8X8). Its data values are shown in the plain text file testdouble_8_8_128.txt. It comes from FLASH_Blast2  at time step 100. 

64 is the size of each data block. For more information about data blocks, please read our paper.

1E-3 is the user-specific error-bound. Currently, SZx only supports absolute error-bound.

--cuda is an optional argument. Use it to enable the GPU-based compression (i.e., cuSZx).

The compressed data files will have a suffix `.szx`. For example, a compressed file named testdfloat_8_8_128.dat.szx will be generated after this compression.

## Decompression

```bash
testfloat_decompress_fastmode2 testfloat_8_8_128.dat.szx 8192 --cuda
```

`Description: `

testfloat_8_8_128.dat.szx is the compressed binary file. 

8192 is the number of elements in the original data. The decompressed data should have the same number of elements.

--cuda is an optional argument. Use it to enable the GPU-based decompression (i.e., cuSZx).

The decompressed data files will have a suffix `.szx.out`. For example, a decompressed file named testdfloat_8_8_128.dat.szx.out will be generated after this decompression.

