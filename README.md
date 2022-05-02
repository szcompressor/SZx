SZx: An Ultrafast Error-Bounded Lossy Compressor for Scientific Datasets
=====
 (C) 2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
       See COPYRIGHT in top-level directory.

* Major Authors: Sheng Di, Xiaodong Yu, Kai Zhao 
* Supervisor: Franck Cappello

## Citations
* SZx: Xiaodong Yu, Sheng Di, Kai Zhao, Jiannan Tian, Dingwen Tao, Xin Liang, and Franck Cappello. "[Ultrafast Error-Bounded Lossy Compression for Scientific Datasets](https:)", In proceedings of the 31st International ACM Symposium on High-Performance Parallel and Distributed Computing (HPDC '22), Minneapolis, Minnesota, USA, 2022.

This document simply introduces how to install and use the SZx compressor.  

## Installation

### Installation way 1:
* ./configure --prefix=[INSTALL_DIR] (Please use --enable-fortran if you need Fortran interface)
* make
* make install

### Installation way 2:
* mkdir build && cd build
* cmake .. -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR]
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and .a and .so libraries in [INSTALL_DIR]/lib

## Testing Examples
--------------------------------------

Examples can be found in the [SZ_PACKAGE]/example

You can use the executable 'sz' command to do the compression/decompression. Please see the user guide or run 'sz --help' for details.

Alternatively, you can also also call our API to do the compression/decompressoin. Here are two examples: testfloat_compress.c and testfloat_decompress.c

## Compression
--------------
* ./test_compress sz.config testdouble_8_8_8_128.dat 8 8 8 128
* ./test_compress sz.config testdouble_8_8_128.dat 8 8 128

`Decription: `

testdouble_8_8_128.dat and testdouble_8_8_8_128.dat are two binary testing files (small-endian format), which contains a 3d array (128X8X8) and a 4d array (128X8X8X8) respectively. Their data values are shown in the two plain text files, testdouble_8_8_128.txt and testdouble_8_8_8_128.txt. These two data files are from FLASH_Blast2 and FLASH_MacLaurin respectively (both are at time step 100). The compressed data files are namely testdouble_8_8_8_128.dat.sz and testdouble_8_8_128.dat.sz respectively.

sz.config is the configuration file. The key settings are errorBoundMode, absErrBound, and relBoundRatio, which are described below.

* absErrBound refers to the absolute error bound, which is to limit the (de)compression errors to be within an absolute error. For example, absErrBound=0.0001 means the decompressed value must be in [V-0.0001,V+0.0001], where V is the original true value.

* relBoundRatio refers to relative bound ratio, which is to limit the (de)compression errors by considering the global data value range size (i.e., taking into account the range size (max_value - min_value)). For example, suppose relBoundRatio is set to 0.01, and the data set is {100,101,102,103,104,...,110}, so the global value range size is 110-100=10, so the error bound will actually be 10*0.01=0.1, from the perspective of "relBoundRatio".

* errorBoundMode is to indicate the error-bounding way in the compression, such as based on absolute error bound, relative error bound, etc. 
The options are shown below.
	* ABS: take only "absolute error bound" into account. That is, relative bound ratio will be ignored.
	* REL: take only "relative bound ratio" into account. That is, absolute error bound will be ignored. 
	* ABS_AND_REL: take both of the two bounds into account. The compression errors will be limited using both absErrBound and relBoundRatio*rangesize. That is, the two bounds must be both met.
	* ABS_OR_REL: take both of the two bounds into account. The compression errors will be limited using either absErrBound or relBoundRatio*rangesize. That is, only one bound is required to be met.
	* PW_REL: take "point-wise relative error bound" in the compression. 

## Decompression

* ./test_decompress testdouble_8_8_8_128.dat.sz 8 8 8 128
* ./test_decompress testdouble_8_8_128.dat.sz 8 8 128

The output files are testdouble_8_8_8_128.dat.sz.out and testdouble_8_8_128.dat.sz.out respectively. You can compare .txt file and .out file for checking the compression errors for each data point. For instance, compare testdouble_8_8_8_128.txt and testdouble_8_8_8_128.dat.sz.out.

