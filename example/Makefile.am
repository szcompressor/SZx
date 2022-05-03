AUTOMAKE_OPTIONS=foreign
AM_CFLAGS = -I../szx/include
AM_FCFLAGS = -I../szx/include -I../szx
if OPENMP
AM_CFLAGS+=-fopenmp
endif


bin_PROGRAMS=testfloat_compress_fastmode1 testfloat_decompress_fastmode1 testfloat_compress_fastmode3 testfloat_decompress_fastmode3 szx

testfloat_compress_fastmode1_SOURCES=testfloat_compress_fastmode1.c
testfloat_compress_fastmode1_LDADD=../szx/.libs/libSZx.a -lm
testfloat_decompress_fastmode1_SOURCES=testfloat_decompress_fastmode1.c
testfloat_decompress_fastmode1_LDADD=../szx/.libs/libSZx.a -lm
testfloat_compress_fastmode3_SOURCES=testfloat_compress_fastmode3.c
testfloat_compress_fastmode3_LDADD=../szx/.libs/libSZx.a -lm
testfloat_decompress_fastmode3_SOURCES=testfloat_decompress_fastmode3.c
testfloat_decompress_fastmode3_LDADD=../szx/.libs/libSZx.a -lm
szx_SOURCES=szx.c
szx_LDADD=../szx/.libs/libSZx.a -lm
