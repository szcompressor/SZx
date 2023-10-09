AUTOMAKE_OPTIONS=foreign
AM_CPPFLAGS = -I../szx/include
AM_FCFLAGS = -I../szx/include -I../szx
if OPENMP
AM_CPPFLAGS+=-fopenmp
endif

bin_PROGRAMS=testfloat_compress_fastmode1 testfloat_decompress_fastmode1 testfloat_compress_fastmode3 testfloat_decompress_fastmode3 szx testfloat_compress_fixed_rate estimateCR

testfloat_compress_fastmode1_SOURCES=testfloat_compress_fastmode1.cc
testfloat_compress_fastmode1_LDADD=../szx/.libs/libSZx.a -lm
testfloat_decompress_fastmode1_SOURCES=testfloat_decompress_fastmode1.cc
testfloat_decompress_fastmode1_LDADD=../szx/.libs/libSZx.a -lm
testfloat_compress_fastmode3_SOURCES=testfloat_compress_fastmode3.cc
testfloat_compress_fastmode3_LDADD=../szx/.libs/libSZx.a -lm
testfloat_decompress_fastmode3_SOURCES=testfloat_decompress_fastmode3.cc
testfloat_decompress_fastmode3_LDADD=../szx/.libs/libSZx.a -lm
szx_SOURCES=szx.cc
szx_LDADD=../szx/.libs/libSZx.a -lm
testfloat_compress_fixed_rate_SOURCES=testfloat_compress_fixed_rate.cc
testfloat_compress_fixed_rate_LDADD=../szx/.libs/libSZx.a -lm
estimateCR_SOURCES=estimateCR.cc
estimateCR_LDADD=../szx/.libs/libSZx.a -lm
