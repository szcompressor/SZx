AUTOMAKE_OPTIONS=foreign
AM_CFLAGS = -I../sz/include -I../zlib -I../zstd $(GSL_CFLAGS) $(GSL_STATIC)
AM_FCFLAGS = -I../sz/include -I../sz -I../zlib -I../zstd $(GSL_CFLAGS) $(GSL_STATIC)


bin_PROGRAMS=sz testfloat_compress testfloat_decompress testfloat_compress_fastmode1 testfloat_compress_fastmode2 testfloat_decompress_fastmode1 testfloat_decompress_fastmode2

sz_SOURCES=sz.c
sz_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testfloat_compress_SOURCES=testfloat_compress.c
testfloat_compress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testfloat_decompress_SOURCES=testfloat_decompress.c
testfloat_decompress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testfloat_compress_fastmode1_SOURCES=testfloat_compress_fastmode1.c
testfloat_compress_fastmode1_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testfloat_compress_fastmode2_SOURCES=testfloat_compress_fastmode2.c
testfloat_compress_fastmode2_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testfloat_decompress_fastmode1_SOURCES=testfloat_decompress_fastmode1.c
testfloat_decompress_fastmode1_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testfloat_decompress_fastmode2_SOURCES=testfloat_decompress_fastmode2.c
testfloat_decompress_fastmode2_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
