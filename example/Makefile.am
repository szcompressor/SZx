AUTOMAKE_OPTIONS=foreign
AM_CFLAGS = -I../sz/include -I../zlib -I../zstd $(GSL_CFLAGS) $(GSL_STATIC)
AM_FCFLAGS = -I../sz/include -I../sz -I../zlib -I../zstd $(GSL_CFLAGS) $(GSL_STATIC)


bin_PROGRAMS=testfloat_compress testfloat_decompress testfloat_compress2 testfloat_decompress2

testfloat_compress_SOURCES=testfloat_compress.c
testfloat_compress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testfloat_decompress_SOURCES=testfloat_decompress.c
testfloat_decompress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testfloat_compress2_SOURCES=testfloat_compress2.c
testfloat_compress2_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
testfloat_decompress2_SOURCES=testfloat_decompress2.c
testfloat_decompress2_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a ../zstd/.libs/libzstd.a -lm
