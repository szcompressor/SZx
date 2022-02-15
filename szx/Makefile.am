#AM_CFLAGS = -I./include -I../zlib
#LDFLAGS=-fPIC -shared

AUTOMAKE_OPTIONS=foreign
include_HEADERS=include/szx.h include/szx_rw.h include/szx_defines.h include/szx_float.h include/szx_BytesToolkit.h include/szx_TypeManager.h include/szx_dataCompression.h include/szxd_float.h include/szx_utility.h

lib_LTLIBRARIES=libSZx.la
libSZx_la_CFLAGS=-I./include
if OPENMP
libSZx_la_CFLAGS+=-fopenmp
endif
libSZx_la_LDFLAGS = -version-info  1:0:0
libSZx_la_SOURCES=src/szx.c src/szx_float.c src/szx_rw.c src/szx_BytesToolkit.c src/szx_TypeManager.c src/szx_dataCompression.c src/szxd_float.c src/szx_utility.c

libSZ_la_LINK= $(AM_V_CC)$(LIBTOOL) --tag=CC --mode=link $(CCLD) $(libSZx_la_CFLAGS) -O3 $(libSZx_la_LDFLAGS) -o $(lib_LTLIBRARIES)
