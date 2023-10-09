#AM_CFLAGS = -I./include -I../zlib
#LDFLAGS=-fPIC -shared

AUTOMAKE_OPTIONS=foreign
include_HEADERS=include/szx.h include/szx_rw.h include/szx_defines.h include/szx_float.h include/szx_BytesToolkit.h include/szx_TypeManager.h include/szx_dataCompression.h include/szxd_float.h include/szx_utility.h include/szx_double.h include/szxd_double.h include/szx_globals.h

lib_LTLIBRARIES=libSZx.la
libSZx_la_CXXFLAGS=-I./include
if OPENMP
libSZx_la_CXXFLAGS+=-fopenmp
endif
libSZx_la_LDFLAGS = -version-info  1:0:0
libSZx_la_SOURCES=src/szx_BytesToolkit.cc src/szx.cc src/szx_float.cc src/szx_rw.cc  src/szx_TypeManager.cc src/szx_dataCompression.cc src/szxd_float.cc src/szx_utility.cc src/szx_double.cc src/szxd_double.cc

libSZx_la_LINK= $(AM_V_CXX)$(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(libSZx_la_CXXFLAGS) -O3 $(libSZx_la_LDFLAGS) -o $(lib_LTLIBRARIES)
