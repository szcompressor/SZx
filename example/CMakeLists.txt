add_executable (testint_compress testint_compress.c)
target_link_libraries (testint_compress SZ)

add_executable (testint_decompress testint_decompress.c)
target_link_libraries (testint_decompress SZ)

add_executable (testfloat_compress testfloat_compress.c)
target_link_libraries (testfloat_compress SZ)

add_executable (testfloat_decompress testfloat_decompress.c)
target_link_libraries (testfloat_decompress SZ)

add_executable (testdouble_compress testdouble_compress.c)
target_link_libraries (testdouble_compress SZ)

add_executable (testdouble_decompress testdouble_decompress.c)
target_link_libraries (testdouble_decompress SZ)

add_executable (testfloat_compress_fastmode1 testfloat_compress_fastmode1.c)
target_link_libraries (testfloat_compress_fastmode1 SZ)

add_executable (testfloat_decompress_fastmode1 testfloat_decompress_fastmode1.c)
target_link_libraries (testfloat_decompress_fastmode1 SZ)

#add_executable (testfloat_compress_fastmode2 testfloat_compress_fastmode2.c)
#target_link_libraries (testfloat_compress_fastmode2 SZ)

add_executable (testfloat_compress_fastmode2 testfloat_compress_fastmode2.cu)
set_property(TARGET testfloat_compress_fastmode2 PROPERTY CUDA_ARCHITECTURES 70 75)
set_target_properties (testfloat_compress_fastmode2
                    PROPERTIES
                    CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries (testfloat_compress_fastmode2 SZ)

#add_executable (testfloat_decompress_fastmode2 testfloat_decompress_fastmode2.c)
add_executable (testfloat_decompress_fastmode2 testfloat_decompress_fastmode2.cu)
set_property(TARGET testfloat_decompress_fastmode2 PROPERTY CUDA_ARCHITECTURES 70 75)
set_target_properties (testfloat_decompress_fastmode2
                    PROPERTIES
                    CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries (testfloat_decompress_fastmode2 SZ)

add_executable (sz sz.c)
target_link_libraries (sz SZ)

install (TARGETS testint_compress testint_decompress testfloat_compress testfloat_decompress testfloat_compress_fastmode1 testfloat_compress_fastmode2
	testfloat_decompress_fastmode1 testfloat_decompress_fastmode2 testdouble_compress testdouble_decompress sz
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})


