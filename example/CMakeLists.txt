add_executable (testfloat_compress_fastmode1 testfloat_compress_fastmode1.c)
target_link_libraries (testfloat_compress_fastmode1 SZx)

add_executable (testfloat_decompress_fastmode1 testfloat_decompress_fastmode1.c)
target_link_libraries (testfloat_decompress_fastmode1 SZx)

add_executable (testfloat_compress_fastmode2 testfloat_compress_fastmode2.cu)
set_property(TARGET testfloat_compress_fastmode2 PROPERTY CUDA_ARCHITECTURES 70 75)
set_target_properties (testfloat_compress_fastmode2
                    PROPERTIES
                    CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries (testfloat_compress_fastmode2 SZx)

add_executable (testfloat_decompress_fastmode2 testfloat_decompress_fastmode2.cu)
set_property(TARGET testfloat_decompress_fastmode2 PROPERTY CUDA_ARCHITECTURES 70 75)
set_target_properties (testfloat_decompress_fastmode2
                    PROPERTIES
                    CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries (testfloat_decompress_fastmode2 SZx)

add_executable (testfloat_compress_fastmode3 testfloat_compress_fastmode3.c)
target_link_libraries (testfloat_compress_fastmode3 SZx)

add_executable (testfloat_decompress_fastmode3 testfloat_decompress_fastmode3.c)
target_link_libraries (testfloat_decompress_fastmode3 SZx)

add_executable (szx szx.c)
target_link_libraries (szx SZx)

install (TARGETS testfloat_compress_fastmode1 testfloat_compress_fastmode2 testfloat_compress_fastmode3
	testfloat_decompress_fastmode1 testfloat_decompress_fastmode2 testfloat_decompress_fastmode3 szx
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})


