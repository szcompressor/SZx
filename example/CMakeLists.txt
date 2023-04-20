#include MPI
find_package(MPI REQUIRED)
set(SZX_TARGETS
		szx
		testfloat_compress_fastmode1
		testfloat_compress_fastmode3
		testfloat_decompress_fastmode1
		testfloat_decompress_fastmode3
		ANL_allreduce
		main
		test_multicomm
		SZX_test
		compare_data
		)
foreach (SZX_TARGET IN LISTS SZX_TARGETS)
	add_executable (${SZX_TARGET} ${SZX_TARGET}.c)
	target_link_libraries (${SZX_TARGET} SZx MPI::MPI_C)
endforeach ()

if (CMAKE_CUDA_COMPILER)
	set(SZX_TARGETS ${SZX_TARGETS}
			testfloat_compress_fastmode2
			testfloat_decompress_fastmode2
			)
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
endif()

install (TARGETS ${SZX_TARGETS}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})