option(SZx_INSTALL_CLI "install SZx CLI" ON)
option(SZx_INSTALL_EXAMPLES "install SZx examples" OFF)

add_executable (szx szx.cc)
target_link_libraries(szx PRIVATE SZx)
if(SZx_INSTALL_CLI)
install(TARGETS szx)
endif()

set(SZX_TARGETS
		testfloat_compress_fastmode1
		testfloat_compress_fastmode3
		testfloat_decompress_fastmode1
		testfloat_decompress_fastmode3
		estimateCR
		)

foreach (SZX_TARGET IN LISTS SZX_TARGETS)
	add_executable (szx_${SZX_TARGET} ${SZX_TARGET}.cc)
    target_link_libraries (szx_${SZX_TARGET} PRIVATE SZx)
    if(SZx_INSTALL_EXAMPLES)
        install(TARGETS szx_${SZX_TARGET})
    endif()
endforeach ()

if (SZx_BUILD_CUDA)
	add_executable (szx_testfloat_compress_fastmode2 testfloat_compress_fastmode2.cu)
    target_link_libraries (szx_testfloat_compress_fastmode2 PRIVATE SZx)

	add_executable (szx_testfloat_decompress_fastmode2 testfloat_decompress_fastmode2.cu)
    target_link_libraries (szx_testfloat_decompress_fastmode2 PRIVATE SZx)

    if(SZx_INSTALL_EXAMPLES)
        install(TARGETS szx_testfloat_compress_fastmode2)
        install(TARGETS szx_testfloat_decompress_fastmode2)
    endif()
endif()

