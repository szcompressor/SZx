add_library (SZx STATIC
	src/szx.c
	src/szx_float.c
	src/szx_rw.c
	src/szx_BytesToolkit.c
	src/szx_TypeManager.c
	src/szx_dataCompression.c
	src/szxd_float.c
	src/szx_utility.c
	src/szx_double.c
	src/szxd_double.c	
)

target_include_directories(SZx
  PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/szx>
  )

target_link_libraries (SZx PUBLIC ${ZLIB_dep} ${ZSTD_dep} m)
if(OpenMP_C_FOUND)
    target_link_libraries(SZx PUBLIC OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()

target_compile_options(SZx
	PRIVATE $<$<CONFIG:Debug>:-Wall -g -Wextra -Wpedantic -Wno-unused-parameter>
	)
target_compile_options(SZx
        PRIVATE -march=native
        )

install (TARGETS SZx ${thirdparty_export} EXPORT SZxConfig
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

install (EXPORT SZxConfig DESTINATION share/SZx/cmake)
install (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/szx)
export(TARGETS SZx FILE SZx.cmake)
