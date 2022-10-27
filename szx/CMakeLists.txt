file(GLOB source_files "src/*.c")
if (CUDAToolkit_FOUND)
    file(GLOB cuda_source_files "src/*.cu")
    list(APPEND source_files ${cuda_source_files})
endif ()

add_library(SZx ${source_files})

target_include_directories(SZx
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/szx>
        )

target_link_libraries(SZx PUBLIC ${ZLIB_dep} ${ZSTD_dep} m)
if (OpenMP_C_FOUND)
    target_link_libraries(SZx PUBLIC OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif ()

target_compile_options(SZx
        PRIVATE $<$<CONFIG:Debug>:-Wall -g -Wextra -Wpedantic -Wno-unused-parameter>
        )
#target_compile_options(SZx
#        PRIVATE -march=native
#        )

install(TARGETS SZx ${thirdparty_export} EXPORT SZxConfig
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )

install(EXPORT SZxConfig DESTINATION share/SZx/cmake)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/szx)
export(TARGETS SZx FILE SZx.cmake)
