cmake_minimum_required (VERSION 3.20)
#define the project
project (szx
        VERSION 1.1.1
        DESCRIPTION "SZx Error Bounded Lossy Compressor"
        )


set (CMAKE_EXPORT_COMPILE_COMMANDS ON)
#disable in-source builds
set (CMAKE_DISABLE_SOURCE_CHANGES ON)
set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)
#set the Compiler ID for clang on macOS to AppleClang
if (POLICY CMP0025)
  cmake_policy (SET CMP0025 NEW)
endif()

include(CTest)
include(GNUInstallDirs)
find_package(OpenMP)

option(BUILD_SHARED_LIBS "build shared libraries over static libraries" ON)
option(SZx_BUILD_CUDA "build the cuda version of SZx" OFF)

#correct was to set a default build type
# https://blog.kitware.com/cmake-and-the-default-build-type/
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type was set. Setting build type to ${default_build_type}.")
  set(CMAKE_BUILD_TYPE ${default_build_type} CACHE
    STRING "Choose the type to build" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

add_library(SZx
    szx/src/szx.cc
    szx/src/szx_BytesToolkit.cc
    szx/src/szx_TypeManager.cc
    szx/src/szx_dataCompression.cc
    szx/src/szx_double.cc
    szx/src/szx_float.cc
    szx/src/szx_rw.cc
    szx/src/szx_utility.cc
    szx/src/szxd_double.cc
    szx/src/szxd_float.cc
    )
target_include_directories(SZx PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/szx/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/szx>
    )
target_compile_features(SZx PUBLIC cxx_std_14 c_std_99)
target_link_libraries(SZx PUBLIC OpenMP::OpenMP_CXX)

if(SZx_BUILD_CUDA)
    enable_language(CUDA)
    target_sources(SZx PRIVATE
        szx/src/cuszx_entry.cu
        szx/src/cuszx_float.cu
        szx/src/cuszxd_float.cu
        szx/src/timingGPU.cu
        )
    target_include_directories(SZx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/szx/src/)
endif()

add_subdirectory(example)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/szx.pc.in
  ${CMAKE_BINARY_DIR}/szx.pc
  @ONLY
)
install(FILES ${CMAKE_BINARY_DIR}/szx.pc DESTINATION ${CMAKE_INSTALL_PREFIX}/share/pkgconfig)

install(TARGETS SZx EXPORT SZxTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
install(DIRECTORY szx/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/szx)
install(EXPORT SZxTargets NAMESPACE SZx:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZx/)
export(EXPORT SZxTargets FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/SZxTargets.cmake NAMESPACE SZx::)
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/SZxConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/SZxConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZx
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/SZxConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY AnyNewerVersion
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SZxConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SZxConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZx
)

