cmake_minimum_required (VERSION 3.10.2)

include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
endif()

#disable in-source builds
set (CMAKE_DISABLE_SOURCE_CHANGES ON)
set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)

#define the project
project (szx
  VERSION 1.0.0
  DESCRIPTION "SZx Error Bounded Lossy Compressor"
  LANGUAGES C CXX CUDA
  )
enable_testing()

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

#set the Compiler ID for clang on macOS to AppleClang
if (POLICY CMP0025)
  cmake_policy (SET CMP0025 NEW)
endif()

#compile with C-99 and standard C++14
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-Wall")

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/szx/include)
find_package(CUDAToolkit  REQUIRED)
include_directories("${CUDAToolkit_INCLUDE_DIRS}")
include_directories (${CUDAToolkit_LIBRARY_ROOT}/samples/common/inc)
#include_directories (/usr/local/cuda/samples/common/inc)

if (NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED True)
endif ()
#generate tags for the project if tags exist
option(BUILD_CTAGS "enable ctags generation target" OFF)
if(BUILD_CTAGS)
  find_program(TAGS ctags)
  if(TAGS)
    add_custom_target(tags ALL
      COMMAND ${TAGS} --exclude=${CMAKE_BINARY_DIR} -f ${CMAKE_BINARY_DIR}/tags --c++-kinds=+p --fields=+iaS -R
      COMMENT Generating Tag files
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      )
  endif()
endif()

option(BUILD_SHARED_LIBS "build shared libraries over static libraries" ON)

#find dependencies
include(GNUInstallDirs)
find_package(PkgConfig)
pkg_search_module(ZSTD  IMPORTED_TARGET libzstd)

#by default pass no 3rd party exports
set(thirdparty_export "")

find_package(OpenMP)

add_subdirectory(szx)
add_subdirectory(example)

option(BUILD_PYTHON_WRAPPER "build python wrapper" OFF)
if(BUILD_PYTHON_WRAPPER)
  add_subdirectory(swig)
endif()

if(BUILD_DOCKER_CONTAINERS)
  foreach(CONTAINER Centos Fedora Ubuntu Travis CentosPackaged)
    set(BuildSentinel ${CMAKE_BINARY_DIR}/${CONTAINER}-built)
    set(Dockerfile docker/Dockerfile-${CONTAINER})
    string(TOLOWER "szx${CONTAINER}" CONTAINER_TAG)
    add_custom_command(OUTPUT ${BuildSentinel}
      COMMAND sudo docker build -t ${CONTAINER_TAG} -f ${Dockerfile} .
      COMMAND touch ${BuildSentinel}
      MAIN_DEPENDENCY ${Dockerfile}
      DEPENDS SZx
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      COMMENT "DOCKER ${Dockerfile}"
      )
    list(APPEND DOCKER_CONTAINERS ${BuildSentinel})
  endforeach()
  add_custom_target(docker DEPENDS ${DOCKER_CONTAINERS} COMMENT "building docker containers")
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/szx.pc.in
  ${CMAKE_BINARY_DIR}/szx.pc
  @ONLY
)
install(FILES ${CMAKE_BINARY_DIR}/szx.pc DESTINATION ${CMAKE_INSTALL_PREFIX}/share/pkgconfig)

