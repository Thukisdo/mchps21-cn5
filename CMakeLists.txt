cmake_minimum_required(VERSION 3.0)

project(CN-TP5)

# Add the path to lapacke.h
if (NOT DEFINED LAPACK_DIR)
  message(FATAL_ERROR "LAPACK_DIR must be set prior to launching CMake")
endif ()

set(BLA_VENDOR Intel10_64lp_seq)
#set(BLA_VENDOR OpenBLAS)

if (CMAKE_C_COMPILER_ID MATCHES "INTEL")
  set(CMAKE_C_FLAGS_RELEASE  "${CMAKE_C_FLAGS_RELEASE} -qmkl=sequential -xHost -fast")
else ()
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
  if (COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -march=native")
  endif ()
endif ()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
add_subdirectory(src)

# Stores scripts for dependency resolution
# and installation
add_subdirectory(extern)

# Enable testing in the root directory
enable_testing()
add_subdirectory(tests)

# Tools and utilities
add_subdirectory(tools)