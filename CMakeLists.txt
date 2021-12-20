cmake_minimum_required(VERSION 3.0)

project(CN-TP5)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/blas.cmake)

try_find_blas(OpenBLAS)

if (NOT DEFINED LAPACKE_DIR)
  set(LAPACKE_DIR /usr/include/lapacke/)
endif()

if (CMAKE_C_COMPILER_ID MATCHES "INTEL")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -qmkl=sequential -xHost -fast")
else ()
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
  if (COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -march=native")
  endif ()
endif ()


set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
add_subdirectory(src)

# Tools and utilities
add_subdirectory(tools)