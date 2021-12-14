cmake_minimum_required(VERSION 3.0)

project(CN-TP5)

# Add the path to lapacke.h
if (NOT DEFINED LAPACK_DIR)
  message(FATAL_ERROR "LAPACK_DIR must be set prior to launching CMake")
endif ()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Am i making it too easy for the teacher ?
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

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