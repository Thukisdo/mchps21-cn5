cmake_minimum_required(VERSION 3.0)

project(CN-TP5)

set(LAPACK_DIR "/usr/include/openblas/")

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
add_subdirectory(src)