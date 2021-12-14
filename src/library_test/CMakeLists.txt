

add_executable(exo_dgbmv exo_dgbmv.c)
target_link_libraries(exo_dgbmv utils m BLAS::BLAS LAPACK::LAPACK poisson1D)
target_include_directories(exo_dgbmv PRIVATE ${INCLUDE_DIR} ${LAPACK_DIR})