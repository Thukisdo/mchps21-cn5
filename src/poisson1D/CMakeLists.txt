

add_library(poisson1D poisson1D.c ${INCLUDE_DIR}/poisson1D.h)
target_include_directories(poisson1D PUBLIC ${INCLUDE_DIR} ${LAPACK_DIR})

add_executable(poisson1D_direct poisson1D_direct.c)
target_link_libraries(poisson1D_direct poisson1D utils LAPACK::LAPACK BLAS::BLAS)

add_executable(poisson1D_iter poisson1D_iter.c)
target_link_libraries(poisson1D_iter poisson1D utils LAPACK::LAPACK BLAS::BLAS)