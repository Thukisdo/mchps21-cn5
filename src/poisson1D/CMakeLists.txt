

add_library(poisson1D poisson1D.c ${INCLUDE_DIR}/poisson1D.h)
target_include_directories(poisson1D PUBLIC ${INCLUDE_DIR} ${LAPACK_DIR})

add_executable(poisson1D_direct poisson1D_direct.c)
target_link_libraries(poisson1D_direct poisson1D utils LAPACK::LAPACK BLAS::BLAS m)
set_target_properties(poisson1D_direct PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(poisson1D_jacobi poisson1D_jacobi.c)
target_link_libraries(poisson1D_jacobi poisson1D utils LAPACK::LAPACK BLAS::BLAS m)
set_target_properties(poisson1D_jacobi PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(poisson1D_gauss_seidel poisson1D_gauss_seidel.c)
target_link_libraries(poisson1D_gauss_seidel poisson1D utils LAPACK::LAPACK BLAS::BLAS m)
set_target_properties(poisson1D_gauss_seidel PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)