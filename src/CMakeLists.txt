

add_library(poisson1D lib_poisson1D.c ${INCLUDE_DIR}/lib_poisson1D.h)
target_include_directories(poisson1D PUBLIC ${INCLUDE_DIR} ${BLAS_DIR} ${LAPACK_DIR})
message(STATUS "CMake: poisson1D: ${BLAS_INCLUDE_DIRS} ${LAPACK_DIR}")
target_link_libraries(poisson1D LAPACK::LAPACK BLAS::BLAS)

add_executable(test_env tp_env.c ${INCLUDE_DIR}/tp_env.h ${INCLUDE_DIR}/blaslapack_headers.h)
target_include_directories(test_env PUBLIC ${INCLUDE_DIR} ${BLAS_DIR} ${LAPACK_DIR})
target_link_libraries(test_env poisson1D LAPACK::LAPACK BLAS::BLAS)

add_executable(poisson1D_direct tp2_poisson1D_direct.c ${INCLUDE_DIR}/blaslapack_headers.h)
target_include_directories(poisson1D_direct PUBLIC ${INCLUDE_DIR} ${BLAS_DIR} ${LAPACK_DIR})
target_link_libraries(poisson1D_direct poisson1D LAPACK::LAPACK BLAS::BLAS)

add_executable(poisson1D_iterat tp2_poisson1D_iter.c ${INCLUDE_DIR}/blaslapack_headers.h)
target_include_directories(poisson1D_iterat PUBLIC ${INCLUDE_DIR} ${BLAS_DIR} ${LAPACK_DIR})
target_link_libraries(poisson1D_iterat poisson1D LAPACK::LAPACK BLAS::BLAS)

add_library(utils utils.c ${INCLUDE_DIR}/utils.h)
target_include_directories(utils PUBLIC ${INCLUDE_DIR})

add_executable(dgbmv_ex4 dgbmv.c)
target_include_directories(dgbmv_ex4 PUBLIC ${INCLUDE_DIR})
target_link_libraries(dgbmv_ex4 poisson1D LAPACK::LAPACK BLAS::BLAS utils)