

function(try_find_blas)

  if (${ARGC} EQUAL 0)
    message(STATUS No BLAS vendor specified, looking for MKL)
    set(BLA_VENDOR Intel10_64lp_seq)
  else()
    message(STATUS "Looking for BLAS vendor ${ARGV0}")
    set(BLA_VENDOR ${ARGV0})
  endif()


  find_package(BLAS QUIET)
  find_package(LAPACK QUIET)

  if (NOT LAPACK_FOUND OR NOT BLAS_FOUND)
    message(STATUS "${ARGV0} vendor not found, looking for another BLAS vendor")
    set(BLA_VENDOR "")
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_CBLAS" PARENT_SCOPE)
  else()
    message(STATUS "${ARGV0} vendor found")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_CBLAS" PARENT_SCOPE)
  endif ()
endfunction()