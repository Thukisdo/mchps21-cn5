

add_library(utils utils.c ${INCLUDE_DIR}/utils.h)
target_include_directories(utils PUBLIC ${INCLUDE_DIR})

add_subdirectory(library_test)
add_subdirectory(poisson1D)