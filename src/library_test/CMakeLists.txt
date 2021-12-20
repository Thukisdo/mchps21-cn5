

add_executable(exo_dgbmv exo_dgbmv.c)
target_include_directories(exo_dgbmv PRIVATE ${INCLUDE_DIR})
target_link_libraries(exo_dgbmv utils poisson1D)

add_executable(exo_ludecomposition exo_ludecomposition.c)
target_link_libraries(exo_ludecomposition utils poisson1D)
target_include_directories(exo_ludecomposition PRIVATE ${INCLUDE_DIR})