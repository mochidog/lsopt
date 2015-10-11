add_library(lsopt lsopt.c)
target_link_libraries(lsopt ${CLAPACK_LIBS} m)

add_executable(benchmark benchmark.c)
target_link_libraries(benchmark lsopt)