#include_directories(/usr/local/inlcude/clapack.h)

#add_executable(test_clapack test_clapack.c)
#target_link_libraries(test_clapack ${CLAPACK_LIBS} m)

#add_executable(test_clapack_eig test_clapack_eig.c)
#target_link_libraries(test_clapack_eig ${CLAPACK_LIBS} m

add_executable(test_atlas test_atlas.c)
target_link_libraries(test_atlas ${ATLAS_LIBS} m)
