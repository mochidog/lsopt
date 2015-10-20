cmake_minimum_required(VERSION 3.2)

set(CMAKE_C_FLAGS "-O3 -Wall -g")

#set(CLAPACK_DIR /usr/local/clapack-3.2.1)
#find_library(LAPACK_LIB lapack HINTS ${CLAPACK_DIR})
#find_library(BLAS_LIB blas HINTS ${CLAPACK_DIR})
##find_library(TMG_LIB tmglib HINTS ${CLAPACK_DIR})
#find_library(F2C_LIB f2c HINTS ${CLAPACK_DIR})
#
#set(CLAPACK_LIBS
#    ${LAPACK_LIB}
#    ${BLAS_LIB}
#    ${F2C_LIB}
#)

#message(${CLAPACK_LIBS})
#
#include_directories(${CLAPACK_DIR}/include)

set(ATLAS_DIR /usr/local/atlas)
find_library(LAPACK_LIB lapack HINTS ${ATLAS_DIR}/lib)
find_library(CBLAS_LIB cblas HINTS ${ATLAS_DIR}/lib)
find_library(ATLAS_LIB atlas HINTS ${ATLAS_DIR}/lib)
set(ATLAS_LIBS
    ${LAPACK_LIB}
    ${CBLAS_LIB}
    ${ATLAS_LIB}
)

message(${ATLAS_LIBS})

include_directories(${ATLAS_DIR}/include)

#add_subdirectory(src)
add_subdirectory(sandbox)
