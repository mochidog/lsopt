cmake_minimum_required(VERSION 3.2)

set(CMAKE_C_FLAGS "-O3 -Wall -g")

set(CLAPACK_DIR /usr/local/clapack-3.2.1)
find_library(LAPACK_LIB lapack HINTS ${CLAPACK_DIR})
find_library(BLAS_LIB blas HINTS ${CLAPACK_DIR})
#find_library(TMG_LIB tmglib HINTS ${CLAPACK_DIR})
find_library(F2C_LIB f2c HINTS ${CLAPACK_DIR})

set(CLAPACK_LIBS
    ${LAPACK_LIB}
    ${BLAS_LIB}
    ${F2C_LIB}
)

message(${CLAPACK_LIBS})

include_directories(${CLAPACK_DIR}/include)

add_subdirectory(src)
#add_subdirectory(sandbox)
