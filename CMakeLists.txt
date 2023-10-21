# sba CMake file; see http://www.cmake.org

PROJECT(SBA)
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

# f2c is sometimes equivalent to libF77 & libI77; in that case, set HAVE_F2C to 0
SET(HAVE_F2C 0 CACHE BOOL "Do we have f2c or F77/I77?" )

# the directory where the lapack/blas/f2c libraries reside
SET(LAPACKBLAS_DIR /usr/lib CACHE PATH "Path to lapack/blas libraries")

# actual names for the lapack/blas/f2c libraries
SET(LAPACK_LIB lapack CACHE STRING "The name of the lapack library")
SET(BLAS_LIB blas CACHE STRING "The name of the blas library")
IF(HAVE_F2C)
  SET(F2C_LIB f2c CACHE STRING "The name of the f2c library")
ENDIF(HAVE_F2C)

########################## NO CHANGES BEYOND THIS POINT ##########################

# sba library source files
ADD_LIBRARY(sba STATIC
  sba_levmar.c sba_levmar_wrap.c sba_lapack.c sba_crsm.c sba_chkjac.c
)

OPTION(BUILD_DEMO "Build demo program?" TRUE)

# demo program
IF(BUILD_DEMO)
  ADD_SUBDIRECTORY(demo)
ENDIF(BUILD_DEMO)
