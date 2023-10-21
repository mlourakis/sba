CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(POSEST)

IF(MSVC)
  ADD_DEFINITIONS(/arch:SSE2)
  # get rid of CRT warnings
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(MSVC)

ADD_LIBRARY(sba sba_levmar.c sba_levmar_wrap.c sba_lapack.c sba_crsm.c sba_chkjac.c
)

OPTION(BUILD_DEMO "Build demo programs?" TRUE)

# demo program
IF(BUILD_DEMO)
  # paths to lapack/blas libraries
  SET(LAPACKBLAS_DIR "/usr/lib" CACHE PATH "Path to lapack/blas libraries")
  MESSAGE(STATUS "LAPACKBLAS_DIR = ${LAPACKBLAS_DIR}")


  # actual names for the lapack/blas/f2c libraries
  SET(LAPACKBLAS_LIB_NAMES "lapack;blas" CACHE STRING "The name of the lapack & blas libraries")
  #SET(F2C_LIB_NAME f2c CACHE STRING "The name of the f2c or F77/I77 library")
  # f2c is sometimes equivalent to libF77 & libI77
  #SET(F2C_LIB_NAME "libF77;libI77" CACHE STRING "The name of the f2c or F77/I77 library")

  SET(LIBS sba ${LAPACKBLAS_LIB_NAMES} ${F2C_LIB_NAME} m)

  LINK_DIRECTORIES(${LAPACKBLAS_DIR})

  INCLUDE_DIRECTORIES(.)
  ADD_EXECUTABLE(eucsbademo demo/eucsbademo.c demo/imgproj.c demo/readparams.c)
  TARGET_LINK_LIBRARIES(eucsbademo ${LIBS})

  MESSAGE(STATUS "demo will be linked against ${LIBS}")

# make sure that the libraries are built before the demo
  ADD_DEPENDENCIES(eucsbademo sba)
ENDIF(BUILD_DEMO)
