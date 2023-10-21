# CMake file for sba's demo program

INCLUDE_DIRECTORIES(..)
LINK_DIRECTORIES(.. ${LAPACKBLAS_DIR})

ADD_EXECUTABLE(eucsbademo eucsbademo.c imgproj.c readparams.c eucsbademo.h readparams.h)
# libraries the demo depends on
IF(HAVE_F2C)
  TARGET_LINK_LIBRARIES(eucsbademo sba ${LAPACK_LIB} ${BLAS_LIB} ${F2C_LIB} m)
ELSE(HAVE_F2C)
  TARGET_LINK_LIBRARIES(eucsbademo sba ${LAPACK_LIB} ${BLAS_LIB} m)
ENDIF(HAVE_F2C)

# make sure that the library is built before the demo
ADD_DEPENDENCIES(eucsbademo sba)
