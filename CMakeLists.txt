cmake_minimum_required(VERSION 3.0.0)
project(toastmm VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/dist)

# liblbfgs
set(LBFGS_USE_SSE ON)
set(LBFGS_lib_TARGET_NAME liblbfgs)
add_subdirectory(extern/liblbfgs)

# Libraries

# Make windows builds slightly quieter for development
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4244 /wd4996 /wd4267 /wd4251")

add_subdirectory(src/libmath)
# add_subdirectory(src/libfe)
# add_subdirectory(src/libstoast)

# Executables
# add_subdirectory(src/supertoast)
# add_subdirectory(src/bintools)

# MATLAB (mex) interface
# add_subdirectory(src/matlab2)



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
