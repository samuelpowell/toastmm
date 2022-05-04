cmake_minimum_required(VERSION 3.0.0)
project(toastmm VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 98)
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/dist)

# liblbfgs
set(LBFGS_USE_SSE ON)
set(LBFGS_lib_TARGET_NAME liblbfgs)
add_subdirectory(extern/liblbfgs)


# Toast architecture specific & rpath settings
#

# Windows
#
if(WIN32)
    if(MSVC)
        # Make builds quieter for development
        add_compile_options(/wd4244 /wd4996 /wd4267 /wd4251 /wd4910 /wd4661)
    endif()
endif()

# MacOS
#
if(APPLE OR UNIX)
    # Make builds quieter for development
    add_compile_options(-Wno-deprecated-register -Wno-format-truncation)
    
    # Fix rpaths on install to point to the same directory
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    if(APPLE)
        set(CMAKE_INSTALL_RPATH "@loader_path/../lib")
    elseif(UNIX)
        set(CMAKE_INSTALL_RPATH "$ORIGIN:$ORIGIN/../lib")
    endif()
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    
endif()


# Linux
#
if(UNIX)

endif()

# Libraries


add_subdirectory(src/libmath)
add_subdirectory(src/libfe)
add_subdirectory(src/libstoast)
add_subdirectory(src/libfdot)

# Executables
add_subdirectory(src/supertoast)
add_subdirectory(src/bintools)

# MATLAB (mex) interface
add_subdirectory(src/matlab2)



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
