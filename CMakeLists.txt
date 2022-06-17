cmake_minimum_required(VERSION 3.0.0)

set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake"
  CACHE STRING "Vcpkg toolchain file")

if(WIN32)
    set(VCPKG_TARGET_TRIPLET "x64-windows"
        CACHE STRING "Vcpkg triplet")
endif()

project(toastmm VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/dist)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(BUILD_SHARED_LIBS OFF)
set(TOAST_THREAD ON)

if(NOT BUILD_SHARED_LIBS)
    add_compile_definitions(TOAST_STATIC)
    # set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

if(TOAST_THREAD)
    add_compile_definitions(TOAST_THREAD)
    find_package(Threads REQUIRED)
    # find_package(OpenMP REQUIRED)
endif()

find_package(suitesparse CONFIG REQUIRED)

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
        if(BUILD_SHARED_LIBS)
            # Shared library builds are ecxeptionally noise due to various unresolved import/export
            # issues, we quieten this for now to concentrate on the bad stuff
            add_compile_options(/wd4244 /wd4996 /wd4267 /wd4251 /wd4910 /wd4661 /wd4305)
        else()
            # Static build is okay but unsafe warnings should be addressed
            add_compile_options(/wd4244 /wd4996)
        endif()
    endif()
endif()

# MacOS and Linux
#
if(APPLE OR UNIX)
    # Make builds quieter for development
    add_compile_options(-Wno-comment -Wno-unused-variable -Wno-unused-result)

    if(NOT APPLE)
        add_compile_options(-Wno-format-truncation)
    endif()
    
    # Fix rpaths on install to point to the lib subdirectory and the same directory
    # to support both normal installations (e.g. MATLAB) and Python packages
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    if(APPLE)
        set(CMAKE_INSTALL_RPATH "@loader_path/../lib")
    elseif(UNIX)
        set(CMAKE_INSTALL_RPATH "$ORIGIN:$ORIGIN/../lib")
    endif()
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    
endif()

# Libraries
add_subdirectory(src/libmath)
add_subdirectory(src/libfe)
add_subdirectory(src/libraster)
add_subdirectory(src/libstoast)
add_subdirectory(src/libfdot)

# Executables
add_subdirectory(src/supertoast)
add_subdirectory(src/bintools)

# Interfaces
#
# We only build the MATLAB and Python interfaces under a static configuration, as this 
# removes significant overhead in managing e.g., rpaths, and bundling of dynamic libraries.
#
if(BUILD_SHARED_LIBS)
  message(STATUS "Shared library build select, MATLAB & Python interfaces will not be built")
else()

    # MATLAB (mex) interface
    add_subdirectory(src/matlab2)

    # Python interface 
    #
    # Note that the python interface is simply an install target, which writes the source
    # files into the distribution directory, compilation is performed using python tools
    # in order to enable binary wheels to be built using, e.g., cibuildwheel. For a local
    # build one may execute, e.g., `python setup.py bdist_wheel` in the output directory.
    add_subdirectory(src/python)

endif()


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
