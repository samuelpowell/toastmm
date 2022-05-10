cmake_minimum_required(VERSION 3.0.0)
project(pytoast VERSION 0.1.0)


# Configure setup.py
#
# We encode the full path to headers and the static libraries in the setup file in order that
# the build can be performed in an isolated environment (on the same host).

# Library headers
#
# Note that we include the libmath includes in addition to the interface includes
# in order to pull in the full root path (required for toastdef.h arch.h)
get_target_property(ROOT_INCLUDES       libmath     INCLUDE_DIRECTORIES)
get_target_property(LIBMATH_INCLUDES    libmath     INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(LIBFE_INCLUDES      libfe       INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(LIBRASTER_INCLUDES  libraster   INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(LIBSTOAST_INCLUDES  libstoast   INTERFACE_INCLUDE_DIRECTORIES)

foreach(dir ${ROOT_INCLUDES} 
            ${LIBMATH_INCLUDES} 
            ${LIBFE_INCLUDES} ${LIBRASTER_INCLUDES}
            ${LIBSTOAST_INCLUDES})

    string(APPEND PYTOAST_INCLUDE_STRING "'${dir}',")
endforeach()

# Library link directories
#
get_target_property(LIBMATH_DIR    libmath      BINARY_DIR)
get_target_property(LIBFE_DIR      libfe        BINARY_DIR)
get_target_property(LIBRASTER_DIR  libraster    BINARY_DIR)
get_target_property(LIBSTOAST_DIR  libstoast    BINARY_DIR)

foreach(dir ${LIBMATH_DIR} 
            ${LIBFE_DIR} 
            ${LIBRASTER_DIR}
            ${LIBSTOAST_DIR})
    if(WIN32)
        # Very ugly hack because I can't access ARCHIVE_CUI
        string(APPEND PYTOAST_LIBRARY_STRING "'${dir}/Release',")
    else()
        string(APPEND PYTOAST_LIBRARY_STRING "'${dir}',")
    endif()
endforeach()


if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(PYTOAST_DEBUG_FLAG "True")
else()
    set(PYTOAST_DEBUG_FLAG "False")
endif()

configure_file(setup.py.in setup.py @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.py DESTINATION ${CMAKE_INSTALL_PREFIX}/python/src)

# Install python source, associated driver scripts, documentation, etc.
#
install(DIRECTORY toast DESTINATION ${CMAKE_INSTALL_PREFIX}/python/src)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/script/python DESTINATION ${CMAKE_INSTALL_PREFIX})       

# Copy files required for build
#
list(APPEND PYTOAST_BUILD_FILES ../../include/arch.h
                                pyproject.toml
                                toastmodule.cc)

# Install python build system files
foreach ( file ${PYTOAST_BUILD_FILES} )
    get_filename_component( dir ${file} DIRECTORY )
    install( FILES ${file} DESTINATION ${CMAKE_INSTALL_PREFIX}/python/src )
endforeach()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

