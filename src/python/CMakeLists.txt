cmake_minimum_required(VERSION 3.9.0)
project(pytoast VERSION 0.1.0)


# Configure setup.py
#
# We encode the full path to headers and the static libraries in the setup file in order that
# the build can be performed in an isolated environment (on the same host).

# Library headers
#
get_target_property(LIBMATH_INCLUDES    libmath     INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(LIBFE_INCLUDES      libfe       INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(LIBRASTER_INCLUDES  libraster   INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(LIBSTOAST_INCLUDES  libstoast   INTERFACE_INCLUDE_DIRECTORIES)

foreach(dir "${CMAKE_SOURCE_DIR}/include/."
            "${CMAKE_SOURCE_DIR}/src/common/."
            ${LIBMATH_INCLUDES} 
            ${LIBFE_INCLUDES} 
            ${LIBRASTER_INCLUDES}
            ${LIBSTOAST_INCLUDES})

    file(RELATIVE_PATH reldir ${CMAKE_SOURCE_DIR} ${dir})
    string(APPEND PYTOAST_INCLUDE_STRING "'${reldir}',")
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
    file(RELATIVE_PATH reldir ${CMAKE_SOURCE_DIR} ${dir})
    if(WIN32)
        # Very ugly hack because I can't access ARCHIVE_CUI
        string(APPEND PYTOAST_LIBRARY_STRING "'${reldir}/Release',")
    else()
        string(APPEND PYTOAST_LIBRARY_STRING "'${reldir}',")
    endif()
endforeach()


if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(PYTOAST_DEBUG_FLAG "True")
else()
    set(PYTOAST_DEBUG_FLAG "False")
endif()

if(TOAST_THREAD)
    set(PYTOAST_THREAD_FLAG "True")
else()
    set(PYTOAST_THREAD_FLAG "False")
endif()

configure_file(setup.py.in setup.py @ONLY)

configure_file(setup.py.in ${CMAKE_SOURCE_DIR}/setup.py @ONLY COPYONLY)
configure_file(pyproject.toml ${CMAKE_SOURCE_DIR}/pyproject.toml COPYONLY)

# Note that we have to place the resultant file in the root, because the CI wheel builder
# requires this to be the case in order to import the entire project. Evidently this is
# sub-optimal since it polutes the source tree. The python project TOML configuration is 
# also required.
#
# install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.py DESTINATION ${CMAKE_INSTALL_PREFIX}/python/src)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.py DESTINATION ${CMAKE_SOURCE_DIR})
install( FILES pyproject.toml DESTINATION ${CMAKE_SOURCE_DIR})

# Install python source, associated driver scripts, documentation, etc.
#
install(DIRECTORY ${CMAKE_SOURCE_DIR}/script/python DESTINATION ${CMAKE_INSTALL_PREFIX}) 

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

