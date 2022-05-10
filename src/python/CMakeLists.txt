cmake_minimum_required(VERSION 3.0.0)
project(pytoast VERSION 0.1.0)

# Configure setup.py
#
# We encode the full path to headers in the setup file in order to ensure that the 
# build is relocatable.
#
# Note that we include the libmath includes in addition to the interface includes
# in order to pull in the full root path (required for toastdef.h arch.h)
get_target_property(ROOT_INCLUDES       libmath     INCLUDE_DIRECTORIES)
get_target_property(LIBMATH_INCLUDES    libmath     INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(LIBFE_INCLUDES      libfe       INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(LIBSTOAST_INCLUDES  libstoast   INTERFACE_INCLUDE_DIRECTORIES)

foreach(dir ${ROOT_INCLUDES} ${LIBMATH_INCLUDES} ${LIBFE_INCLUDES} ${LIBSTOAST_INCLUDES})
    string(APPEND PYTOAST_INCLUDE_STRING "'${dir}',")
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

# Install libraries 
#
# We place the built toast libraries in the module in order that seuptools will install
# these are project files. 
#
if(WIN32)
    list(APPEND PYTOAST_LIBS ${CMAKE_INSTALL_PREFIX}/bin/libmath.dll
                             ${CMAKE_INSTALL_PREFIX}/bin/libfe.dll
                             ${CMAKE_INSTALL_PREFIX}/bin/libstoast.dll)

    foreach ( file ${PYTOAST_LIBS} )
        install( FILES ${file} DESTINATION ${CMAKE_INSTALL_PREFIX}/python/src/toast )
    endforeach()

    list(APPEND PYTOAST_IMPORT_LIBS ${CMAKE_INSTALL_PREFIX}/lib/libmath.lib
                                    ${CMAKE_INSTALL_PREFIX}/lib/libfe.lib
                                    ${CMAKE_INSTALL_PREFIX}/lib/libstoast.lib)

    foreach ( file ${PYTOAST_IMPORT_LIBS} )
        install( FILES ${file} DESTINATION ${CMAKE_INSTALL_PREFIX}/python/src/toast/lib )
    endforeach()

endif()

if(APPLE OR UNIX)
    if(APPLE)
        list(APPEND PYTOAST_LIBS ${CMAKE_INSTALL_PREFIX}/lib/libmath.dylib
                                 ${CMAKE_INSTALL_PREFIX}/lib/libfe.dylib
                                 ${CMAKE_INSTALL_PREFIX}/lib/libstoast.dylib)
    else()
        list(APPEND PYTOAST_LIBS ${CMAKE_INSTALL_PREFIX}/lib/libmath.so
                                 ${CMAKE_INSTALL_PREFIX}/lib/libfe.so
                                 ${CMAKE_INSTALL_PREFIX}/lib/libstoast.so)
    endif()

    foreach ( file ${PYTOAST_LIBS} )
        install( FILES ${file} DESTINATION ${CMAKE_INSTALL_PREFIX}/python/src/toast/lib )
    endforeach()
endif()


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

