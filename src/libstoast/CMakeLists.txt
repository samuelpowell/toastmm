cmake_minimum_required(VERSION 3.9.0)
project(libstoast VERSION 0.1.0)

if(BUILD_SHARED_LIBS)
        add_library(libstoast SHARED)
else()
        add_library(libstoast STATIC)
endif()

set_target_properties(libstoast PROPERTIES PREFIX "")

target_include_directories(libstoast PRIVATE ../../include)            # For blasnames, toastdef, arch
target_include_directories(libstoast PUBLIC ../../extern/eigen-3.4.0)  # Eigen for solvers, etc.
target_include_directories(libstoast INTERFACE .)                      # Provide headers for dependents

# Standard build
#
# Removed: fwdsolver_zslu.cc, fwdsolver_cslu.cc (reimplemented with Eigen)
#
target_sources(libstoast PRIVATE    source.cc fwdsolver.cc 
                                    fwdsolver_mw.cc jacobian.cc of.cc regul.cc
                                    pscaler.cc pparse.cc lsearch.cc
                                    projector.cc camera.cc)


target_link_libraries(libstoast libfe libmath libraster)
target_link_libraries(libstoast SuiteSparse::suitesparseconfig
                                SuiteSparse::amd
                                SuiteSparse::btf
                                SuiteSparse::camd
                                SuiteSparse::ccolamd
                                SuiteSparse::colamd
                                SuiteSparse::cholmod)

if(TOAST_THREAD)
        # Enable threading in eigen
        # target_link_libraries(libstoast OpenMP::OpenMP_CXX)
endif()

install(TARGETS libstoast
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

