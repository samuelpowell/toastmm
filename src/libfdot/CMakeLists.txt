cmake_minimum_required(VERSION 3.10.0)
project(libfdot VERSION 0.1.0)

if(BUILD_SHARED_LIBS)
        add_library(libfdot SHARED)
else()
        add_library(libfdot STATIC)
endif()

set_target_properties(libfdot PROPERTIES PREFIX "")

target_link_libraries(libfdot PRIVATE libmath)
target_link_libraries(libfdot PRIVATE libfe)
target_link_libraries(libfdot PRIVATE libstoast)

target_include_directories(libfdot PRIVATE ../../include)            # For blasnames, toastdef, arch
target_include_directories(libmath INTERFACE .)                      # Provide headers for dependents

# Standard build
#
target_sources(libfdot PRIVATE FDOTFwd.cc MLEMSolver.cc	matrixFreeSolver.cc muaSolver.cc util.cc)

install(TARGETS libfdot
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
