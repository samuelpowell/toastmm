cmake_minimum_required(VERSION 3.9.0)
project(libmath VERSION 0.1.0)

if(BUILD_SHARED_LIBS)
        add_library(libmath SHARED)
else()
        add_library(libmath STATIC)
endif()

set_target_properties(libmath PROPERTIES PREFIX "")

target_include_directories(libmath PRIVATE ../../include)     # For blasnames, toastdef, arch
target_include_directories(libmath INTERFACE .)               # Provide headers for dependents

# Standard build
#
# Removed: crmatrix.cc
#
target_sources(libmath PRIVATE  mathdef.cc util.cc timing.cc vector.cc gsmatrix.cc 
                                scrmatrix.cc crmatrix_cm.cc precon.cc nr.cc error.cc cr_cholesky.cc 
                                fourn.cc ilutoast.cc util.cc timing.cc mathdef.cc)

# Otional 
#
if(TOAST_THREAD)
        target_sources(libmath PRIVATE task.cc)
        target_link_libraries(libmath PUBLIC Threads::Threads)
        # Enable threading in Eigen
        target_link_libraries(libmath PUBLIC OpenMP::OpenMP_CXX)
endif()

# Assorted header defines to resolve
#
# COMPUTE_FLOPS

install(TARGETS libmath
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

