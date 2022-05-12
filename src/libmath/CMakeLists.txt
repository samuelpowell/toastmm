cmake_minimum_required(VERSION 3.0.0)
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

# Unsupported optional 
#
# ENABLE_PTHREAD: task.cc vector_MT.cc

# Assorted header defines to resolve
#
# TOAST_THREAD
# TOAST_PARALLEL
# COMPUTE_FLOPS

install(TARGETS libmath
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

