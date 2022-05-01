cmake_minimum_required(VERSION 3.0.0)
project(libmath VERSION 0.1.0)

add_library(libmath SHARED)

target_compile_definitions(libmath PRIVATE NEED_EXPLICIT_INSTANTIATION)

target_include_directories(libmath PRIVATE ../../)            # For toastdef
target_include_directories(libmath PRIVATE ../../include)     # For blasnames for toastdef

# Standard build
#
# Removed: arpach.cc, crmatrix.cc
#
target_sources(libmath PRIVATE  mathdef.cc util.cc timing.cc vector.cc gsmatrix.cc 
                                scrmatrix.cc crmatrix_cm.cc precon.cc nr.cc error.cc cr_cholesky.cc 
                                fourn.cc ilutoast.cc util.cc timing.cc mathdef.cc)

# Unsupported optional 
#
# USE_BLZPACK: eigpair.cc
# ENABLE_PTHREAD: task.cc vector_MT.cc
# ENABLE_MPI crmatrix_mpi.cc dnsmatrix_mpi.cc toast_mpi.cc
# ENABLE_CUDA vector_cusp.cu crmatrix_cusp.cc

# Assorted header defines to resolve
#
# TOAST_THREAD
# TOAST_MPI
# TOAST_PARALLEL
# COMPUTE_FLOPS

install(TARGETS libmath
        CONFIGURATIONS Release RelWithDebInfo
        LIBRARY DESTINATION release/lib
        RUNTIME DESTINATION release/bin)

install(TARGETS libmath
        CONFIGURATIONS Debug
        LIBRARY DESTINATION debug/lib
        RUNTIME DESTINATION debug/bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
