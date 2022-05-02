cmake_minimum_required(VERSION 3.0.0)
project(libstoast VERSION 0.1.0)

add_library(libstoast SHARED)

# Manged by headers
# target_compile_definitions(libstoast PRIVATE NEED_EXPLICIT_INSTANTIATION)

target_include_directories(libstoast PRIVATE ../../)                               # For toastdef
target_include_directories(libstoast PRIVATE ../../include)                        # For blasnames for toastdef
target_include_directories(libstoast PRIVATE ../../extern/eigen-3.4.0)             # Eigen for solvers, etc.

target_include_directories(libstoast INTERFACE .)     # Provide headers for dependents

# Standard build
#
# Removed: fwdsolver_zslu.cc, fwdsolver_cslu.cc (reimplemented with Eigen)
#
target_sources(libstoast PRIVATE    solution.cc mwsolution.cc source.cc fwdsolver.cc 
                                    fwdsolver_mw.cc jacobian.cc of.cc regul.cc
                                    pscaler.cc pparse.cc lsearch.cc
                                    raster.cc raster_px.cc raster_cp.cc raster_bl.cc raster_gb.cc
                                    raster_bb.cc raster_hb.cc raster_rb.cc raster_sb.cc
                                    raster2.cc raster_px2.cc raster_px2_tri.cc raster_px2_tet.cc
                                    raster_cpx.cc raster_cpx_tri.cc raster_cpx_tet.cc
                                    raster_blob2.cc raster_blob2_tri.cc raster_blob2_tet.cc raster_rb2.cc
                                    raster_bb2.cc raster_sb2.cc raster_hb2.cc raster_gb2.cc
                                    raster_cpx_tree.cc raster_cpx_tree_tri.cc
                                    sh.cc tetsplit.cc projector.cc camera.cc GLProjector.cc)

                                    
# Unsupported optional 
#
# ENABLE_MPI: jacobian_mpi.cc
#                          

target_link_libraries(libstoast libfe libmath)

if(OpenMP_CXX_FOUND)
    target_link_libraries(libstoast OpenMP::OpenMP_CXX)
endif()

install(TARGETS libstoast
        CONFIGURATIONS Release
        LIBRARY DESTINATION release/lib
        RUNTIME DESTINATION release/bin)

install(TARGETS libstoast
        CONFIGURATIONS Debug
        LIBRARY DESTINATION debug/lib
        RUNTIME DESTINATION debug/bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

