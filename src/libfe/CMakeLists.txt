cmake_minimum_required(VERSION 3.0.0)
project(libfe VERSION 0.1.0)

add_library(libfe SHARED)

target_link_libraries(libfe PRIVATE libmath)

# Managed by headers
# target_compile_definitions(libfe PRIVATE NEED_EXPLICIT_INSTANTIATION)

target_include_directories(libfe PRIVATE ../libmath)
target_include_directories(libfe PRIVATE ../../)                # For toastdef
target_include_directories(libfe PRIVATE ../../include)         # For blasnames for toastdef

# Standard build
#
target_sources(libfe PRIVATE    point.cc node.cc ndlist.cc element.cc ellist.cc surface.cc
                                ptsource.cc timespec.cc mesh.cc param.cc qmmesh.cc 
                                tri3old.cc tri3.cc tri6.cc tri6_ip.cc tri10.cc tri10_ip.cc 
                                quad4.cc pix4.cc vox8.cc vox27.cc 
                                tet4.cc tet10.cc tet10_ip.cc wdg6.cc wdg18inf.cc 
                                tri3D3.cc tri3D6.cc line2d2.cc
                                lin_qr.cc tri_qr.cc tet_qr.cc 
                                ordmmd.cc rcm.cc meshopt.cc nim.cc 
                                refine_mesh.cc)

# Unsupported optional 
#
# ENABLE_MPI: mesh_mpi.cc
#

# Assorted header defines to resolve
#
# TOAST_MPI

install(TARGETS libfe
        CONFIGURATIONS Release
        LIBRARY DESTINATION release/lib
        RUNTIME DESTINATION release/bin)

install(TARGETS libfe
        CONFIGURATIONS Debug
        LIBRARY DESTINATION debug/lib
        RUNTIME DESTINATION debug/bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
