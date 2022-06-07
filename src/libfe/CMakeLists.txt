cmake_minimum_required(VERSION 3.9.0)
project(libfe VERSION 0.1.0)

if(BUILD_SHARED_LIBS)
        add_library(libfe SHARED)
else()
        add_library(libfe STATIC)
endif()

set_target_properties(libfe PROPERTIES PREFIX "")

target_link_libraries(libfe PRIVATE libmath)

target_include_directories(libfe PRIVATE ../../include) # For blasnames, toastdef
target_include_directories(libfe INTERFACE .)           # Provide headers for dependents

# Standard build
#
target_sources(libfe PRIVATE    point.cc node.cc ndlist.cc element.cc ellist.cc surface.cc
                                ptsource.cc timespec.cc mesh.cc param.cc qmmesh.cc 
                                tri_qr.cc tri3old.cc tri3.cc tri6.cc tri6_ip.cc tri10.cc tri10_ip.cc 
                                pix4.cc vox8.cc
                                tet_qr.cc tet4.cc tet10.cc tet10_ip.cc wdg6.cc 
                                lin_qr.cc tri3D3.cc tri3D6.cc line2d2.cc
                                nim.cc 
                                refine_mesh.cc)


install(TARGETS libfe
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
