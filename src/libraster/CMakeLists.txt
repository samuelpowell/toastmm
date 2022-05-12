cmake_minimum_required(VERSION 3.0.0)
project(libraster VERSION 0.1.0)

if(BUILD_SHARED_LIBS)
        add_library(libraster SHARED)
else()
        add_library(libraster STATIC)
endif()

set_target_properties(libraster PROPERTIES PREFIX "")

target_include_directories(libraster PRIVATE ../../include)            # For blasnames, toastdef, arch
target_include_directories(libraster INTERFACE .)                      # Provide headers for dependents

# Standard build
#
target_sources(libraster PRIVATE    solution.cc mwsolution.cc sh.cc tetsplit.cc
                                    raster.cc raster_px.cc raster_cp.cc raster_bl.cc raster_gb.cc
                                    raster_bb.cc raster_hb.cc raster_rb.cc raster_sb.cc
                                    raster2.cc raster_px2.cc raster_px2_tri.cc raster_px2_tet.cc
                                    raster_cpx.cc raster_cpx_tri.cc raster_cpx_tet.cc
                                    raster_blob2.cc raster_blob2_tri.cc raster_blob2_tet.cc raster_rb2.cc
                                    raster_bb2.cc raster_sb2.cc raster_hb2.cc raster_gb2.cc
                                    raster_cpx_tree.cc raster_cpx_tree_tri.cc)


target_link_libraries(libraster libfe libmath)

install(TARGETS libraster
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

