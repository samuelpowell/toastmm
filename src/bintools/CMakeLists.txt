cmake_minimum_required(VERSION 3.0.0)
project(bintools VERSION 0.1.0)

add_executable(cmpvec cmpvec.cc)
target_link_libraries(cmpvec PRIVATE libmath)
target_include_directories(cmpvec PRIVATE ../../ ../../include) 

add_executable(quad2linmesh quad2linmesh.cc)
target_link_libraries(quad2linmesh PRIVATE libmath libfe)
target_include_directories(quad2linmesh PRIVATE ../../ ../../include) 

add_executable(opt2 opt2.cc)
target_link_libraries(opt2 PRIVATE libmath libfe)
target_include_directories(opt2 PRIVATE ../../ ../../include) 

add_executable(mkmesh_sphere mkmesh_sphere.cc)
target_link_libraries(mkmesh_sphere PRIVATE libmath libfe)
target_include_directories(mkmesh_sphere PRIVATE ../../ ../../include) 

add_executable(mkmesh_cyl mkmesh_cyl.cc)
target_link_libraries(mkmesh_cyl PRIVATE libmath libfe)
target_include_directories(mkmesh_cyl PRIVATE ../../ ../../include)

add_executable(mkmesh_brick mkmesh_brick.cc)
target_link_libraries(mkmesh_brick PRIVATE libmath libfe)
target_include_directories(mkmesh_brick PRIVATE ../../ ../../include)

add_executable(meshinfo meshinfo.cc)
target_link_libraries(meshinfo PRIVATE libmath libfe)
target_include_directories(meshinfo PRIVATE ../../ ../../include) 

add_executable(makeqm3d makeqm3d.cc)
target_link_libraries(makeqm3d PRIVATE libmath libfe)
target_include_directories(makeqm3d PRIVATE ../../ ../../include) 

add_executable(makeqm makeqm.cc)
target_link_libraries(makeqm PRIVATE libmath libfe)
target_include_directories(makeqm PRIVATE ../../ ../../include) 

add_executable(lin2quadmesh lin2quadmesh.cc)
target_link_libraries(lin2quadmesh PRIVATE libmath libfe)
target_include_directories(lin2quadmesh PRIVATE ../../ ../../include) 

add_executable(lin2cubicmesh lin2cubicmesh.cc)
target_link_libraries(lin2cubicmesh PRIVATE libmath libfe)
target_include_directories(lin2cubicmesh PRIVATE ../../ ../../include) 

add_executable(getnim getnim.cc)
target_link_libraries(getnim PRIVATE libmath)
target_include_directories(getnim PRIVATE ../../ ../../include) 

install(TARGETS cmpvec quad2linmesh opt2 mkmesh_sphere mkmesh_cyl mkmesh_brick meshinfo makeqm3d makeqm lin2quadmesh lin2cubicmesh getnim
        CONFIGURATIONS Release
        LIBRARY DESTINATION release/lib
        RUNTIME DESTINATION release/bin)

        install(TARGETS cmpvec quad2linmesh opt2 mkmesh_sphere mkmesh_cyl mkmesh_brick meshinfo makeqm3d makeqm lin2quadmesh lin2cubicmesh getnim
        CONFIGURATIONS Debug
        LIBRARY DESTINATION debug/lib
        RUNTIME DESTINATION debug/bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
