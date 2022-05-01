cmake_minimum_required(VERSION 3.0.0)
project(libfdot VERSION 0.1.0)

add_library(libfdot SHARED)

target_link_libraries(libfdot PRIVATE libmath)
target_link_libraries(libfdot PRIVATE libfe)
target_link_libraries(libfdot PRIVATE libstoast)

# Managed by headers
# target_compile_definitions(libmath PRIVATE NEED_EXPLICIT_INSTANTIATION)

target_include_directories(libfdot PRIVATE ../../)                      # For toastdef
target_include_directories(libfdot PRIVATE ../../include)               # For blasnames for toastdef
target_include_directories(libfdot PRIVATE ../libstoast)                # For stoastlib.h
target_include_directories(libfdot PRIVATE ../libfe)                    # For felib.h
target_include_directories(libfdot PRIVATE ../libmath)                  # For mathlib.h
target_include_directories(libfdot PRIVATE ../../extern/eigen-3.4.0)    # Eigen for solvers, etc.

# Standard build
#
target_sources(libfdot PRIVATE FDOTFwd.cc MLEMSolver.cc	matrixFreeSolver.cc muaSolver.cc util.cc)

install(TARGETS libfdot
        CONFIGURATIONS Release RelWithDebInfo
        LIBRARY DESTINATION release/lib
        RUNTIME DESTINATION release/bin)

install(TARGETS libfdot
        CONFIGURATIONS Debug
        LIBRARY DESTINATION debug/lib
        RUNTIME DESTINATION debug/bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
