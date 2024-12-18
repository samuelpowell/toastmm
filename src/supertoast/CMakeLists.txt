cmake_minimum_required(VERSION 3.10.0)
project(supertoast VERSION 0.1.0)

# fwdfem
add_executable(fwdfem fwdfem.cc)
target_include_directories(fwdfem PRIVATE ../../include)
target_link_libraries(fwdfem libstoast libfe libmath)

# fwdfem_cw
add_executable(fwdfem_cw fwdfem_cw.cc)
target_include_directories(fwdfem_cw PRIVATE ../../include)
target_link_libraries(fwdfem_cw libstoast libfe libmath)

# fwdfem_tpsf
add_executable(fwdfem_tpsf fwdfem_tpsf.cc)
target_include_directories(fwdfem_tpsf PRIVATE ../../include)
target_link_libraries(fwdfem_tpsf libstoast libfe libmath)
        
# supertoast
add_executable(supertoast supertoast.cc supertoast_util.cc util.cc solver.cc solverlm.cc
                          solverpcg.cc solverlin.cc
                          solverbfgs.cc solverlbfgs.cc solverart.cc solverblockart.cc)
target_include_directories(supertoast PRIVATE ../../include)
target_link_libraries(supertoast liblbfgs libstoast libfe libmath)

# supertoast_mw
add_executable(supertoast_mw supertoast_mw.cc supertoast_util.cc util.cc solver_mw.cc solverlm_mw.cc
                             solverpcg_mw.cc 
                             solverlbfgs_mw.cc solverlm2_mw.cc)
target_include_directories(supertoast_mw PRIVATE ../../include)
target_link_libraries(supertoast_mw liblbfgs libstoast libfe libmath)

# supertoast_cw_mw
add_executable(supertoast_cw_mw supertoast_cw_mw.cc supertoast_util.cc util.cc solver_cw.cc solverlm_cw_mw.cc
                                solverpcg_cw_mw.cc 
                                solverlbfgs_cw_mw.cc solverlm2_cw_mw.cc)
target_include_directories(supertoast_cw_mw PRIVATE ../../include)
target_link_libraries(supertoast_cw_mw liblbfgs libstoast libfe libmath)


install(TARGETS fwdfem fwdfem_cw fwdfem_tpsf supertoast supertoast_mw supertoast_cw_mw 
        RUNTIME DESTINATION bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

