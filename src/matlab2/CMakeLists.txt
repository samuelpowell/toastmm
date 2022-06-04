cmake_minimum_required(VERSION 3.9.0)
project(toastmex VERSION 0.1.0)

find_package(Matlab 7.16)

if(Matlab_FOUND)

  message(STATUS "MATLAB located, toast MATLAB interface will be built")

  matlab_add_mex(NAME toastmex    SHARED R2018a 
                                  SRC toastmex.cc 
                                  matlabtoast.cc mexutil.cc util.cc 
                                  mtMesh.cc mtElement.cc mtQMvec.cc mtBasis.cc mtField.cc mtGradient.cc
                                  mtJacobian.cc mtRegul.cc toastlbfgs.cc mtSysmat_basis.cc
                                  matlabfdot.cc mtProjector.cc
                                  ../common/calc_jacobian.cc
                                  ../common/calc_gradient.cc
                                  ../common/calc_mesh.cc )                               

  target_include_directories(toastmex PRIVATE Matlab_INCLUDE_DIRS)                                        
  target_include_directories(toastmex PRIVATE ../../include)    # For blasnames, toastdef, arch

  # Add link libraries (note that plain version must be used as it has been in the FindMatlab)
  target_link_libraries(toastmex libmath)
  target_link_libraries(toastmex libfe)
  target_link_libraries(toastmex libraster)
  target_link_libraries(toastmex libstoast)
  target_link_libraries(toastmex liblbfgs)
  target_link_libraries(toastmex libfdot)
  target_link_libraries(toastmex ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY})

  # Install the output dynamic library
  install(TARGETS toastmex
          LIBRARY DESTINATION matlab/toast
          RUNTIME DESTINATION matlab/toast)

  # Install associated driver scripts, documentation, etc.
  install(DIRECTORY ${CMAKE_SOURCE_DIR}/script/matlab DESTINATION ${CMAKE_INSTALL_PREFIX})                

  # Install test meshes used by examples
  install(DIRECTORY ${CMAKE_SOURCE_DIR}/test/2D/meshes/ DESTINATION ${CMAKE_INSTALL_PREFIX}/matlab/meshes/2D)     
  install(DIRECTORY ${CMAKE_SOURCE_DIR}/test/3D/meshes/ DESTINATION ${CMAKE_INSTALL_PREFIX}/matlab/meshes/3D)     

  set(CPACK_PROJECT_NAME ${PROJECT_NAME})
  set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
  include(CPack)

else()

    message(STATUS "MATLAB not found, toast MATLAB interface will not be built")

endif()
