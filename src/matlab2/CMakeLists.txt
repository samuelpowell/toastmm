cmake_minimum_required(VERSION 3.0.0)
project(toastmex VERSION 0.1.0)

find_package(Matlab 7.16)

if(Matlab_FOUND)

        message(STATUS "MATLAB located, toast MATLAB interface will be built")

        matlab_add_mex(NAME toastmex    SHARED R2017b
                                        SRC toastmex.cc 
                                        matlabtoast.cc mexutil.cc util.cc 
                                        mtMesh.cc mtElement.cc mtQMvec.cc mtBasis.cc mtField.cc mtGradient.cc
                                        mtJacobian.cc mtRegul.cc toastlbfgs.cc mtSysmat_basis.cc
                                        matlabfdot.cc mtProjector.cc )

        target_include_directories(toastmex PRIVATE Matlab_INCLUDE_DIRS)                                        
        target_include_directories(toastmex PRIVATE ../../)                            # For toastdef
        target_include_directories(toastmex PRIVATE ../../include)                     # For blasnames for toastdef
        target_include_directories(toastmex PRIVATE ../libmath)                        # For mathlib.h
        target_include_directories(toastmex PRIVATE ../libstoast)                      # For stoastlib.h
        target_include_directories(toastmex PRIVATE ../libfe)                          # For felib.h
        target_include_directories(toastmex PRIVATE ../libfdot)                        # For FDOTFwd.h
        target_include_directories(toastmex PRIVATE ../../extern/liblbfgs/include)   # For lbfgs
        target_include_directories(toastmex PRIVATE ../../extern/eigen-3.4.0)

        # Add link libraries (note that plain version must be used as it has been in the FindMatlab)
        target_link_libraries(toastmex libstoast)
        target_link_libraries(toastmex liblbfgs)
        target_link_libraries(toastmex libfdot)
        target_link_libraries(toastmex ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY})

        install(TARGETS toastmex
                LIBRARY DESTINATION lib
                RUNTIME DESTINATION bin)

        set(CPACK_PROJECT_NAME ${PROJECT_NAME})
        set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
        include(CPack)

else()

        message(STATUS "MATLAB not found, toast MATLAB interface will not be built")

endif()

