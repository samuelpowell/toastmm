# Notes

Capture intent and changes applied from the initial import of the toast++.

## Build system, dependencies, and you

Aim: working one-liner build, without external dependencies, on all platforms

 - Create CMake build system for all modules required for Matlab/Python interfaces
 - Replace any numerics which don't easily conform with Eigen
 - Disable components which cannot be build (e.g. requiring mesa), note build flags

## Optional numerics / support libraries

 - reintroduce external solvers and BLAS as options
 - use system provided builds where possible, else e.g., `vcpkg` to avoid maintenance

## Threading

 - manage threading within toast, and optional BLAS

## Remove features and unused code

Aim: simplify codebase, ensuring that all code is builds and is used

 - CUDA, MPI

# Changelog (build system)

 - Bump liblbfgs to master for CMake support
 - Add Eigen
 - Remove LU decomposition from TCompRowMatrix, implementation (`crmatrix.cc`) removed from build and decleration commented, remove SuperLU headers from `crmatrix_cm.cc`.
 - Implement libmath CMakeLists.txt for default build without options
 - Add `#define FELIB_IMPLEMENTATION` to `wdg18inf.cc` for proper symbol import/export on Windows
 - Implement libfe CMakeLists.txt for default build without options
 - Remove complex<T> SuperLU solver engines from build (`fwdsolver_*slu.cc`), and exclude headers
 - Replace functionalisty of SuperLU solver engines with Eigen templated LU
 - Fix headers that do not actually require SuperLU
 - Implement CMakeLists.txt for libstoast
 - Implement CMakeLists.txt for libfdot
 - Implement CMakeLists.txt for matlab2 MEX build
 - Provide interface headers for subprojects, pick up headers automatically
 - Fix gitignore exclusion of liblbfgs and Eigen and reimport
 - Fixup nr.cc for MacOS Clang by removing definition of `vector` and using `dvector` instead
 - Need explicit instantiation required on MacOS
 - Standard CXX 98 required for MacOS (owing to narrowing in matlab interface) which requires spacing on double template args
 - Implement CMakeLists.txt for all bintools targets
 - Enable explicit instantiation in fwdsolver_mw.h  (linux build failed with this, perhaps this could be a Clang thing, added to TODO)
 - Add STOASTLIB (dll import/export) to MWsolution class
 - Implement CMakeLists.txt for all buildable supertoast targets (NB: some Makefile targets no longer exist or reference old code)
 - Purge unused numerics, configuration, projects
 - Purge MPI support
 - Purge CUDA support
 - Cleanup core headers
   - Make all defines check for `_WIN32` not `WIN32`, remove setting from CMake
   - Remove unused `TOASTLIB` define
   - Remove expiry functions (not used) and defines
   - Remove USE_SPLBLAS (unued)



# TODO

 - Check Makefile / Xcode projects configuration for shared library flags, e.g. -fPIC required on Linux. Consider `FLAGS`, `SHLIB_CFLAGS`
 - Look at building static libraries for intermediary components
 - Investigate `warning C4910` on MSVC - some declspec conflict
 - Look at MEX config as per first item
 - Check `MESA_SUPPORT`
 - Make individual libraries properly CMake with interface exports, etc. to avoid replicating includes in e.g. matlab2#
 - Look at fwdsolver_mw.h instantiation requirements, determine appropriate preprocessor gaurd (e.g. Clang?)
 - Remove MESA based projection
 - Resolve `TOAST_THREAD`, `TOAST_THREAD_LEVEL`, `TOAST_THREAD_MATLAB_GRADIENT`, `TOAST_THREAD_MATLAB_QMVEC`, `TOAST_THREAD_ASSEMBLE`
- 

