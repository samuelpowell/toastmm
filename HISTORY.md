# Toast-- Development History


## v0.9.1

 - Add build for Apple Silicon
 - Bump CI

## v0.9.0

Initial release:
 - MATLAB interface compatible with Toast++ 2.0.2
 - Python interface unstable

## Development from Toast++

 - Bump liblbfgs to master for CMake support
 - Add Eigen
 - Remove LU decomposition from TCompRowMatrix, implementation (`crmatrix.cc`) removed from build and decleration commented, remove SuperLU headers from `crmatrix_cm.cc`.
 - Implement CMakeLists.txt for libmath, libstoast, libfdot, matlab MEX
 - Add `#define FELIB_IMPLEMENTATION` to `wdg18inf.cc` for proper symbol import/export on Windows
 - Implement libfe CMakeLists.txt for default build without options
 - Remove complex<T> SuperLU solver engines from build (`fwdsolver_*slu.cc`), and exclude headers
 - Replace functionalisty of SuperLU solver engines with Eigen templated LU
 - Fix headers that do not actually require SuperLU
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
 - Python interface development
   - Complete refactor of python build with view to CI builds of wheels
   - CMake target configures a setuptools script with full paths, installs dynamic and import libraries in the tree for automatic installation
   - CMake install target provides a tree that can be built using `python -m pip wheel` or `python -m build` or `python setup.py bdist_wheel`. The former are isolated builds.
   - Update all reconstruction scripts, noting
     - (silent) errors when basis mapping input are not rank-1
     - Jacobians take handles to bases, rathter than just taking the basis
     - Arugement checking throughout module
     - Move from np.matrix to np.ndarray through module and in examples
     - Reorganise examples
 - Clean up part two
   - Refactor headers
   - Remove unused eigpair.cc, arpack.cc, crmatric.cc (LU)
   - Python interface integer refactor
 - Refactor raster out of libstoast and into its own module (libraster)
 - Provide for static build of libraries to output single mex and toast modules without dependency 'challenges'
 - Multithreading
    - Rework coarse parallelism to use C++ threads
    - Enable `TOAST_THREAD` at `THREAD_LEVEL_2` by default, providing parallel assembly, source/meas construction, some Jacobian calculations, iterative solvers over QM
    - Remove `THREAD_LEVEL_1` (which uses the thread pool, not ported), and make `THREAD_LEVEL_2 == TOAST_THREAD`, remove vector_MT.
    - Enable multi-threading in Python interface
    - Thread pool based parallelism removed as used in limited places. Will consider reimplementation with OpenMP or standard library implementation: `TOAST_PARALLEL` define, parallel CG implementation, thread pool implementation.
  - Refactor Python build configuration to be standalone, and run from the source root, in order to use cibuildwheel in CI
  - Enale CI for Python wheels and MATLAB mex builds
  - Remove remaining MESA based functionality
  - Remove mesh reordering functionality due to LGPL license
  - Remove ILUPACK bindings
  - Remove unused solver implementations
  - Fix swap bug with internal Cholesky implementation, use Eigen Cholesky for forward solvers
  - Split factorise/analyse steps in forward solver
  - toastFields reimplemented with both MATLAB and TOAST solvers, the former used by default for direct solves
  - toastJacobianCW reimplemented to use toastFields
  - toastJacobianCW mesh integrals multithreaded, 6x speedup on 8x threads
  - Update MEX interface to interleaved complex (>R2018a) format with runtime improvements
  - Remove old MATLAB interface scripts
  - Cache sparsity structure inside mesh
  - Use C++ stdlib parallel sort for sparsity structure computation
  - Compound matrix assembly works over all parameters at once, reducing time for local summation
  - Make single precision (for CUDA) and the new raster interface optional compile time defines
