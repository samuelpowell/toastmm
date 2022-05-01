# Notes

Capture intent and changes applied from the initial import of the toast++.

## Build system, dependencies, and you

Aim: working one-liner build, without external dependencies, on all platforms

 - Create CMake build system for all modules required for Matlab/Python interfaces
 - Replace any numerics which don't easily conform with Eigen
 - Disable components which cannot be build (e.g. requiring mesa), note build flags

## Optional numerics / support libraries

 - reintroduce external solvers and BLAS
 - use system provided builds where possible, else e.g., `vcpkg` to avoid maintenance

## Remove unused code

Aim: simplify codebase, ensuring that all code is builds and is used


# Changelog

 - Bump liblbfgs to master for CMake support