# Toast-- Development Notes

## TODO

 - Investigate `warning C4910` on MSVC - some declspec conflict
 - Look at `fwdsolver_mw.h` instantiation requirements, determine appropriate preprocessor gaurd (e.g. Clang?)
 - MATLAB QM are always complex, will crash on forward solve with real (NB interleaved API motivation)
 - MATLAB gmsh reader out of date
 - MATLAB flourescence example performance
 - Interleaved API permits shallow copy of various inputs in the MATLAB interface to reduce round-trip, exploit
 - Review element types in libfe, some contain unfinished defintions of operators and constants, remove
 - Move semantics for mathlib vectors and matrices
 - Improved initialisation for element entries
 - Check propensity for structural nonzeros viz. direct solvers
 - Python interface build assumes Release paths on Windows
 - Default link list after make mesh appears arbitrary, resulting in enormous linklist/qmvec

## Perfromance

- Bottlenecks
  - Mesh sparsity calculation heapsort (single-threaded), called when computing the system matrix for fields
  - Solvers
    - Fast direct solvers require supernodal + BLAS implementation. Use of e.g. CHOLMOD for direct solve when computing
      forward and adjoint fields for an HD problem is optimum (c. N=200k, nQM = 60).
    - MKL PARDISO less competitive than CHOLMOD.
    - Simplicial solvers such as Eigen LLT, and legacy Cholesky implementation are not competitive with iterative 
      solvers. 
    - Block Krylov methods don't appear to offer significant speedup and are reliant upon fast matrix solves thus
      indirectly require a decent BLAS.
    - Iterative solvers (CG, BICGSTAB) readily parallelised and within an order of mangnitude of direct solvers, hot path
      is Sparse-Dense Ax & Cholesky substitution. No memory issues. SpMv improvements using different CSR structures
      have shown limited improvement.
  - Jacobian computation, fast in basis, slow in mesh. Mesh path is dominated by IntFG cost, which in turn uses a
    virtual method call to element IntFG in hot loop. Experiments show 50% speedup possible by extracting this call
    and working over all RHS, and / or precomputing element integrals to avoid scaling and indexing cost.

  