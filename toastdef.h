// -*-C++-*-
// ==========================================================================
// A list of DEFINEs which fine tune TOAST compilation options
// ==========================================================================

#ifndef __TOASTDEF_H
#define __TOASTDEF_H

#include "arch.h"

// All platforms appreciated explicit instantiation
#define NEED_EXPLICIT_INSTANTIATION 

#if defined(_WIN32) || defined(_WIN64)
#define DLLEXPORT __declspec(dllexport)
#define DLLIMPORT __declspec(dllimport)
#define DLLEXTERN extern
#else
#define DLLEXPORT
#define DLLIMPORT
#define DLLEXTERN
#endif

// #define NEED_FRIEND_PT

#define VERBOSE_LEVEL 1	// Log file verbosity level


// Threading
//
#ifdef TOAST_THREAD
#define THREAD_LEVEL 2 // 0=none, 1=fine-grain, 2=coarse-grain
#define TOAST_THREAD_MATLAB_GRADIENT  // parallelise Matlab toastGradient
#define TOAST_THREAD_MATLAB_QMVEC     // parallelise Matlab Mvec
#define TOAST_THREAD_ASSEMBLE         // parallelise system matrix assembly
#else
#define THREAD_LEVEL 0
#endif

#include "blasnames.h"

#if (!defined(_WIN32)) && (!defined(_WIN64))
//#define USE_BLAS_LEVEL1
// Use external BLAS level 1 routines (vector-vector operations)
// instead of local C++ implementations where possible.
// This is not guaranteed to improve performance, since the function
// call overhead may outweigh any savings from optimised BLAS.

#endif

//#define USE_BLAS_LEVEL2
// Use external BLAS level 2 routines (matrix-vector operations)
// instead of local C++ implementations where possible

//#define USE_BLAS_LEVEL3
// Use external BLAS level 3 routines (matrix-matrix operations)
// instead of local C++ implementations where possible

#define TRI6IP_STORE_COORDS
#define TRI10IP_STORE_COORDS
#define TET10IP_STORE_COORDS
#define TRI10_STORE_COORDS

// The number of sub-samples (in 1-D) used for numerical integration
// over elements. The total number of samples is n(n+1)/2 in a triangle, and
// ((n+3)*n+2)*n/6 in a tetrahedron, where n=NSUBSAMPLE
#define NSUBSAMPLE 50

// GCC version number
#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 \
		     + __GNUC_MINOR__ * 100 \
		     + __GNUC_PATCHLEVEL__)
#else
#define GCC_VERSION 0
#endif

#endif // !__TOASTDEF_H