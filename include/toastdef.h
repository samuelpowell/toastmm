// -*-C++-*-
// ==========================================================================
// A list of DEFINEs which fine tune TOAST compilation options
// ==========================================================================

#ifndef __TOASTDEF_H
#define __TOASTDEF_H

#include "arch.h"
#include "blasnames.h"

// All platforms appreciated explicit instantiation
//
#define NEED_EXPLICIT_INSTANTIATION 

// Windows DLL import/export declerations
//
#if (defined(_WIN32) || defined(_WIN64)) && !defined(TOAST_STATIC)
#define DLLEXPORT __declspec(dllexport)
#define DLLIMPORT __declspec(dllimport)
#define DLLEXTERN extern
#else
#define DLLEXPORT
#define DLLIMPORT
#define DLLEXTERN
#endif

// #define NEED_FRIEND_PT

// Verbosity
//
#define VERBOSE_LEVEL 1	// Log file verbosity level

// Threading
//
#ifdef TOAST_THREAD
#define TOAST_THREAD_ASSEMBLE         	// parallelise system matrix assembly
#else
#define THREAD_LEVEL 0
#endif

// BLAS
//
// Use external BLAS level 1 routines (vector-vector operations) instead of local C++ 
// implementations where possible. This is not guaranteed to improve performance, since the
// function call overhead may outweigh any savings from optimised BLAS.
//
//#define USE_BLAS_LEVEL1

// Use external BLAS level 2 routines (matrix-vector operations) instead of local C++
// implementations where possible
//
//#define USE_BLAS_LEVEL2

// Use external BLAS level 3 routines (matrix-matrix operations) instead of local C++
// implementations where possible
//
//#define USE_BLAS_LEVEL3

// libfe options
//
// The number of sub-samples (in 1-D) used for numerical integration over elements. The
// total number of samples is n(n+1)/2 in a triangle, and ((n+3)*n+2)*n/6 in a tetrahedron,
// where n=NSUBSAMPLE
//
#define NSUBSAMPLE 50

#define TRI6IP_STORE_COORDS
#define TRI10IP_STORE_COORDS
#define TET10IP_STORE_COORDS
#define TRI10_STORE_COORDS

#endif // !__TOASTDEF_H