// -*-C++-*-
// ==========================================================================
// Module mathlib
// File mathlib.h
// General maths declarations and inclusion of maths headers
// ==========================================================================

#ifndef __MATHLIB_H
#define __MATHLIB_H

#include <math.h>
#include <stdlib.h>

// General toast flags
#include "toastdef.h"

// Symbol import/export direction
#ifdef MATHLIB_IMPLEMENTATION
#define MATHLIB DLLEXPORT
#else
#define MATHLIB DLLIMPORT
#endif

#ifndef NO_MATH_HEADERS
#include "util.h"
#include "mathdef.h"
#include "error.h"
#include "nr.h"

#include <complex>

#include "vector.h"
#include "matrix.h"
#include "dgmatrix.h"
#include "dnsmatrix.h"
#include "symatrix.h"
#include "gsmatrix.h"
#include "cdmatrix.h"
#include "crmatrix.h"
#include "crmatrix_cm.h"
#include "scrmatrix.h"
#include "cr_cholesky.h"
#include "precon.h"
#include "gmres.h"
#include "fourn.h" // only needed if using C version rather than C++

#ifdef TOAST_THREAD
#include "task.h"
#endif

#ifdef COMPUTE_FLOPS
extern unsigned int flops_add;
extern unsigned int flops_mul;
#endif

#endif // !NO_MATH_HEADERS

#ifndef SQR
#define SQR(_X) ((_X)*(_X))
#endif
#ifndef CUBE
#define CUBE(_X) ((_X)*(_X)*(_X))
#endif

// The template implementation headers
#include "vector_imp.hpp"
#include "gmres_imp.hpp"
#include "precon_imp.hpp"
#include "matrix_imp.hpp"
#include "dnsmatrix_imp.hpp"
#include "symatrix_imp.hpp"
#include "gsmatrix_imp.hpp"
#include "dgmatrix_imp.hpp"
#include "crmatrix_imp.hpp"
#include "cdmatrix_imp.hpp"

#endif // !__MATHLIB_H
