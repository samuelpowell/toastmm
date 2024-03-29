// ==========================================================================
// Module stoastlib
// File stoastlib.h
// General inclusion of supertoast lib headers
// ==========================================================================

#ifndef __STOASTLIB_H
#define __STOASTLIB_H

// General toast flags
#include "toastdef.h"

// Symbol import/export direction
#ifdef STOASTLIB_IMPLEMENTATION
#define STOASTLIB DLLEXPORT
#else
#define STOASTLIB DLLIMPORT
#endif

#include "felib.h"
#include "rasterlib.h"
#include "pparse.h"
#include "fwdsolver.h"
#include "fwdsolver_mw.h"
#include "source.h"
#include "of.h"
#include "jacobian.h"
#include "regul.h"
#include "pscaler.h"
#include "lsearch.h"

#endif // !__STOASTLIB_H
