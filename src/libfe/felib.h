// ==========================================================================
// Module felib
// File felib.h
// General inclusion of FE lib headers
// ==========================================================================

#ifndef __FELIB_H
#define __FELIB_H

// General toast flags
#include "toastdef.h"

// Symbol import/export direction
#ifdef FELIB_IMPLEMENTATION
#define FELIB DLLEXPORT
#else
#define FELIB DLLIMPORT
#endif

#include <math.h>
#include <stdlib.h>

#include "mathlib.h"
#include "point.h"
#include "ptsource.h"
#include "node.h"
#include "ndlist.h"
#include "element.h"
#include "tri3old.h"
#include "tri3.h"
#include "tri6.h"
#include "tri3D3.h"
#include "tri3D6.h"
#include "line2d2.h"
#include "tri6_ip.h"
#include "tri10.h"
#include "tri10_ip.h"
#include "tet4.h"
#include "tet10.h"
#include "tet10_ip.h"
#include "pix4.h"
#include "vox8.h"
#include "wdg6.h"
#include "ellist.h"
#include "surface.h"
#include "param.h"
#include "mesh.h"
#include "qmmesh.h"
#include "nim.h"
#include "timespec.h"
#include "refine_mesh.h"

#endif // !__FELIB_H
