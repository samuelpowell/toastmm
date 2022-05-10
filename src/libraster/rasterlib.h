// ==========================================================================
// Module rasterlib
// File rasterlib.h
// General inclusion of raster lib headers
// ==========================================================================

#ifndef __RASTERLIB_H
#define __RASTERLIB_H

// General toast flags
#include "toastdef.h"

// Symbol import/export direction
#ifdef RASTERLIB_IMPLEMENTATION
#define RASTERLIB DLLEXPORT
#else
#define RASTERLIB DLLIMPORT
#endif

#include "felib.h"
#include "solution.h"
#include "mwsolution.h"
#include "raster.h"
#include "raster_px.h"
#include "raster_cp.h"
#include "raster_bl.h"
#include "raster_gb.h"
#include "raster_bb.h"
#include "raster_hb.h"
#include "raster_rb.h"
#include "raster_sb.h"
#include "raster2.h"
#include "raster_px2.h"
#include "raster_cpx.h"
#include "raster_blob2.h"
#include "raster_rb2.h"
#include "raster_bb2.h"
#include "raster_sb2.h"
#include "raster_hb2.h"
#include "raster_gb2.h"
#include "raster_cpx_tree.h"

#endif // !__RASTERLIB_H
