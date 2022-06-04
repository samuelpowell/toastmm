// ========================================================================
// Prototypes of Qvec/Mvec methods for script interfaces
// ========================================================================

#ifndef __CALC_QMVEC_H
#define __CALC_QMVEC_H

#include "toastdef.h"
#include "stoastlib.h"


void CalcQvec(QMMesh *mesh, SourceMode qtype,
              SRC_PROFILE qprof, double qwidth, CCompRowMatrix *qvec);

void CalcMvec(QMMesh *mesh, SRC_PROFILE mprof, double mwidth,
              RVector *ref, bool apply_c2a, CCompRowMatrix *mvec);

#endif
