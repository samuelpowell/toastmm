
// =========================================================================
// Prototypes of gradient functions
// =========================================================================

#ifndef __CALC_GRADIENT_H
#define __CALC_GRADIENT_H

#include "toastdef.h"
#include "stoastlib.h"

void AddDataGradientReal (QMMesh *mesh, Raster *raster, const RFwdSolver &FWS,
    const RVector &proj, const RVector &data, const RVector &sd, RVector *dphi,
    const RCompRowMatrix &mvec, RVector &grad);

void AddDataGradientCplx (QMMesh *mesh, Raster *raster, const CFwdSolver &FWS,
    const RVector &proj, const RVector &data, const RVector &sd, CVector *dphi,
    const CCompRowMatrix &mvec, RVector &grad);

void GetGradientReal (QMMesh *mesh, Raster *raster, RFwdSolver &FWS,
    const RVector &mua, const RVector &mus, const RVector &ref,
    const RVector &data, const RVector &sd,
    const RCompRowMatrix &qvec, const RCompRowMatrix &mvec,
    RVector &grad, RVector *phi=0, RVector *proj=0);

void GetGradientCplx (QMMesh *mesh, Raster *raster, CFwdSolver &FWS,
    const RVector &mua, const RVector &mus, const RVector &ref, double freq,
    const RVector &data, const RVector &sd,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    RVector &grad, CVector *phi=0, RVector *proj=0);

#endif