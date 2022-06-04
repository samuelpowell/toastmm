// ========================================================================
// Prototypes of Jacobian methods for script interfaces
// ========================================================================

#include "stoastlib.h"

void CalcJacobian (QMMesh *mesh, Raster *raster,
    const CVector *dphi, const CVector *aphi,
    const CVector *proj, DataScale dscale, RDenseMatrix &J);

void CalcJacobian (QMMesh *mesh, Raster *raster,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    const RVector &mua, const RVector &mus, const RVector &ref,
    double freq, char *solver, double tol, RDenseMatrix &J);

void CalcJacobianCW (QMMesh *mesh, Raster *raster,
    const RVector *dphi, const RVector *aphi,
    const RVector *proj, DataScale dscale, RDenseMatrix &J);

void CalcJacobianCW (QMMesh *mesh, Raster *raster,
    const RCompRowMatrix &qvec, const RCompRowMatrix &mvec,
    const RVector &mua, const RVector &mus, const RVector &ref,
    char *solver, double tol, RDenseMatrix &J);