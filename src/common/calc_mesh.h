// ========================================================================
// Prototypes of Mesh methods for script interfaces
// ========================================================================

#ifndef __CALC_MESH_H
#define __CALC_MESH_H

#include "toastdef.h"
#include "stoastlib.h"

// =========================================================================
// Mesh
// =========================================================================

void BuildMesh(QMMesh **mesh, int nvtx, int nel, int dim, int nnd0,
               double *vtx, int *idx, int *eltp);

// =========================================================================
// Fields
// =========================================================================

// Frequency domain (complex)
void CalcFields (QMMesh *mesh, Raster *raster,
    const CCompRowMatrix &qvec, const RVector &mua, const RVector &mus,
    const RVector &ref, double freq, char *solver, double tol,
    CVector **dphi);

// Continuous wave (real)
void CalcFields (QMMesh *mesh, Raster *raster,
    const RCompRowMatrix &qvec, const RVector &mua, const RVector &mus,
    const RVector &ref, double freq, char *solver, double tol,
    RVector **dphi);    

#endif
