// ========================================================================
// Implementation of Jacobian methods for script interfaces
// ========================================================================

#include "toastdef.h"
#include "calc_jacobian.h"

// ==========================================================================
// Calculate CW Jacobian from given optical parameters.
// This version calculates the direct and adjoint fields, and boundary
// projection data on the fly from the provided optical parameters

void CalcJacobianCW (QMMesh *mesh, Raster *raster,
    const RCompRowMatrix &qvec, const RCompRowMatrix &mvec,
    const RVector &mua, const RVector &mus, const RVector &ref,
    char *solver, double tol, RDenseMatrix &J)
{
    const double c0 = 0.3;
    int i, n, dim, nQ, nM, slen;

    n    = mesh->nlen();
    dim  = mesh->Dimension();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    slen = (raster ? raster->SLen() : n);

    RVector *dphi, *aphi;
    RFwdSolver FWS (mesh, solver, tol);

    // Solution in mesh basis
    Solution msol(OT_NPARAM, n);
    msol.SetActive (OT_CMUA, true);
    msol.SetActive (OT_CKAPPA, true);
    

    // Set optical coefficients
    msol.SetParam (OT_CMUA, mua*c0/ref);
    msol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    RVector c2a(n);
    for (i = 0; i < n; i++)
	c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));
    msol.SetParam (OT_C2A, c2a);

    // build the field vectors
    dphi = new RVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New (n);
    aphi = new RVector[nM];
    for (i = 0; i < nM; i++) aphi[i].New (n);

    // Calculate direct and adjoint fields
    FWS.Allocate ();
    FWS.Reset (msol, 0);
    FWS.CalcFields (qvec, dphi);
    FWS.CalcFields (mvec, aphi);

    // Calculate Jacobian
    int ndat = mesh->nQM;
    int nprm = slen;
    J.New(ndat,nprm);
    GenerateJacobian_cw (raster, mesh, mvec, dphi, aphi, DATA_LOG, &J);

    delete []dphi;
    delete []aphi;
}                                                                              

// ==========================================================================
// Calculate CW Jacobian from given direct and adjoint fields and boundary
// projection data


void CalcJacobianCW (QMMesh *mesh, Raster *raster,
    const RVector *dphi, const RVector *aphi,
    const RVector *proj, DataScale dscale, RDenseMatrix &J)
{
    int nQM  = mesh->nQM;
    int slen = (raster ? raster->SLen() : mesh->nlen());
    int ndat = nQM;
    int nprm = slen;

    J.New(ndat,nprm);
    GenerateJacobian_cw (raster, mesh, dphi, aphi, proj, DATA_LOG, &J);
}



// ==========================================================================
// Calculate Jacobian from given direct and adjoint fields and boundary
// projection data

void CalcJacobian (QMMesh *mesh, Raster *raster,
    const CVector *dphi, const CVector *aphi,
    const CVector *proj, DataScale dscale, RDenseMatrix &J)
{
    int nQM, slen, ndat, nprm;
    nQM  = mesh->nQM;
    slen = (raster ? raster->SLen() : mesh->nlen());
    ndat = nQM * 2;
    nprm = slen * 2;

    J.New(ndat,nprm);
    GenerateJacobian (raster, mesh, dphi, aphi, proj, dscale, J);
}

// ==========================================================================
// Calculate Jacobian from given optical parameters.
// This version calculates the direct and adjoint fields, and boundary
// projection data on the fly from the provided optical parameters

void CalcJacobian (QMMesh *mesh, Raster *raster,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    const RVector &mua, const RVector &mus, const RVector &ref,
    double freq, char *solver, double tol, RDenseMatrix &J)
{
    const double c0 = 0.3;
    int i, n, dim, nQ, nM, nQM, slen;

    n    = mesh->nlen();
    dim  = mesh->Dimension();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    nQM  = mesh->nQM;
    slen = (raster ? raster->SLen() : n);

    CVector *dphi, *aphi;
    CFwdSolver FWS (mesh, solver, tol);

    // Solution in mesh basis
    Solution msol(OT_NPARAM, n);
    msol.SetActive (OT_CMUA, true);
    msol.SetActive (OT_CKAPPA, true);

    // Set optical coefficients
    msol.SetParam (OT_CMUA, mua*c0/ref);
    msol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    RVector c2a(n);
    for (i = 0; i < n; i++)
	c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));
    msol.SetParam (OT_C2A, c2a);

    FWS.SetDataScaling (DATA_LOG);

    double omega = freq * 2.0*Pi*1e-6; // convert from MHz to rad

    // build the field vectors
    dphi = new CVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New (n);
    aphi = new CVector[nM];
    for (i = 0; i < nM; i++) aphi[i].New (n);

    // Calculate direct and adjoint fields
    FWS.Allocate ();
    FWS.Reset (msol, omega);
    FWS.CalcFields (qvec, dphi);
    FWS.CalcFields (mvec, aphi);

    // Calculate projections if required
    CVector *proj = 0;
    DataScale dscale = FWS.GetDataScaling();
    if (dscale == DATA_LOG) {
    	proj = new CVector(nQM);
	    *proj = FWS.ProjectAll (mvec, dphi, DATA_LIN);
    }

    // Calculate Jacobian
    CalcJacobian (mesh, raster, dphi, aphi, proj, dscale, J);

    delete []dphi;
    delete []aphi;
    if (proj) delete proj;
}                                                                              

