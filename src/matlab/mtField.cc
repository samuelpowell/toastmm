// ========================================================================
// Implementation of class MatlabToast
// field-related methods
// ========================================================================

#include "matlabtoast.h"
#include "mexutil.h"

// Computation routines
#include "../common/calc_mesh.h"

using namespace std;

// =========================================================================

void MatlabToast::Fields (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // mesh
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");
    int i, j, idx, nQ;
    int n = mesh->nlen();

    // raster
    Raster *raster = GetBasis(prhs[1], 0, true);
    int slen = (raster ? raster->SLen() : n);

    // nodal optical parameters
    const RVector mua (n, mxGetPr (prhs[3]));
    const RVector mus (n, mxGetPr (prhs[4]));
    const RVector ref (n, mxGetPr (prhs[5]));

    // modulation frequency
    double freq = mxGetScalar (prhs[6]);

    // linear solver parameters
    char solver[128];
    double tol = 1e-10;
    mxGetString (prhs[7], solver, 128);
    if (nrhs >= 9) tol = mxGetScalar (prhs[8]);

    if(freq == 0) {

        // Real CW fields
        RCompRowMatrix qvec;
        RVector *dphi;
        RVector sphi(slen);
        CopyTMatrix (qvec, prhs[2]);
        nQ = qvec.nRows();

        CalcFields (mesh, raster, qvec, mua, mus, ref, freq, solver, tol, &dphi);

        // Map or copy to output
        mxArray *dmx = mxCreateDoubleMatrix (slen, nQ, mxREAL);
        double *pr = mxGetPr (dmx);

        for (i = idx = 0; i < nQ; i++) {
        
            if (raster) {
                raster->Map_MeshToSol (dphi[i], sphi);
            } else {
                sphi = dphi[i];
            }
            for (j = 0; j < slen; j++) {
                pr[idx] = sphi[j];
                idx++;
            }
        }
        delete []dphi;
        plhs[0] = dmx;

    } else {

        // Complex FD fields
        CCompRowMatrix qvec;
        CVector *dphi;
        CVector sphi(slen);
        CopyTMatrix (qvec, prhs[2]);
        nQ = qvec.nRows();

        CalcFields (mesh, raster, qvec, mua, mus, ref, freq, solver, tol, &dphi);
        
        // Map or copy to output
        mxArray *dmx = mxCreateDoubleMatrix (slen, nQ, mxCOMPLEX);

        #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble *pc = mxGetComplexDoubles(dmx);
        #else
        double *pr = mxGetPr (dmx);
        double *pi = mxGetPi (dmx);
        #endif

        for (i = idx = 0; i < nQ; i++) {
        if (raster) raster->Map_MeshToSol (dphi[i], sphi);
        else        sphi = dphi[i];
        for (j = 0; j < slen; j++) {
            #if MX_HAS_INTERLEAVED_COMPLEX
            pc[idx].real = sphi[j].real();
            pc[idx].imag = sphi[j].imag();
            idx++;
            #else
            pr[idx] = sphi[j].real();
            pi[idx] = sphi[j].imag();
            idx++;
            #endif
        }
        }
        delete []dphi;
        plhs[0] = dmx;

    }
}
