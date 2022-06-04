// ========================================================================
// Implementation of class MatlabToast
// Gradient-related methods
// ========================================================================

#include "matlabtoast.h"
#include "mexutil.h"

// Computation routines
#include "../common/calc_gradient.h"

#ifdef TOAST_THREAD
#include "task.h"
#endif

using namespace std;

// =========================================================================

void MatlabToast::Gradient (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // See if we need to do a real or complex computation
    double freq = mxGetScalar (prhs[7]);

    if (!freq) GradientReal (nlhs, plhs, nrhs, prhs);
    else       GradientCplx (nlhs, plhs, nrhs, prhs);
}

// =========================================================================

void MatlabToast::GradientReal (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j;

    // mesh
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");

    int n = mesh->nlen();

    Raster *raster = GetBasis(prhs[1]);
    ASSERTARG(raster, 2, "Basis not found");

    RVector grad(raster->SLen()*2);

    // source vectors
    RCompRowMatrix qvec;
    CopyTMatrix (qvec, prhs[2]);

    // measurement vectors
    RCompRowMatrix mvec;
    CopyTMatrix (mvec, prhs[3]);

    // nodal optical parameters
    RVector mua (n, mxGetPr (prhs[4]));
    RVector mus (n, mxGetPr (prhs[5]));
    RVector ref (n, mxGetPr (prhs[6]));

    // data & sd vectors (truncate phase parts if present)
    RVector data (mesh->nQM, mxGetPr (prhs[8]));
    RVector sd (mesh->nQM, mxGetPr (prhs[9]));
    
    char solver[128] = "direct";
    double tol = 1e-10;
    RVector *phi = 0;
    RVector *proj = 0;

    // read optional parameters as tupels
    for (i = 10; i < nrhs; i++) {
        char label[128];
        mxGetString (prhs[i], label, 128);

	if (!strcasecmp(label, "Method")) {     // linear solver method
	    mxGetString (prhs[++i], solver, 128);
	} else if (!strcasecmp(label, "Tolerance")) { // linear solver tolerance
  	    tol = mxGetScalar (prhs[++i]);
	} else if (!strcasecmp(label, "Fields")) { // fields for all sources
	    i++;
	    int q, nq = mesh->nQ, nlen = mesh->nlen();
	    mwSize m = mxGetM(prhs[i]);
	    mwSize n = mxGetN(prhs[i]);
	    double *pr = mxGetPr(prhs[i]);
	    xASSERT(m == nlen && n == nq, "Parameter phi wrong dimension");
	    phi = new RVector[nq];
	    for (q = 0; q < nq; q++) {
	        phi[q].New (nlen);
		for (j = 0; j < nlen; j++)
		    phi[q][j] = *pr++;
	    }
	} else if (!strcasecmp(label,"Projections")) {
	    proj = new RVector(mesh->nQM, mxGetPr (prhs[++i]));
	} else if (!strcasecmp(label,"Unwrap")) {
	    // N/A
	} else {
	    mexErrMsgTxt("Error parsing arguments");
	}
    }

    RFwdSolver FWS (mesh, solver, tol);
    FWS.SetDataScaling (DATA_LOG);

    GetGradientReal (mesh, raster, FWS, mua, mus, ref, data, sd,
        qvec, mvec, grad, phi, proj);
    CopyVector (&plhs[0], grad);

    if (phi) delete []phi;
    if (proj) delete proj;
}

// =========================================================================

void MatlabToast::GradientCplx (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j;

    // mesh
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");

    int n = mesh->nlen();

    Raster *raster = GetBasis(prhs[1]);
    ASSERTARG(raster, 2, "Basis not found");

    RVector grad(raster->SLen()*2);

    // source vectors
    CCompRowMatrix qvec;
    CopyTMatrix (qvec, prhs[2]);

    // measurement vectors
    CCompRowMatrix mvec;
    CopyTMatrix (mvec, prhs[3]);

    // nodal optical parameters
    RVector mua (n, mxGetPr (prhs[4]));
    RVector mus (n, mxGetPr (prhs[5]));
    RVector ref (n, mxGetPr (prhs[6]));

    // modulation frequency
    double freq = mxGetScalar (prhs[7]);
    
    // data & sd vectors
    RVector data (mesh->nQM*2, mxGetPr (prhs[8]));
    RVector sd (mesh->nQM*2, mxGetPr (prhs[9]));
    
    char solver[128] = "direct";
    double tol = 1e-10;
    bool unwrap = false;
    CVector *phi = 0;
    RVector *proj = 0;

    // read optional parameters as tupels
    for (i = 10; i < nrhs; i++) {
        char label[128];
        mxGetString (prhs[i], label, 128);

	if (!strcasecmp(label, "Method")) {     // linear solver method
	    mxGetString (prhs[++i], solver, 128);
	} else if (!strcasecmp(label, "Tolerance")) { // linear solver tolerance
  	    tol = mxGetScalar (prhs[++i]);
	} else if (!strcasecmp(label, "Fields")) { // fields for all sources
	    i++;
	    int q, nq = mesh->nQ, nlen = mesh->nlen();
	    mwSize m = mxGetM(prhs[i]);
	    mwSize n = mxGetN(prhs[i]);

        #if MX_HAS_INTERLEAVED_COMPLEX
	    mxComplexDouble *pc = mxGetComplexDoubles(prhs[1]);
	    #else
	    double *pr = mxGetPr(prhs[i]);
	    double *pi = mxGetPi(prhs[i]);
        #endif

	    xASSERT(m == nlen && n == nq, "Parameter phi wrong dimension");
	    phi = new CVector[nq];
	    for (q = 0; q < nq; q++) {
	        phi[q].New (nlen);
		for (j = 0; j < nlen; j++)  {
            #if MX_HAS_INTERLEAVED_COMPLEX
            phi[q][j] = std::complex<double> ((*pc).real, (*pc).imag);
            pc++;
            #else
		    phi[q][j] = std::complex<double> (*pr++, *pi++);
            #endif
        }
	    }
	} else if (!strcasecmp(label,"Projections")) {
	    proj = new RVector(mesh->nQM*2, mxGetPr (prhs[++i]));
	} else if (!strcasecmp(label,"Unwrap")) {
	    unwrap = mxIsLogicalScalarTrue (prhs[++i]);	    
	} else {
	    mexErrMsgTxt("Error parsing arguments");
	}
    }

	int nth;
#ifdef TOAST_THREAD
	nth = Task::GetThreadCount();
#else
	nth = 1;
#endif

    CFwdSolver FWS (mesh, solver, tol, nth);
    FWS.SetDataScaling (DATA_LOG);
    FWS.SetPhaseUnwrap (unwrap);

    GetGradientCplx (mesh, raster, FWS, mua, mus, ref, freq, data, sd,
		 qvec, mvec, grad, phi, proj);
    CopyVector (&plhs[0], grad);

    if (phi) delete []phi;
    if (proj) delete proj;
}
