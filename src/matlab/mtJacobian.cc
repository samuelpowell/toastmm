// ========================================================================
// Implementation of class MatlabToast
// Jacobian-related methods
// ========================================================================

#include "matlabtoast.h"
#include "mexutil.h"

// Computation routines
#include "../common/calc_jacobian.h"

using namespace std;

unsigned int verb;


// =========================================================================
// Matlab interface
// =========================================================================

// =========================================================================
// Frequency-domain Jacobian

void MatlabToast::Jacobian (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    verb = verbosity;

    // mesh
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");

    int n   = mesh->nlen();
    int nq  = mesh->nQ;
    int nm  = mesh->nM;
    int nqm = mesh->nQM;

    // raster
    Raster *raster = GetBasis(prhs[1], 0, true);

    // Output 
    RDenseMatrix J;

    if (nrhs == 5) {

	// this is the version that provides fields and projections directly
	int i, j;

	#if MX_HAS_INTERLEAVED_COMPLEX
    mxComplexDouble *pc;
    #else
    double *pr, *pi;
    #endif

	// copy fields
	const mxArray *mx_dphi = prhs[2];
	ASSERTARG(mxGetM(mx_dphi) == n, 3, "Unexpected number of rows");
	ASSERTARG(mxGetN(mx_dphi) == nq, 3, "Unexpected number of columns");
    ASSERTARG(mxIsComplex (mx_dphi), 3, "Must be complex");
    #if MX_HAS_INTERLEAVED_COMPLEX
    pc = mxGetComplexDoubles(mx_dphi);
    #else
	pr  = mxGetPr (mx_dphi);
	pi  = mxGetPi (mx_dphi);
    #endif
	CVector *dphi = new CVector[nq];
	for (i = 0; i < nq; i++) {
	    dphi[i].New (n);
	    std::complex<double> *v = dphi[i].data_buffer();
	    for (j = 0; j < n; j++) {
            #if MX_HAS_INTERLEAVED_COMPLEX
            *v++ = std::complex<double> ((*pc).real, (*pc).imag);
            pc++;
            #else
	        *v++ = std::complex<double> (*pr++, *pi++);
            #endif
        }
	}
	// copy adjoint fields
	const mxArray *mx_aphi = prhs[3];
	ASSERTARG(mxGetM(mx_aphi) == n, 4, "Unexpected number of rows");
	ASSERTARG(mxGetN(mx_aphi) == nm, 4, "Unexpected number of columns");
	ASSERTARG(mxIsComplex (mx_aphi), 4, "Must be complex");
    #if MX_HAS_INTERLEAVED_COMPLEX
    pc = mxGetComplexDoubles(mx_aphi);
    #else
	pr  = mxGetPr (mx_aphi);
	pi  = mxGetPi (mx_aphi);
    #endif
	CVector *aphi = new CVector[nm];
	for (i = 0; i < nm; i++) {
	    aphi[i].New (n);
	    std::complex<double> *v = aphi[i].data_buffer();
	    for (j = 0; j < n; j++) {
            #if MX_HAS_INTERLEAVED_COMPLEX
            *v++ = std::complex<double> ((*pc).real, (*pc).imag);
            pc++;
            #else
	        *v++ = std::complex<double> (*pr++, *pi++);
            #endif
        }
	}
	// copy projections
	const mxArray *mx_proj = prhs[4];
	ASSERTARG(mxGetM(mx_proj)*mxGetN(mx_proj) == nqm, 5,"Unexpected size");
	ASSERTARG(mxIsComplex(mx_proj), 5, "Must be complex");
	CVector proj(nqm);    
    #if MX_HAS_INTERLEAVED_COMPLEX
    pc = mxGetComplexDoubles(mx_proj);
    #else
	pr  = mxGetPr (mx_proj);
	pi  = mxGetPi (mx_proj);
    #endif
	std::complex<double> *v = proj.data_buffer();
	for (i = 0; i < nqm; i++) {
        #if MX_HAS_INTERLEAVED_COMPLEX
        *v++ = std::complex<double> ((*pc).real, (*pc).imag);
        pc++;
        #else
        *v++ = std::complex<double> (*pr++, *pi++);
        #endif
    }

	CalcJacobian (mesh, raster, dphi, aphi, &proj, DATA_LOG, J);              

    } else {

	// this is the version that calculates fields on the fly

	// source vectors
	CCompRowMatrix qvec;
	CopyTMatrix (qvec, prhs[2]);

	// measurement vectors
	CCompRowMatrix mvec;
	CopyTMatrix (mvec, prhs[3]);

    // bounds check on optical parameters
    ASSERTARG(mxGetM(prhs[4]) == n, 5, "Unexpected size");
    ASSERTARG(mxGetM(prhs[5]) == n, 6, "Unexpected size");
    ASSERTARG(mxGetM(prhs[6]) == n, 7, "Unexpected size");
        
	// nodal optical parameters
	RVector mua (n, mxGetPr (prhs[4]));
	RVector mus (n, mxGetPr (prhs[5]));
	RVector ref (n, mxGetPr (prhs[6]));

	// modulation frequency
	double freq = mxGetScalar (prhs[7]);

	// linear solver parameters
	char solver[128];
	double tol = 1e-10;
	mxGetString (prhs[8], solver, 128);
	if (nrhs >= 10) tol = mxGetScalar (prhs[9]);
	
	CalcJacobian (mesh, raster, qvec, mvec, mua, mus, ref, freq, solver, tol, J);
    }    

    CopyMatrix (&plhs[0], J);
}

// ==========================================================================
// CW Jacobian

void MatlabToast::JacobianCW (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    verb = verbosity;

    // mesh
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");

    int n   = mesh->nlen();
    int nq  = mesh->nQ;
    int nm  = mesh->nM;
    int nqm = mesh->nQM;

    // raster
    Raster *raster = GetBasis(prhs[1], 0, true);
    
    // Output 
    RDenseMatrix J;

    if (nrhs == 5) {

    // compute jacobian from fields and projection data       
	int i, j;
	double *pr;
	
	// copy fields
	const mxArray *mx_dphi = prhs[2];
	ASSERTARG(mxGetM(mx_dphi) == n, 3, "Unexpected number of rows");
	ASSERTARG(mxGetN(mx_dphi) == nq, 3, "Unexpected number of columns");
	pr  = mxGetPr (mx_dphi);
	RVector *dphi = new RVector[nq];
	for (i = 0; i < nq; i++) {
	    dphi[i].New (n);
	    double *v = dphi[i].data_buffer();
	    for (j = 0; j < n; j++) {
	        *v++ = *pr++;
        }
	}
	// copy adjoint fields
	const mxArray *mx_aphi = prhs[3];
	ASSERTARG(mxGetM(mx_aphi) == n, 4, "Unexpected number of rows");
	ASSERTARG(mxGetN(mx_aphi) == nm, 4, "Unexpected number of columns");
	pr = mxGetPr (mx_aphi);
	RVector *aphi = new RVector[nm];
	for (i = 0; i < nm; i++) {
	    aphi[i].New (n);
	    double *v = aphi[i].data_buffer();
	    for (j = 0; j < n; j++) {
	        *v++ = *pr++;
        }
	}
	// copy projections
	const mxArray *mx_proj = prhs[4];
	ASSERTARG(mxGetM(mx_proj)*mxGetN(mx_proj) == nqm, 5,"Unexpected size");
	RVector proj(nqm);
	pr = mxGetPr (mx_proj);
	double *v = proj.data_buffer();
	for (i = 0; i < nqm; i++)
	    *v++ = *pr++;

    CalcJacobianCW (mesh, raster, dphi, aphi, &proj, DATA_LOG, J);
    
    } else {

        // source vectors
        RCompRowMatrix qvec;
        CopyTMatrix (qvec, prhs[2]);

        // measurement vectors
        RCompRowMatrix mvec;
        CopyTMatrix (mvec, prhs[3]);
        
        // bounds check on optical parameters
        ASSERTARG(mxGetM(prhs[4]) == n, 5, "Unexpected size");
        ASSERTARG(mxGetM(prhs[5]) == n, 6, "Unexpected size");
        ASSERTARG(mxGetM(prhs[6]) == n, 7, "Unexpected size");

        // nodal optical parameters
        RVector mua (n, mxGetPr (prhs[4]));
        RVector mus (n, mxGetPr (prhs[5]));
        RVector ref (n, mxGetPr (prhs[6]));
        
        
        // linear solver parameters
        char solver[128];
        double tol = 1e-10;
        mxGetString (prhs[7], solver, 128);
        if (nrhs >= 9) tol = mxGetScalar (prhs[8]);
        
        CalcJacobianCW (mesh, raster, qvec, mvec, mua, mus, ref, solver, tol, J);
    }

    CopyMatrix (&plhs[0], J);
}


// ==========================================================================
// Implementation: see common module
// ==========================================================================
