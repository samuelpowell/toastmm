// ========================================================================
// Implementation of class MatlabToast
// Qvec/Mvec related methods
// ========================================================================

#include "matlabtoast.h"
#include "mexutil.h"

// Computation routines
#include "../common/calc_qmvec.h"

using namespace std;

// ----------------------------------------------------------------------------

void MatlabToast::Qvec (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char typestr[256] = "";
    char profstr[256] = "";
    double w = 0.0;

    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    int idx, n = mesh->nlen(), nQ = mesh->nQ;

    // read source parameters from function parameters
    if (mxIsChar(prhs[1])) mxGetString (prhs[1], typestr, 256);
    if (mxIsChar(prhs[2])) mxGetString (prhs[2], profstr, 256);
    if (nrhs >= 4) {
	idx = 3;
	if (mxIsNumeric(prhs[idx]) && mxGetNumberOfElements(prhs[idx])==1){
	    w = mxGetScalar (prhs[idx]);
	    idx++;
	}
	// additional optional parameters to go here
    }

    SourceMode qtype;
    if      (!strcasecmp (typestr, "Neumann"))   qtype = SRCMODE_NEUMANN;
    else if (!strcasecmp (typestr, "Isotropic")) qtype = SRCMODE_ISOTROPIC;
    else    mexErrMsgTxt ("toastQvec: Invalid source type");

    SRC_PROFILE qprof;
    if      (!strcasecmp (profstr, "Point"))     qprof = PROF_POINT;
    else if (!strcasecmp (profstr, "Gaussian"))  qprof = PROF_GAUSSIAN;
    else if (!strcasecmp (profstr, "Cosine"))    qprof = PROF_COSINE;
    else if (!strcasecmp (profstr, "TrigBasis")) qprof = PROF_COMPLETETRIG;
    else    mexErrMsgTxt ("toastQvec: Invalid source profile");

    double qwidth;
    if (qprof != PROF_POINT) {
	if   (w > 0) qwidth = w;
	else mexErrMsgTxt ("toastQvec: Invalid source width");
    }

    CCompRowMatrix qvec;

    // build the source vectors, return as columns
    CalcQvec(mesh, qtype, qprof, qwidth, &qvec);
    CopyTMatrix (&plhs[0], qvec);
}

// =========================================================================


void MatlabToast::Mvec (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char cbuf[256];

    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    SRC_PROFILE mprof;
    mxGetString (prhs[1], cbuf, 256);
    if      (!strcasecmp (cbuf, "Gaussian"))  mprof = PROF_GAUSSIAN;
    else if (!strcasecmp (cbuf, "Cosine"))    mprof = PROF_COSINE;
    else if (!strcasecmp (cbuf, "TrigBasis")) mprof = PROF_COMPLETETRIG;
    else mexErrMsgTxt ("Invalid measurement profile");

    double mwidth = mxGetScalar (prhs[2]);

    int n, nM;

    n = mesh->nlen();
    nM = mesh->nM;
    CCompRowMatrix mvec;
    RVector ref(n);
    bool apply_c2a = true;

    if (nrhs >= 4) {
		int len = (int) (mxGetM(prhs[3])*mxGetN(prhs[3]));
		if ((len != 1 && len != n) || !mxIsDouble(prhs[3])) {
			char cbuf[256];
			sprintf (cbuf, "Mvec: parameter 3: expected double scalar or double vector of length %d", n);
			mexErrMsgTxt (cbuf);
		}
		if (len == 1) {
			double ref_homog = mxGetScalar(prhs[3]);
			if (ref_homog) ref = mxGetScalar(prhs[3]);
			else apply_c2a = false;
		} else
			CopyVector (ref, prhs[3]);
    } else {
	mexErrMsgTxt("Mvec: no refractive index values supplied");
    }
    
    // build the measurement vectors, return as matrix columns to matlab
    CalcMvec(mesh, mprof, mwidth, &ref, apply_c2a, &mvec);
    CopyTMatrix (&plhs[0], mvec);
}
