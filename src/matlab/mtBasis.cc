// ========================================================================
// Implementation of class MatlabToast
// Basis-related methods
// ========================================================================

#include "matlabtoast.h"
#include "mexutil.h"

using namespace std;

// =========================================================================
// Matlab interface
// =========================================================================

void MatlabToast::SetBasis (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char optionstr[256];
    int i, basistp = 0;
    RDenseMatrix *bb = 0;
    double blobrad = 1.0;
    double blobarg = 1.0;
    double maptol = 1e-10;
    double dgscale = 0.1;
    int npad = 0;
    
    // mesh
    ASSERTARG(nrhs > 0, 0, "Too few parameters provided");
    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    // solution basis dimensions
    ASSERTARG(nrhs > 1, 1, "Too few parameters provided");
    int dim = (int) mxGetNumberOfElements (prhs[1]);
    ASSERTARG(dim == mesh->Dimension(), 3, "Invalid basis dimensions");
    IVector bdim(dim);
    for (i = 0; i < dim; i++)
	bdim[i] = (int)mxGetPr (prhs[1])[i];
    IVector gdim(bdim);

    // optional parameters
    if (nrhs > 2 && mxIsCell (prhs[2])) {
        mwSize narg = mxGetNumberOfElements (prhs[2]);
	for (i = 0; i < narg; i++) {
	    mxArray *cell = mxGetCell(prhs[2], i);
	    if (mxIsChar (cell)) { 
	        // basis type
		mxGetString (cell, optionstr, 256);
		if (!strcasecmp (optionstr, "LINEAR")) {
		    basistp = 0;
		} else if (!strcasecmp (optionstr, "CUBIC")) {
		    basistp = 1;
		} else if (!strcasecmp (optionstr, "BLOB_GAUSS")) {
		    basistp = 2;
		} else if (!strcasecmp (optionstr, "BLOB_BESSEL")) {
		    basistp = 3;
		} else if (!strcasecmp (optionstr, "BLOB_HANNING")) {
		    basistp = 4;
		} else if (!strcasecmp (optionstr, "BLOB_RAMP")) {
		    basistp = 5;
		} else if (!strcasecmp (optionstr, "BLOB_SPLINE")) {
		    basistp = 6;
		} else if (!strcasecmp (optionstr, "RADIUS")) {
		    blobrad = mxGetScalar(mxGetCell(prhs[2], ++i));
		} else if (!strcasecmp (optionstr, "SIGMA")) {
		    blobarg = mxGetScalar(mxGetCell(prhs[2], ++i));
		} else if (!strcasecmp (optionstr, "ALPHA")) {
		    blobarg = mxGetScalar(mxGetCell(prhs[2], ++i));
		} else if (!strcasecmp (optionstr, "LINEAR_V2")) {
		    basistp = 7;
		} else if (!strcasecmp (optionstr, "CONST")) {
		    basistp = 8;
		} else if (!strcasecmp (optionstr, "BLOB_RAMP2")) {
		    basistp = 9;
		} else if (!strcasecmp (optionstr, "BLOB_BESSEL2")) {
		    basistp = 10;
		} else if (!strcasecmp (optionstr, "BLOB_SPLINE2")) {
		    basistp = 11;
		} else if (!strcasecmp (optionstr, "BLOB_HANNING2")) {
		    basistp = 12;
		} else if (!strcasecmp (optionstr, "BLOB_GAUSS2")) {
		    basistp = 13;
		} else if (!strcasecmp (optionstr, "CONST_TREE")) {
		    basistp = 14;
		} else if (!strcasecmp (optionstr, "MAPTOL")) {
		    maptol = mxGetScalar(mxGetCell(prhs[2], ++i));
		} else if (!strcasecmp (optionstr, "DIAGSCALE")) {
		    dgscale = mxGetScalar(mxGetCell(prhs[2], ++i));
		} else if (!strcasecmp (optionstr, "PADDING")) {
		    npad = (int)(mxGetScalar (mxGetCell(prhs[2], ++i))+0.5);
		} else {
		    ASSERTARG(false, 1, "unexpected value");
		}
	    } else if (mxGetNumberOfElements (cell) == dim) {
	        // intermediate grid dimension
	        for (i = 0; i < dim; i++)
		    gdim[i] = (int)mxGetPr (cell)[i];
	    } else if (mxGetM (cell) == dim && mxGetN (cell) == 2) {
	        // grid bounding box
	        bb = new RDenseMatrix (2,dim);
		CopyTMatrix (*bb, cell);
	    }
	}
    }

    Raster *raster = 0;
    switch (basistp) {
    case 0:
	raster = new Raster_Pixel (bdim, gdim, mesh, bb);
	break;
    case 1:
	raster = new Raster_CubicPixel (bdim, gdim, mesh, bb);
	break;
    case 2:
	raster = new Raster_GaussBlob (bdim, gdim, mesh, blobarg, blobrad, bb);
	break;
    case 3:
	raster = new Raster_BesselBlob (bdim, gdim, mesh, blobarg, blobrad,
	    bb);
	break;
    case 4:
	raster = new Raster_HanningBlob (bdim, gdim, mesh, blobrad, bb);
	break;
    case 5:
	raster = new Raster_RampBlob (bdim, gdim, mesh, blobrad, bb);
	break;
    case 6:
	raster = new Raster_SplineBlob (bdim, gdim, mesh, blobrad, bb);
	break;
	#ifdef TOAST_FEATURE_RASTER2
    case 7:
	raster = Raster2::Create<Raster_Pixel2> (bdim, bdim, mesh, bb,
	    maptol);
	break;
    case 8:
	raster = Raster2::Create<Raster_CPixel> (bdim, bdim, mesh, bb,
            maptol);
	break;
    case 9:
	raster = Raster_Blob2::Create<Raster_Blob2_RB> (bdim, bdim, mesh,
	    blobrad, blobarg, dgscale, bb, maptol, npad);
	break;
    case 10:
	raster = Raster_Blob2::Create<Raster_Blob2_BB> (bdim, bdim, mesh,
            blobrad, blobarg, dgscale, bb, maptol, npad);
	break;
    case 11:
	raster = Raster_Blob2::Create<Raster_Blob2_SB> (bdim, bdim, mesh,
            blobrad, blobarg, dgscale, bb, maptol, npad);
	break;
    case 12:
	raster = Raster_Blob2::Create<Raster_Blob2_HB> (bdim, bdim, mesh,
	    blobrad, blobarg, dgscale, bb, maptol, npad);
	break;
    case 13:
	raster = Raster_Blob2::Create<Raster_Blob2_GB> (bdim, bdim, mesh,
	    blobrad, blobarg, dgscale, bb, maptol, npad);
	break;
    case 14:
	raster = Raster2::Create<Raster_CPixel_Tree> (bdim, bdim, mesh, bb,
            maptol);
	break;
	#endif
    }
    if (raster)
	mexLock(); // prevent mex file unloading while basis is allocated

    plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT64_CLASS, mxREAL);
    uint64_T *ptr = (uint64_T*)mxGetData (plhs[0]);
    *ptr = Ptr2Handle(raster);

    if (bb) delete bb;
}

// =========================================================================

void MatlabToast::ClearBasis (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GETBASIS_SAFE(0);
    delete raster;
    mexUnlock();
    if (verbosity >= 1)
        mexPrintf ("<Basis object deleted>\n");
}

// =========================================================================

void MatlabToast::GetBasisSize (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GETBASIS_SAFE(0);

    int i, dim = raster->Dim();
    if (nlhs >= 1) {
	RVector tmp(dim);
	for (i = 0; i < dim; i++) tmp[i] = raster->BDim()[i];
	CopyVector (&plhs[0], tmp);
    }
    if (nlhs >= 2) {
	RVector tmp(dim);
	for (i = 0; i < dim; i++) tmp[i] = raster->GDim()[i];
	CopyVector (&plhs[1], tmp);
    }
}

// =========================================================================

void MatlabToast::GetBasisNLen (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GETBASIS_SAFE(0);
    if (nlhs > 0)
        plhs[0] = mxCreateDoubleScalar (raster->mesh().nlen());
}

// =========================================================================

void MatlabToast::GetBasisBLen (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GETBASIS_SAFE(0);
    if (nlhs > 0)
        plhs[0] = mxCreateDoubleScalar (raster->BLen());
}

// =========================================================================

void MatlabToast::GetBasisSLen (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GETBASIS_SAFE(0);
    if (nlhs > 0)
        plhs[0] = mxCreateDoubleScalar (raster->SLen());
}

// =========================================================================

#ifdef TOAST_FEATURE_RASTER2

void MatlabToast::GetBasisBuu (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GETBASIS_SAFE(0);
    if (nlhs > 0) {
	Raster2 *raster2 = dynamic_cast<Raster2*>(raster);
	if (raster2) {
	    CopyMatrix (&plhs[0], *raster2->GetBuu());
	} else {
	    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	}
    }
}

// =========================================================================

void MatlabToast::GetBasisBvv (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GETBASIS_SAFE(0);
    if (nlhs > 0) {
	Raster2 *raster2 = dynamic_cast<Raster2*>(raster);
	if (raster2) {
	    CopyMatrix (&plhs[0], *raster2->GetBvv());
	} else {
	    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	}
    }
}

// =========================================================================

void MatlabToast::GetBasisBuv (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GETBASIS_SAFE(0);
    if (nlhs > 0) {
	Raster2 *raster2 = dynamic_cast<Raster2*>(raster);
	if (raster2) {
	    CopyMatrix (&plhs[0], *raster2->GetBuv());
	} else {
	    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	}
    }
}

// =========================================================================

void MatlabToast::GetBasisBvw (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GETBASIS_SAFE(0);
    IVector wdim(raster->Dim());
    for (int i = 0; i < raster->Dim(); i++)
        wdim[i] = (int)mxGetPr (prhs[1])[i];

    if (nlhs > 0) {
	Raster2 *raster2 = dynamic_cast<Raster2*>(raster);
	if (raster2) {
	    CopyMatrix (&plhs[0], *raster2->GetBvw(wdim));
	} else {
	    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	}
    }
}

// =========================================================================

void MatlabToast::GetBasisDuu (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GETBASIS_SAFE(0);
    if (nlhs > 0) {
	Raster2 *raster2 = dynamic_cast<Raster2*>(raster);
	if (raster2) {
	    CopyMatrix (&plhs[0], *raster2->GetDuu());
	} else {
	    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	}
    }
}

// =========================================================================

void MatlabToast::GetBasisDvv (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GETBASIS_SAFE(0);
    if (nlhs > 0) {
	Raster2 *raster2 = dynamic_cast<Raster2*>(raster);
	if (raster2) {
	    CopyMatrix (&plhs[0], *raster2->GetDvv());
	} else {
	    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	}
    }
}

// =========================================================================

void MatlabToast::GetBasisDuv (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GETBASIS_SAFE(0);
    if (nlhs > 0) {
	Raster2 *raster2 = dynamic_cast<Raster2*>(raster);
	if (raster2) {
	    CopyMatrix (&plhs[0], *raster2->GetDuv());
	} else {
	    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	}
    }
}

#endif // TOAST_FEATURE_RASTER2

// =========================================================================

void MatlabToast::BasisValue (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GETBASIS_SAFE(0);
    RVector pt;
    CopyVector (pt, prhs[1]);

    mwSize m = mxGetM(prhs[2]), n = mxGetN(prhs[2]);
    if (m*n == 1) { // scalar: assume this is basis index
	bool is_solidx = true;
	int idx = (int)mxGetScalar(prhs[2]) - 1;
	if (nrhs > 3 && mxIsChar(prhs[3])) {
	    char cbuf[256];
	    mxGetString(prhs[3], cbuf, 256);
	    if (!strcasecmp(cbuf,"FULLIDX"))
		is_solidx = false;
	}
	plhs[0] = mxCreateDoubleScalar(raster->Value(pt, idx, is_solidx));
    } else { // vector: assume this is coefficient vector
	RVector coeff;
	CopyVector (coeff, prhs[2]);
	bool mask = true;
	if (nrhs > 3 && mxIsChar(prhs[3])) {
	    char cbuf[256];
	    mxGetString(prhs[3], cbuf, 256);
	    if (!strcasecmp(cbuf,"NOMASK"))
		mask = false;
	}
	plhs[0] = mxCreateDoubleScalar(raster->Value(pt, coeff, mask));
    }
}

// =========================================================================

void MatlabToast::GetBasisSupportArea (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // basis handle
    Raster *raster = GETBASIS_SAFE(0);

    int idx = (int)(mxGetScalar(prhs[1])-0.5);
    RDenseMatrix sa = raster->SupportArea(idx);
    CopyMatrix (&plhs[0], sa);
}

// =========================================================================

void MatlabToast::BasisRefine (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // basis handle
    Raster *raster = GETBASIS_SAFE(0);

    RVector ridx;
    CopyVector (ridx, prhs[1]);

    int nidx = ridx.Dim();
    int *idx = new int[nidx];
    for (int i = 0; i < nidx; i++)
	idx[i] = (int)(ridx[i]-0.5);

    raster->Refine (idx, nidx);
    delete []idx;
}

// =========================================================================

void MatlabToast::MapBasis (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // basis handle
    Raster *raster = GETBASIS_SAFE(0);

    // source and target basis representations
    ASSERTARG(mxIsChar (prhs[1]), 2, "Expected string.");

    char basisstr[256], srcid, tgtid;
    mxGetString (prhs[1], basisstr, 256);
    ASSERTARG(strlen(basisstr) == 4, 2, "Format not recognised.");
    if (!strncmp (basisstr+1, "->", 2)) {
	srcid = toupper(basisstr[0]);
	tgtid = toupper(basisstr[3]);
    } else if (!strncmp (basisstr+1, "<-", 2)) {
	srcid = toupper(basisstr[3]);
	tgtid = toupper(basisstr[0]);
    } else {
	ASSERTARG(0, 2, "Format not recognised.");
    }
    
    // source field
    mwIndex m = mxGetM (prhs[2]);
    mwIndex n = mxGetN (prhs[2]);
    int nsrc = (int)(m*n);

    if (mxIsComplex (prhs[2])) {
	
	#if MX_HAS_INTERLEAVED_COMPLEX
	mxComplexDouble *pc = mxGetComplexDoubles(prhs[2]);
	#else
	double *pr = mxGetPr(prhs[2]);
	double *pi = mxGetPi(prhs[2]);
	#endif
	CVector src(nsrc), tgt;
	std::complex<double> *vptr = src.data_buffer();
	for (int i = 0; i < nsrc; i++) {
		#if MX_HAS_INTERLEAVED_COMPLEX
		vptr[i] = std::complex<double> (pc[i].real, pc[i].imag);
		#else
	    vptr[i] = std::complex<double> (pr[i], pi[i]);
		#endif
		}
	switch (srcid) {
	case 'M':
	    ASSERTARG(raster->mesh().nlen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'G':
		raster->Map_MeshToGrid (src, tgt);
		break;
	    case 'B':
		raster->Map_MeshToBasis (src, tgt);
		break;
	    case 'S':
		raster->Map_MeshToSol (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	case 'G':
	    ASSERTARG(raster->GLen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'M':
		raster->Map_GridToMesh (src, tgt);
		break;
	    case 'B':
		raster->Map_GridToBasis (src, tgt);
		break;
	    case 'S':
		raster->Map_GridToSol (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	case 'B':
	    ASSERTARG(raster->BLen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'G':
		raster->Map_BasisToGrid (src, tgt);
		break;
	    case 'S':
		raster->Map_BasisToSol (src, tgt);
		break;
	    case 'M':
		raster->Map_BasisToMesh (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	case 'S':
	    ASSERTARG(raster->SLen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'B':
		raster->Map_SolToBasis (src, tgt);
		break;
	    case 'G':
		raster->Map_SolToGrid (src, tgt);
		break;
	    case 'M':
		raster->Map_SolToMesh (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	default:
	    ASSERTARG(0, 2, "Source id not recognised.");
	    return;
	}

	CopyVector (&plhs[0], tgt);

    } else {

	RVector src (nsrc, mxGetPr (prhs[2]));
	RVector tgt;

	switch (srcid) {
	case 'M':
	    ASSERTARG(raster->mesh().nlen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'G':
		raster->Map_MeshToGrid (src, tgt);
		break;
	    case 'B':
		raster->Map_MeshToBasis (src, tgt);
		break;
	    case 'S':
		raster->Map_MeshToSol (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	case 'G':
	    ASSERTARG(raster->GLen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'M':
		raster->Map_GridToMesh (src, tgt);
		break;
	    case 'B':
		raster->Map_GridToBasis (src, tgt);
		break;
	    case 'S':
		raster->Map_GridToSol (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	case 'B':
	    ASSERTARG(raster->BLen() == nsrc, 3,
	        "Source vector unexpected length.");
	    switch (tgtid) {
	    case 'G':
		raster->Map_BasisToGrid (src, tgt);
		break;
	    case 'S':
		raster->Map_BasisToSol (src, tgt);
		break;
	    case 'M':
		raster->Map_BasisToMesh (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	case 'S':
	    ASSERTARG(raster->SLen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'B':
		raster->Map_SolToBasis (src, tgt);
		break;
	    case 'G':
		raster->Map_SolToGrid (src, tgt);
		break;
	    case 'M':
		raster->Map_SolToMesh (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	default:
	    ASSERTARG(0, 2, "Source id not recognised.");
	    return;
	}

	CopyVector (&plhs[0], tgt);
    }    
}

// =========================================================================

void MatlabToast::SolutionMask (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");
    int slen = raster->SLen();

    plhs[0] = mxCreateDoubleMatrix (1, slen, mxREAL);
    double *pr = mxGetPr (plhs[0]);

    for (int i = 0; i < slen; i++)
	pr[i] = raster->Sol2Basis (i) + 1;
        // "+1" to adjust to matlab's 1-based array indices
}

// =========================================================================

void MatlabToast::GridElref (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    int i, n = raster->GLen();
    int *elref = raster->Elref();

    plhs[0] = mxCreateDoubleMatrix (n,1,mxREAL);
    double *pr = mxGetPr(plhs[0]);
    for (i = 0; i < n; i++) {
      pr[i] = elref[i] + 1;
        // make 1-based. Gridpoints outside mesh support are set to zero
    }    
}

// =========================================================================

void MatlabToast::SampleBasis (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, n;

    // raster
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    // for now, assume real-valued argument
    RVector svec;
    CopyVector (svec, prhs[1]);

    RVector grd;
    CopyVector (grd, prhs[2]);
    IVector igrd(grd.Dim());
    for (i = 0, n = 1; i < grd.Dim(); i++) {
	igrd[i] = (int)(grd[i]+0.5);
	n *= igrd[i];
    }
    RVector img(n);
    raster->Sample(svec,igrd,img);
    CopyVector (&plhs[0],img);
}

void MatlabToast::ImageGradient (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // raster
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    const IVector &gdim = raster->GDim();
    const RVector &gsize = raster->GSize();
    const int *elref = raster->Elref();
    int glen = raster->GLen();
    int d, dim = gdim.Dim();

    if (mxIsComplex (prhs[1])) {
	CVector img;
	CVector *grad = new CVector[dim];
	for (d = 0; d < dim; d++) grad[d].New (glen);
	CopyVector (img, prhs[1]);
	::ImageGradient (gdim, gsize, img, grad, elref);
	CDenseMatrix g (dim, glen);
	for (d = 0; d < dim; d++) g.SetRow (d, grad[d]);
	delete []grad;
	CopyMatrix (&plhs[0], g);
    } else {
	RVector img;
	RVector *grad = new RVector[dim];
	for (d = 0; d < dim; d++) grad[d].New (glen);
	CopyVector (img, prhs[1]);
	::ImageGradient (gdim, gsize, img, grad, elref);
	RDenseMatrix g (dim, glen);
	for (d = 0; d < dim; d++) g.SetRow (d, grad[d]);
	delete []grad;
	CopyMatrix (&plhs[0], g);
    }
}
