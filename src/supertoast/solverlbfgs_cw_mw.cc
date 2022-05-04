// ==========================================================================
// SolverLBFGS: LBFGS method (Limited-memory Broyden-Fletcher-Goldfarb-Shanno)
// This version interfaces to the external liblbfgs library
// ==========================================================================

#include "stoastlib.h"
#include "util.h"
#include "mwsolution.h"
#include "solverlbfgs_cw_mw.h"
#include "lbfgs.h"
#include "timing.h"

using namespace std;

// ==========================================================================
// external references

extern int g_imgfmt;
extern char g_prefix[256];
extern double clock0;

// ==========================================================================
// local prototypes

static void MW_get_gradient (const Raster &raster, RFwdSolverMW &FWS,
		      RVector *dphi,
		      const RCompRowMatrix &qvec, const RCompRowMatrix &mvec,
		      const MWsolution *msol, RVector &grad,
		      const RVector &data, const RVector &sd);

// ==========================================================================
// Data structure to provide context for LBFGS library callback functions

struct LBFGS_DATA {
    const ObjectiveFunction *of;
    const Raster *raster;
    RFwdSolverMW *fws;
    Regularisation *reg;
    const RCompRowMatrix *qvec, *mvec;
    const Scaler *pscaler;
    MWsolution *msol;
    Solution *bsol;
    const RVector *data;
    const RVector *sd;
};

// ==========================================================================

SolverLBFGS_CW_MW::SolverLBFGS_CW_MW (ParamParser *_pp): Solver_CW (_pp)
{
    itmax = 50;
    epsilon = 1e-6;
    delta = 1e-5;
    history = 5;
}

// ==========================================================================
// Callback function for LBFGS library solver:
// Evaluate objective function and gradient

static lbfgsfloatval_t evaluate (void *instance,
    const lbfgsfloatval_t *x_lbfgs,
    lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{
    // Calculate objective function and gradient for given solution x
    LBFGS_DATA *ldata = (LBFGS_DATA*)instance;
    const ObjectiveFunction *of = ldata->of;
    const Raster *raster = ldata->raster;
    RFwdSolverMW *fws = ldata->fws;
    Regularisation *reg = ldata->reg;
    const RCompRowMatrix *qvec = ldata->qvec;
    const RCompRowMatrix *mvec = ldata->mvec;
    const Scaler *pscaler = ldata->pscaler;
    MWsolution *msol = ldata->msol;
    Solution *bsol = ldata->bsol;
    const RVector *data = ldata->data;
    const RVector *sd = ldata->sd;
    const QMMesh *mesh = fws->MeshPtr();
    int i, nlen = mesh->nlen();
    int ndat = data->Dim();
    int nofwavel = msol->nofwavel;

    RVector x(n, (double*)x_lbfgs);
    RVector r(n);
    RVector proj(ndat);
    double f;

    RVector *dphi = new RVector[mesh->nQ];
    for (i = 0; i < mesh->nQ; i++) dphi[i].New (nlen);

    bsol->SetActiveParams (pscaler->Unscale(x));
    raster->Map_SolToMesh (*bsol, *msol, true);
    msol->RegisterChange();
    for (i = 0; i < nofwavel; i++) {
	fws->Reset (*msol->swsol[i], 0);
	fws->CalcFields (*qvec, dphi);
	RVector proj_i(proj, i*mesh->nQM, mesh->nQM);
	proj_i = fws->ProjectAll (*mvec, dphi);
    }

    if (g) {
        MW_get_gradient (*raster, *fws, dphi, *qvec, *mvec, msol, r, *data,
			 *sd);
	pscaler->ScaleGradient (bsol->GetActiveParams(), r);
	r += reg->GetGradient(x);
	memcpy (g, r.data_buffer(), n*sizeof(double));
    }

    f = ObjectiveFunction::get_value (*data, proj, *sd);
    f += reg->GetValue (x);

    delete []dphi;

    return f;
}

// ==========================================================================
// Callback function for LBFGS library solver:
// iteration progress output (write images and echo objective function)

static int progress (void *instance, const lbfgsfloatval_t *x_lbfgs,
    const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step, int n, int k, int ls)
{
    int i;
    LBFGS_DATA *ldata = (LBFGS_DATA*)instance;
    MWsolution *msol = ldata->msol;
    Solution *bsol = ldata->bsol;
    const Raster *raster = ldata->raster;
    const Scaler *pscaler = ldata->pscaler;
    Regularisation *reg = ldata->reg;
    int blen = raster->BLen();

    RVector x(n, (double*)x_lbfgs);
    bsol->SetActiveParams (pscaler->Unscale(x));
    raster->Map_SolToMesh (*bsol, *msol, true);

    if (g_imgfmt != IMGFMT_RAW) {
	for (i = 0; i < msol->nParam(); i++) {
	    char fname[256];
	    if (msol->IsActive(i)) {
		if (i < msol->nmuaChromo) 
		    snprintf (fname, sizeof(fname), "%sreconChromophore_%d.nim",
		        g_prefix,i+1);
		else if (i == msol->nmuaChromo)
		    snprintf (fname, sizeof(fname), "%sreconScatPrefactor_A.nim",
			g_prefix);
		else if (i == msol->nmuaChromo + 1) 
		    snprintf (fname, sizeof(fname), "%sreconScatPower_b.nim",
			g_prefix);
		msol->WriteImgGeneric (k, fname, i);
	    }
	}
    }
    if (g_imgfmt != IMGFMT_NIM) {
        Solution gsol(msol->nParam(), blen);
	raster->Map_SolToBasis (*bsol, gsol, true);
	for (i = 0; i < msol->nParam(); i++) {
	    char fname[256];
	    if (msol->IsActive(i)) {
	        if (i < msol->nmuaChromo)
		    snprintf (fname, sizeof(fname),  "%sreconChromophore_%d.raw",
		        g_prefix,i+1);
		else if (i == msol->nmuaChromo)
		    snprintf (fname, sizeof(fname),  "%sreconScatPrefactor_A.raw",
		        g_prefix);
		else if (i == msol->nmuaChromo+1)
		    snprintf (fname, sizeof(fname), "%sreconScatPower_b.raw",
			g_prefix);
		else
		  xERROR("Invalid parameter index during output");
		gsol.WriteImgGeneric (k, fname, i);
	    }
	}
    }
    double pr = reg->GetValue (x);
    LOGOUT("Iteration %d  CPU %f  OF %g [LH %g PR %g]",
	   k, toc(clock0), fx, fx-pr, pr);
    return 0;
}

// ==========================================================================

void SolverLBFGS_CW_MW::Solve (RFwdSolverMW &FWS, const Raster &raster,
    const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
    const RVector &sd, Solution &bsol, MWsolution &msol,
    const RCompRowMatrix &qvec, const RCompRowMatrix &mvec)
{
    // Limited memory BFGS solver
    // Uses identity or diag(JTJ) as initial Hessian

    // Ref: Byrd et. al. "Representations of Quasi-Newton matrices and
    // their use in limited memory methods", Mathematical Programming 63(4),
    // pp. 129-156 (1996)

    LOGOUT ("SOLVER: L-BFGS");

    int n = bsol.ActiveDim();
    int ret;
    lbfgs_parameter_t param;

    lbfgs_parameter_init (&param);
    param.m = history;
    param.epsilon = epsilon;
    param.delta = delta;
    param.past = 1;
    param.max_iterations = itmax;

    RVector x0 = pscaler->Scale (bsol.GetActiveParams());
    Regularisation *reg = Regularisation::Create (pp, &x0, &raster);

    LBFGS_DATA ldata;
    ldata.of = &OF;
    ldata.raster = &raster;
    ldata.fws = &FWS;
    ldata.reg = reg;
    ldata.qvec = &qvec;
    ldata.mvec = &mvec;
    ldata.pscaler = pscaler;
    ldata.msol = &msol;
    ldata.bsol = &bsol;
    ldata.data = &data;
    ldata.sd = &sd;

    // find initial guess for minimizer
    RVector x = pscaler->Scale (bsol.GetActiveParams());

    lbfgsfloatval_t *x_lbfgs = (lbfgsfloatval_t*)x.data_buffer();
    lbfgsfloatval_t fx;

    // initial objective function
    double of = (double)evaluate (&ldata, x_lbfgs, 0, x.Dim(), 0);
    double fp = reg->GetValue (x);
    LOGOUT("Iteration 0  CPU %f  OF %g [LH %g PR %g]", toc(clock0),
	   of, of-fp, fp);

    ret = lbfgs (n, x_lbfgs, &fx, evaluate, progress, &ldata, &param);

    delete reg;
}

// ==========================================================================

void SolverLBFGS_CW_MW::ReadParams (ParamParser &pp)
{
    // === CONVERGENCE LIMIT ===
    if (!(pp.GetReal ("LBFGS_TOL", epsilon) ||
	  pp.GetReal ("NONLIN_TOL", epsilon)) || epsilon <= 0) {
        cout << "\nL-BFGS nonlinear solver --------------------------------\n";
        cout << "Enter convergence criterion epsilon:\n";
	cout << "|g(x)| / max(1,|x|) < epsilon\n";
	cout << "with parameter vector x and gradient g.\n";
	    do {
		cout << "\nNONLIN_TOL (float, >0):\n>> ";
		cin >> epsilon;
	    } while (epsilon <= 0);
    }

    // === STOPPING CRITERION ===
    if (!pp.GetReal ("NONLIN_DELTA", delta)) {
        cout << "\nL-BFGS nonlinear solver --------------------------------\n";
        cout << "Enter stopping criterion delta. This stops the\n";
	cout << "reconstruction if [f(x_(k-1))-f(x)]/f(x) < delta\n";
        do {
	    cout << "\nNONLIN_DELTA (float, >=0):\n>> ";
	    cin >> delta;
	} while (delta < 0.0);
    }

    // === MAX ITERATION COUNT ===
    if (!(pp.GetInt ("LBFGS_ITMAX", itmax) ||
	  pp.GetInt ("NONLIN_ITMAX", itmax)) || itmax <= 0) {
        cout << "\nL-BFGS nonlinear solver --------------------------------\n";
        cout << "Enter the maximum number of iterations to be performed if\n";
	cout << "the convergence and tolerance limits are not satisfied:\n";
	do {
	    cout << "\nNONLIN_ITMAX (int, >0):\n";
	    cout << ">> ";
	    cin >> itmax;
	} while (itmax <= 0);
    }

    // === NUMBER OF PREVIOUS VECTORS TO STORE ===
    if (!pp.GetInt ("LBFGS_HISTORY", history) || history <= 0) {
        cout << "\nL-BFGS nonlinear solver --------------------------------\n";
        cout << "Enter the number of corrections to approximate the inverse\n";
	cout << "Hessian matrix. Larger values require more memory, but can\n";
	cout << "improve convergence rates. Values less than 3 are not\n";
	cout << "recommended.\n";
	do {
	    cout << "\nLBFGS_HISTORY (int, >=3):\n";
	    cout << ">> ";
	    cin >> history;
	} while (history <= 0);
    }
}

// ==========================================================================

void SolverLBFGS_CW_MW::WriteParams (ParamParser &pp)
{
    pp.PutString ("SOLVER", "LBFGS");
    pp.PutInt ("NONLIN_ITMAX", itmax);
    pp.PutReal ("NONLIN_TOL", epsilon);
    pp.PutReal ("NONLIN_DELTA", delta);
    pp.PutInt ("LBFGS_HISTORY", history);
}


// ==========================================================================

static RVector single_gradient_data (const Raster &raster,
    RFwdSolverMW &FWS, const RVector &proj, RVector *dphi,
    const RCompRowMatrix &mvec, const RVector &data,
    const RVector &sd)
{
    const QMMesh &mesh = *FWS.meshptr;
    int i, j, q, m, n, idx, ofs_mod, ofs_arg, nQ = mesh.nQ;
    int glen = raster.GLen();
    int slen = raster.SLen();
    int dim  = raster.Dim();
    double term;
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    const int *elref = raster.Elref();
    RVector grad(slen*2);
    RVector wqa (mesh.nlen());
    RVector wqb (mesh.nlen());
    RVector dgrad (slen);
    ofs_mod = 0;         // data offset for Mod data
    ofs_arg = mesh.nQM;  // data offset for Arg data
    RVector grad_cmua(grad, 0, slen);       // mua part of grad
    RVector grad_ckappa (grad, slen, slen); // kappa part of grad
    
    double tm_mesh2grid = 0.0;
    double tm_grid2sol = 0.0;
    double tm_gradient = 0.0;
    double tm_innerloop = 0.0;

    LOGOUT1_INIT_PROGRESSBAR ("AddGradient", 50, mesh.nQ);
    for (q = 0; q < mesh.nQ; q++) {

        // expand field and gradient
        RVector cdfield (glen);
        RVector *cdfield_grad = new RVector[dim];
	tic();
	raster.Map_MeshToGrid (dphi[q], cdfield);
	tm_mesh2grid += toc();
	tic();
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);
	tm_gradient += toc();

        n = mesh.nQMref[q];

	RVector y_mod (data, ofs_mod, n);
	RVector s_mod (sd, ofs_mod, n);
	RVector ypm_mod (proj, ofs_mod, n);
	RVector b_mod(n);
	b_mod = (y_mod-ypm_mod)/s_mod;

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	RVector cproj(n);
	cproj = FWS.ProjectSingle (q, mvec, dphi[q], DATA_LIN);
	wqa = 0.0;
	wqb = 0.0;

	tic();
	for (m = idx = 0; m < mesh.nM; m++) {
	    if (!mesh.Connected (q, m)) continue;
	    const RVector qs = mvec.Row(m);
	    double rp = cproj[idx];
	    double dn = 1.0/(rp*rp);

	    // amplitude term
	    term = -/* 2.0 * */ b_mod[idx] / (ype[idx]*s_mod[idx]);
	    wqa += qs * (term*rp*dn);

	    //wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}
	tm_innerloop += toc();

	// adjoint field and gradient
	RVector wphia (mesh.nlen());
	FWS.CalcField (wqa, wphia);

	RVector cafield(glen);
	RVector *cafield_grad = new RVector[dim];
	tic();
	raster.Map_MeshToGrid (wphia, cafield);
	tm_mesh2grid += toc();
	tic();
	ImageGradient (gdim, gsize, cafield, cafield_grad, elref);
	tm_gradient += toc();

	// absorption contribution
	tic();
	raster.Map_GridToSol (cdfield * cafield, dgrad);

	tm_grid2sol += toc();
	grad_cmua -= dgrad;

	// diffusion contribution
	// multiply complex field gradients
	RVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	tic();
	raster.Map_GridToSol (gk, dgrad);

	tm_grid2sol += toc();
	grad_ckappa -= dgrad;

	ofs_mod += n; // step to next source
	ofs_arg += n;

	delete []cdfield_grad;
	delete []cafield_grad;
	LOGOUT1_PROGRESS (q);
    }
    return grad;
}

// ==========================================================================

void MW_get_gradient (const Raster &raster, RFwdSolverMW &FWS,
		      RVector *dphi,
		      const RCompRowMatrix &qvec, const RCompRowMatrix &mvec,
		      const MWsolution *msol, RVector &grad,
		      const RVector &data, const RVector &sd)
{
    int i, j;
    const QMMesh *mesh = FWS.MeshPtr();
    RDenseMatrix extcoef = msol->extcoef;
    int nqm     = mesh->nQM;
    int slen    = raster.SLen();
    int nchromo = extcoef.nCols();
    int nlambda = extcoef.nRows();
    int nprm    = msol->nParam();
    bool bFactorA = msol->IsActive(nchromo);
    bool bPowerb  = msol->IsActive(nchromo+1);
    bool bScatter = bFactorA || bPowerb;
    RVector proj_i(mesh->nQM);
    RVector data_i(mesh->nQM);
    RVector sd_i(mesh->nQM);
    RVector sA, sb;

    grad.Clear();
    for (i = 0; i < nlambda; i++) {
	FWS.Reset (*msol->swsol[i], 0);
	FWS.CalcFields (qvec, dphi);
	proj_i = FWS.ProjectAll (mvec, dphi);
	data_i = RVector(data,i*nqm,nqm);
	sd_i = RVector(sd,i*nqm,nqm);
	//RVector data_i(data,i*nqm,nqm);
	//RVector sd_i(sd,i*nqm,nqm);
	RVector sgrad = single_gradient_data (raster, FWS, proj_i, dphi, mvec, data_i, sd_i);
	RVector sgrad_mua(sgrad,0,slen);
	RVector sgrad_kap(sgrad,slen,slen);

	// chromophore contributions
	for (j = 0; j < nchromo; j++) {
	    RVector grad_ch(grad, j*slen, slen);
	    grad_ch += sgrad_mua*extcoef(i,j);
	}

	// scattering prefactor
	if (bFactorA) {
	    raster.Map_MeshToSol (msol->GetJacobianCoeff_A(i), sA);
	    RVector grad_A(grad, nchromo*slen, slen);
	    grad_A += sgrad_kap*sA;
	}

	// scattering power
	if (bPowerb) {
	    raster.Map_MeshToSol (msol->GetJacobianCoeff_b(i), sb);
	    RVector grad_b(grad, (nchromo+1)*slen, slen);
	    grad_b += sgrad_kap*sb;
	}
    }
}
