#define __FWDSOLVER_MW_CC
#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "fwdsolver_mw.h"
#include "timing.h"

using namespace std;

// =========================================================================

template<class T>
TFwdSolverMW<T>::TFwdSolverMW (const QMMesh *mesh, LSOLVER linsolver,
    double tol)
    : TFwdSolver<T> (mesh, linsolver, tol)
{
    Setup();
}

// =========================================================================

template<class T>
TFwdSolverMW<T>::TFwdSolverMW (const QMMesh *mesh, char *solver, double tol)
    : TFwdSolver<T> (mesh, solver, tol)
{
    Setup();
}

// =========================================================================

template<class T>
TFwdSolverMW<T>::TFwdSolverMW (const QMMesh *mesh, ParamParser &pp)
    : TFwdSolver<T> (mesh, pp)
{
    Setup();
}

// =========================================================================

template<class T>
TFwdSolverMW<T>::~TFwdSolverMW ()
{
    Cleanup();
}

// =========================================================================

template<class T>
TVector<T> TFwdSolverMW<T>::ProjectAll_wavel (const TCompRowMatrix<T> &qvec,
    const TCompRowMatrix<T> &mvec, const MWsolution &sol, double omega,
    DataScale scl)
{
    const QMMesh *mesh = this->meshptr;
    int i, n = mesh->nlen(), nq = mesh->nQ;
    int nofwavel = sol.nofwavel;
    int w;

    TVector<T> proj(nofwavel*mesh->nQM);

    TVector<T> *phi = new TVector<T>[nq];
    for (i = 0; i < nq; i++) phi[i].New(n);

    for (w = 0; w < nofwavel; w++) {
	this->Reset (*sol.swsol[w], omega);
	this->CalcFields (qvec, phi);
	TVector<T> proj_w (proj, w*mesh->nQM, mesh->nQM);
	proj_w = this->ProjectAll (mvec, phi, scl);
    }
    delete []phi;

    return proj;
}

// =========================================================================

template<>
RVector RFwdSolverMW::ProjectAll_wavel_real (const RCompRowMatrix &qvec,
    const RCompRowMatrix &mvec, const MWsolution &sol, double omega,
    DataScale scl)
{
    // real case: nothing to do
    return ProjectAll_wavel (qvec, mvec, sol, omega, scl);
}

template<>
RVector CFwdSolverMW::ProjectAll_wavel_real (const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, const MWsolution &sol, double omega,
    DataScale scl)
{
    // complex case: split real and imaginary parts
    CVector cproj = ProjectAll_wavel (qvec, mvec, sol, omega, scl);

    const QMMesh *mesh = this->meshptr;
    int nofwavel = sol.nofwavel;
    int nqm = mesh->nQM;
    int n = cproj.Dim();
    int i;
    RVector proj(n*2);
    for (i = 0; i < nofwavel; i++) {
        CVector cproj_i(cproj, i*nqm, nqm);
	RVector proj_re (proj, i*2*nqm, nqm);
	RVector proj_im (proj, (i*2+1)*nqm, nqm);
	proj_re = Re(cproj_i);
	proj_im = Im(cproj_i);
    }
    return proj;
}

// =========================================================================

template<class T>
void TFwdSolverMW<T>::Setup ()
{
}

// =========================================================================

template<class T>
void TFwdSolverMW<T>::Cleanup ()
{
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class STOASTLIB TFwdSolverMW<double>;
template class STOASTLIB TFwdSolverMW<std::complex<double> >;

#endif // NEED_EXPLICIT_INSTANTIATION
