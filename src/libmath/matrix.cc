// ==========================================================================
// Module mathlib
// File matrix.cc
// Definition of template class TMatrix
// ==========================================================================

#define MATHLIB_IMPLEMENTATION

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

#include "mathlib.h"

using namespace std;
using namespace toast;

// ==========================================================================
// typedefs for specific instances of `TMatrix'
// These are local because TMatrix should not be used publicly anyway
// ==========================================================================
typedef TMatrix<double>  RMatrix;	// 'real'
typedef TMatrix<toast::complex> CMatrix;	// 'toast::complex'
typedef TMatrix<int>     IMatrix;	// 'integer'
#ifdef TOAST_FEATURE_SINGLEPREC
typedef TMatrix<float>   FMatrix;	// 'float'
#endif
// ==========================================================================

// ==========================================================================
// member definitions

template<class MT>
MATHLIB TSymMatrix<MT> AAT (const TMatrix<MT> &A)
{
    int i, j, nr = A.nRows(), nc = A.nCols();
    TSymMatrix<MT> aat(nr);
    for (i = 0; i < nr; i++) {
        TVector<MT> row = A.Row(i);
	for (j = 0; j <= i; j++)
	    aat(i,j) = row & A.Row(j);
    }
    return aat;
}

template<class MT>
TVector<MT> ATA_diag (const TMatrix<MT> &A)
{
    int i, j;
    int n = A.nRows();
    int m = A.nCols();
    MT sum, Ai;
    TVector<MT> diag(m);

    for (j = 0; j < m; j++) {
        sum = (MT)0;
	for (i = 0; i < n; i++) {
	    Ai = A(i,j);
	    sum += Ai*Ai;
	}
	diag[j] = sum;
    }
    return diag;
}


// ==========================================================================
// friend definitions

template<class MT>
MATHLIB ostream &operator<< (ostream &os, const TMatrix<MT> &mat)
{
    os << '[';
    for (int i = 0; i < mat.rows; i++) {
	if (i) os << '\n';
	//os << mat.Row(i);

	os << '[';
	for (int j = 0; j < mat.cols; j++) {
	    if (j) os << ' ';
	    os << mat.Get (i, j);
	}
	os << ']';
    }
    os << ']';
    return os;
}


// ==========================================================================
// Preconditioned conjugate gradients (PCG)
// ==========================================================================
// ==========================================================================

template<class MT>
MATHLIB int PCG (const TMatrix<MT> &A, const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon, int maxit)
{
    dASSERT(A.rows == A.cols, "Matrix not square");
    dASSERT(b.Dim() == A.rows, "Dimension mismatch");
    dASSERT(x.Dim() == A.rows, "Dimension mismatch");

    double dnew, dold, alpha, beta;
    double err, bnorm;
    int niter, dim = x.Dim();
    if (!maxit) maxit = dim+1;
    TVector<MT> r(dim), d(dim), q(dim);
    A.Ax (x, r);
    r = b - r;
    if (precon) precon->Apply (r, d);
    else d = r;
    dnew = r & d;
    bnorm = l2norm (b);

    for (niter = 0; niter < maxit; niter++) {
        A.Ax (d, q);
	alpha = dnew / (d & q);
	
	// check failure
	if (!alpha || std::isnan(alpha)) {
	    cerr << "*** PCG fails to converge" << endl;
	    break;
	}

	x += d * alpha;
	r -= q * alpha;
	dold = dnew;
	if (precon) {
	    precon->Apply (r, q);
	    dnew = r & q;
	    beta = dnew/dold;
	    d *= beta;
	    d += q;
	} else {
	    dnew = r & r;
	    beta = dnew/dold;
	    d *= beta;
	    d += r;
	}

	// check failure
	if (!beta || std::isnan(beta)) {
	    cerr << "*** PCG fails to converge" << endl;
	    break;
	}

	// check convergence
	err = l2norm (r);
	//cerr << niter << " " << err/bnorm << endl;
	if (err <= tol*bnorm) break;
    }
    tol = err/bnorm;
    return niter;
}

template<class MT>
MATHLIB void PCG (const TMatrix<MT> &A, const TVector<MT> *b, TVector<MT> *x,
    int nrhs, double tol, int maxit, TPreconditioner<MT> *precon,
    IterativeSolverResult *res)
{
    int i, it;
    double tol0 = tol;    
    if (res) res->it_count = 0;

    for (i = 0; i < nrhs; i++) {
        tol = tol0;
	it = PCG (A, b[i], x[i], tol, precon, maxit);
	if (res) {
	    res->it_count += it;
	    if (!i || tol > res->rel_err) res->rel_err = tol;
	}
    }
}

// ==========================================================================
// Preconditioned conjugate gradients (PCG) with callback parameter
// ==========================================================================

template<class MT>
MATHLIB int PCG (TVector<MT> (*Mv_clbk)( const TVector<MT> &v,
    void *context), void * context, const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon, int maxit)
{
    double dnew, dold, alpha, beta;
    double err, bnorm;
    int niter, dim = x.Dim();
    if (!maxit) maxit = dim+1;
    TVector<MT> r(dim), d(dim), q(dim);
    r = b - ((*Mv_clbk)(x, context));
    if (precon) precon->Apply (r, d);
    else d = r;
    dnew = r & d;
    bnorm = l2norm (b);

    for (niter = 0; niter < maxit; niter++) {
        q = (*Mv_clbk)(d, context);
	//A.Ax (d, q);
	alpha = dnew / (d & q);
	x += d * alpha;
	r -= q * alpha;
	dold = dnew;
	if (precon) {
	    precon->Apply (r, q);
	    dnew = r & q;
	    beta = dnew/dold;
	    d *= beta;
	    d += q;
	} else {
	    dnew = r & r;
	    beta = dnew/dold;
	    d *= beta;
	    d += r;
	}
	// check convergence
	err = l2norm (r);
	cout << "PCG: iteration " << niter << " err = " << err/bnorm << endl;
	cerr << niter << " " << err/bnorm << endl;
	if (err <= tol*bnorm) break;
    }
    tol = err/bnorm;
    return niter;
}

// ==========================================================================
// Preconditioned bi-conjugate gradients (BiPCG)
// from: netlib (http://www.netlib.org/linalg/html_templates/node32.html)
// ==========================================================================

template<class MT>
int BiPCG (const TMatrix<MT> &A, const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon, int maxit)
{
    int i, niter, dim = x.Dim();
    double rho, alpha, aden, beta, bden, err, bnorm;
    TVector<MT> r(dim), rd(dim), z(dim), zd(dim);
    TVector<MT> p(dim), pd(dim), q(dim), qd(dim);

    if (!maxit) maxit = dim+1;
    r = rd = b - (A*x);
    bnorm = l2norm(b);

    for (niter = 0; niter < maxit; niter++) {
        if (precon) {
	    precon->Apply (r, z);
	    precon->Apply (rd, zd);
	    // note this assumes symmetric preconditioner!
	    // generally we need to solve M^T zd = rd, rather than M zd = rd
	} else {
	    z = r;
	    zd = rd;
	}

	bden = rho;
	rho = z & rd;
	xASSERT (rho != 0.0, "BiCG solver fails");
	if (!niter) {
	    p = z;
	    pd = zd;
	} else {
	    beta = rho/bden;
	    for (i = 0; i < dim; i++) {
	        p[i]  *= beta; p[i]  += z[i];
		pd[i] *= beta; pd[i] += zd[i];
	    }
	}
	A.Ax (p, q);
	A.ATx (pd, qd);
	aden = pd & q;
	alpha = rho/aden;
	for (i = 0; i < dim; i++) {
	    x[i]  += alpha * p[i];
	    r[i]  -= alpha * q[i];
	    rd[i] -= alpha * qd[i];
	}

	// check convergence
	err = l2norm(r);
	//cerr << "BiCG tol=" << err/bnorm << endl;
	if (err <= tol*bnorm) break;
    }
    tol = err/bnorm;
    return niter;
}

template<> // specialisation: complex
// should try to merge this into the templated version
int BiPCG<toast::complex> (const CMatrix &A, const CVector &b, CVector &x,
    double &tol, CPreconditioner *cprecon, int maxit)
{
    int i, niter, dim = x.Dim();
    double rho=0, alpha, aden, beta, bden, err, bnorm;
    if (!maxit) maxit = dim+1;
    CVector r(dim), rd(dim), z(dim), zd(dim), p(dim), pd(dim), q(dim), qd(dim);
    r = rd = b - (A*x);
    bnorm = l2norm(b);

    for (niter = 0; niter < maxit; niter++) {
        if (cprecon) {
	    cprecon->Apply (r, z);
	    cprecon->Apply (rd, zd);
	    // note this assumes symmetric preconditioner!
	    // generally we need to solve M^T zd = rd, rather than M zd = rd
	} else {
	    z = r;
	    zd = rd;
	}
	bden = rho;
	for (rho = 0.0, i = 0; i < dim; i++)
	    rho += z[i].re*rd[i].re + z[i].im*rd[i].im;
	    // is this equivalent to rho = re(conj(z) & rd) ?
	xASSERT(rho != 0.0, "BiCG solver fails");
	if (!niter) {
	    p = z;
	    pd = zd;
	} else {
	    beta = rho/bden;
	    for (i = 0; i < dim; i++) {
	        p[i].re  *= beta;  p[i].re += z[i].re;
		p[i].im  *= beta;  p[i].im += z[i].im;
		pd[i].re *= beta;  pd[i].re += zd[i].re;
		pd[i].im *= beta;  pd[i].im += zd[i].im;
	    }
	}
	A.Ax (p, q);
	A.ATx (pd, qd);
	for (aden = 0.0, i = 0; i < dim; i++)
	    aden += pd[i].re*q[i].re + pd[i].im*q[i].im;
	alpha = rho/aden;
	for (i = 0; i < dim; i++) {
	    x[i].re  += alpha*p[i].re;
	    x[i].im  += alpha*p[i].im;
	    r[i].re  -= alpha*q[i].re;
	    r[i].im  -= alpha*q[i].im;
	    rd[i].re -= alpha*qd[i].re;
	    rd[i].im -= alpha*qd[i].im;
	}

	// check convergence
	for (err = 0.0, i = 0; i < dim; i++)
	    err += r[i].re*r[i].re + r[i].im*r[i].im;
	err = sqrt (err);
	//cerr << "BiCG tol=" << err/bnorm << endl;
	if (err <= tol*bnorm) break;
    }
    tol = err/bnorm;
    return niter;
}

// ==========================================================================
// Preconditioned bi-conjugate gradient Stabilised (BiCGSTAB) solver
// from: netlib (http://www.netlib.org/linalg/html_templates/node41.html)
// ==========================================================================

template<class MT>
MATHLIB int BiCGSTAB (const TMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon, int maxit)
{
    int i, niter, dim = x.Dim();
    double rho=0, alpha, aden, beta, bden, omega, err, bnorm;
    TVector<MT> r(dim), rd(dim), p(dim), pd(dim), s(dim), sd(dim), t(dim), v(dim);
    if (!maxit) maxit = dim+1;
    r = rd = b - Ax(A,x);
    bnorm = l2norm (b);

    for (niter = 0; niter < maxit; niter++) {
        bden = rho;
        rho = rd & r;
	if (!rho) {
	    cerr << "BiCGSTAB solver fails to converge" << endl;
	    break;
	}
	if (!niter) {
	    p = r;
	} else {
	    beta = (rho*alpha)/(bden*omega);
	    for (i = 0; i < dim; i++)
	        p[i] = (MT)(r[i] + beta * (p[i] - omega*v[i]));
	}
	if (precon) precon->Apply (p, pd);
	else        pd = p;
	v = Ax(A,pd);
	aden = rd & v;
	alpha = rho/aden;
	for (i = 0; i < dim; i++)
	    s[i] = (MT)(r[i] - alpha*v[i]);

	// check break condition 1
	err = l2norm(s);
	if (err < tol*bnorm) {
	    for (i = 0; i < dim; i++) x[i] += (MT)(alpha*pd[i]);
	    break;
	}

	if (precon) precon->Apply (s, sd);
	else        sd = s;
	t = Ax(A,sd);
	omega = (t & s) / (t & t);
	for (i = 0; i < dim; i++) {
	    x[i] += (MT)(alpha * pd[i] + omega * sd[i]);
	    r[i] = (MT)(s[i] - omega * t[i]);
	}

	// check convergence
	err = l2norm (r);
	if (err <= tol*bnorm) break;
	if (omega == 0.0) {
	    cerr << "BiCGSTAB solver fails to converge" << endl;
	    break;
	}
	//cerr << "BICGSTAB: err = " << err/bnorm << endl;
    }
    if (toastVerbosity > 1) {
        if (err > tol*bnorm)
	    cout << "BiCGSTAB residual " << err/bnorm << " after " << niter
		 << " iterations" << endl;
	else
	    cout << "BiCGSTAB converged after " << niter
		 << " iterations (Residual " << err/bnorm << ")" << endl;
    }
    tol = err/bnorm;
    return niter;
}

template<> // specialisation: complex
MATHLIB int BiCGSTAB<toast::complex> (const CMatrix &A, const CVector &b,
    CVector &x, double &tol, CPreconditioner *precon, int maxit)
{
    int i, niter, dim = x.Dim();
    double err, bnorm, rho=0, rhop, alpha, beta, omega, aden, onum, oden;
    if (!maxit) maxit = dim+1;
    CVector r(dim), rd(dim), p(dim), pd(dim), v(dim), s(dim), sd(dim), t(dim);
    r = rd = b - Ax(A,x);
    bnorm = l2norm (b);

    for (niter = 0; niter < maxit; niter++) {
        rhop = rho;
        for (rho = 0.0, i = 0; i < dim; i++)
	    rho += rd[i].re*r[i].re + rd[i].im*r[i].im;
	xASSERT(rho != 0.0, "Bi-CGSTAB fails");
	if (!niter) {
	    p = r;
	} else {
	    beta = (rho*alpha)/(rhop*omega);
	    for (i = 0; i < dim; i++) {
	        p[i].re = r[i].re + beta * (p[i].re-omega*v[i].re);
		p[i].im = r[i].im + beta * (p[i].im-omega*v[i].im);
	    }
	}

	if (precon) precon->Apply (p, pd);
	else        pd = p;

	v = Ax(A,pd);
	for (aden = 0.0, i = 0; i < dim; i++)
	    aden += rd[i].re*v[i].re + rd[i].im*v[i].im;
	alpha = rho/aden;
	for (i = 0; i < dim; i++) {
	    s[i].re = r[i].re - alpha*v[i].re;
	    s[i].im = r[i].im - alpha*v[i].im;
	}
	if (precon) precon->Apply (s, sd);
	else        sd = s;
	t = Ax(A,sd);
	for (onum = oden = 0.0, i = 0; i < dim; i++) {
	    onum += t[i].re*s[i].re + t[i].im*s[i].im;
	    oden += t[i].re*t[i].re + t[i].im*t[i].im;
	}
	omega = onum/oden;
	for (i = 0; i < dim; i++) {
	    x[i].re += alpha*pd[i].re + omega*sd[i].re;
	    x[i].im += alpha*pd[i].im + omega*sd[i].im;
	    r[i].re = s[i].re - omega*t[i].re;
	    r[i].im = s[i].im - omega*t[i].im;
	}

	// check convergence
	err = l2norm (r);
	if (err <= tol*bnorm) break;
	if (omega == 0.0) break; // no convergence
	//cerr << "BICGSTAB: err = " << err/bnorm << endl;
    }
    if (toastVerbosity > 1) {
        if (err > tol*bnorm)
	    cout << "BiCGSTAB residual " << err/bnorm << " after " << niter
		 << " iterations" << endl;
	else
	    cout << "BiCGSTAB converged after " << niter
		 << " iterations (Residual " << err/bnorm << ")" << endl;
    }
    tol = err/bnorm;
    return niter;
}

#ifdef TOAST_FEATURE_SINGLEPREC
template<> // specialisation: single complex
MATHLIB int BiCGSTAB<scomplex> (const SCMatrix &A, const SCVector &b,
    SCVector &x, double &tol, SCPreconditioner *precon, int maxit)
{
    int i, niter, dim = x.Dim();
    float err, bnorm, rho=0, rhop, alpha, beta, omega, aden, onum, oden;
    if (!maxit) maxit = dim+1;
    SCVector r(dim), rd(dim), p(dim), pd(dim), v(dim), s(dim), sd(dim), t(dim);
    r = rd = b - Ax(A,x);
    bnorm = l2norm (b);

    for (niter = 0; niter < maxit; niter++) {
        rhop = rho;
        for (rho = 0.0, i = 0; i < dim; i++)
	    rho += rd[i].re*r[i].re + rd[i].im*r[i].im;
	xASSERT(rho != 0.0, "Bi-CGSTAB fails");
	if (!niter) {
	    p = r;
	} else {
	    beta = (rho*alpha)/(rhop*omega);
	    for (i = 0; i < dim; i++) {
	        p[i].re = r[i].re + beta * (p[i].re-omega*v[i].re);
		p[i].im = r[i].im + beta * (p[i].im-omega*v[i].im);
	    }
	}
	if (precon) precon->Apply (p, pd);
	else        pd = p;

	v = Ax(A,pd);
	for (aden = 0.0, i = 0; i < dim; i++)
	    aden += rd[i].re*v[i].re + rd[i].im*v[i].im;
	alpha = rho/aden;
	for (i = 0; i < dim; i++) {
	    s[i].re = r[i].re - alpha*v[i].re;
	    s[i].im = r[i].im - alpha*v[i].im;
	}
	if (precon) precon->Apply (s, sd);
	else        sd = s;
	t = Ax(A,sd);
	for (onum = oden = 0.0, i = 0; i < dim; i++) {
	    onum += t[i].re*s[i].re + t[i].im*s[i].im;
	    oden += t[i].re*t[i].re + t[i].im*t[i].im;
	}
	omega = onum/oden;
	for (i = 0; i < dim; i++) {
	    x[i].re += alpha*pd[i].re + omega*sd[i].re;
	    x[i].im += alpha*pd[i].im + omega*sd[i].im;
	    r[i].re = s[i].re - omega*t[i].re;
	    r[i].im = s[i].im - omega*t[i].im;
	}

	// check convergence
	err = l2norm (r);
	if (err <= tol*bnorm) break;
	if (omega == 0.0) break; // no convergence
	//cerr << "BICGSTAB: err = " << err/bnorm << endl;
    }
    if (toastVerbosity > 1) {
        if (err > tol*bnorm)
	    cout << "BiCGSTAB residual " << err/bnorm << " after " << niter
		 << " iterations" << endl;
	else
	    cout << "BiCGSTAB converged after " << niter
		 << " iterations (Residual " << err/bnorm << ")" << endl;
    }
    tol = err/bnorm;
    return niter;
}
#endif

// ==========================================================================

template<class MT>
MATHLIB void BiCGSTAB (const TMatrix<MT> &A, const TVector<MT> *b,
    TVector<MT> *x,  int nrhs, double tol, int maxit,
    TPreconditioner<MT> *precon, IterativeSolverResult *res)
{
    int i, it;
    double tol0 = tol;    
    if (res) res->it_count = 0;

    for (i = 0; i < nrhs; i++) {
        tol = tol0;
	it = BiCGSTAB (A, b[i], x[i], tol, precon, maxit);
	if (res) {
	    res->it_count += it;
	    if (!i || tol > res->rel_err) res->rel_err = tol;
	}
    }
}

// ==========================================================================

// In this case MVM Matrix vector multiply returns a vector.
template<class MT>
MATHLIB int BiCGSTAB (TVector<MT> (*Mv_clbk)( const TVector<MT> &v,
    void *context), void *context, const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon, int maxit)
{
    int i, niter, dim = x.Dim();
    MT rho = (MT)0, aden, bden, alpha, beta, omega;
    double err, bnorm;
    TVector<MT> r(dim), rd(dim), p(dim), pd(dim), s(dim), sd(dim), t(dim), v(dim);
    if (!maxit) maxit = dim+1;
    //    r = rd = b - (A*x);
    r =  (*Mv_clbk)(x, context);  // replace r <- A*x;
    r = rd = b - r;
    bnorm = l2norm (b);

    for (niter = 0; niter < maxit; niter++) {
        bden = rho;
        rho = rd & r;
	if (!rho) {
	    cerr << "BiCGSTAB solver fails to converge" << endl;
	    break;
	}
	if (!niter) {
	    p = r;
	} else {
	    beta = (rho*alpha)/(bden*omega);
	    for (i = 0; i < dim; i++)
	        p[i] = r[i] + beta * (p[i] - omega*v[i]);
	}
	if (precon) precon->Apply (p, pd);
	else        pd = p;
	//	A.Ax (pd, v);
	v = (*Mv_clbk) (pd, context); //  multiply
	aden = rd & v;
	alpha = rho/aden;
	for (i = 0; i < dim; i++)
	    s[i] = r[i] - alpha*v[i];

	// check break condition 1
	err = l2norm(s);
	if (err < tol*bnorm) {
	    for (i = 0; i < dim; i++) x[i] += alpha*pd[i];
	    break;
	}

	if (precon) precon->Apply (s, sd);
	else        sd = s;
	//	A.Ax (sd, t);
	t = (*Mv_clbk)(sd, context);
	omega = (t & s) / (t & t);
	for (i = 0; i < dim; i++) {
	    x[i] += alpha * pd[i] + omega * sd[i];
	    r[i] = s[i] - omega * t[i];
	}

	// check convergence
	err = l2norm (r);
	if (err <= tol*bnorm) break;
	if (omega == 0.0) {
	    cerr << "BiCGSTAB solver fails to converge" << endl;
	    break;
	}
	//cerr << "BICGSTAB: err = " << err/bnorm << endl;
    }
    //cerr << "BICGSTAB: " << niter << " iterations" << endl;
    if (err > tol*bnorm)
        cerr << "BiCGSTAB residual " << err/bnorm << " after " << niter
	     << " iterations" << endl;
    tol = err/bnorm;
    return niter;
}

// This version uses pointer to function
// MVM is "in place" Matrix vector multiply.
template<class MT>
MATHLIB int BiCGSTAB (void (* MVM)( TVector<MT> &),  const TVector<MT> &b, 
    TVector<MT> &x, double &tol, TPreconditioner<MT> *precon, int maxit)
{
    int i, niter, dim = x.Dim();
    double rho=0, alpha, aden, beta, bden, omega, err, bnorm;
    TVector<MT> r(dim), rd(dim), p(dim), pd(dim), s(dim), sd(dim), t(dim), v(dim);
    if (!maxit) maxit = dim+1;
    //    r = rd = b - (A*x);
    r = x;       // first copy x.
    (* MVM)(r);  // replace r <- A*x;
    r = rd = b - r;
    bnorm = l2norm (b);

    for (niter = 0; niter < maxit; niter++) {
        bden = rho;
        rho = rd & r;
	if (!rho) {
	    cerr << "BiCGSTAB solver fails to converge" << endl;
	    break;
	}
	if (!niter) {
	    p = r;
	} else {
	    beta = (rho*alpha)/(bden*omega);
	    for (i = 0; i < dim; i++)
	        p[i] = r[i] + beta * (p[i] - omega*v[i]);
	}
	if (precon) precon->Apply (p, pd);
	else        pd = p;
	//	A.Ax (pd, v);
	v = pd; (*MVM) (v); // copy, and multiply
	aden = rd & v;
	alpha = rho/aden;
	for (i = 0; i < dim; i++)
	    s[i] = r[i] - alpha*v[i];

	// check break condition 1
	err = l2norm(s);
	if (err < tol*bnorm) {
	    for (i = 0; i < dim; i++) x[i] += alpha*pd[i];
	    break;
	}

	if (precon) precon->Apply (s, sd);
	else        sd = s;
	//	A.Ax (sd, t);
	t = sd; (*MVM)(t);
	omega = (t & s) / (t & t);
	for (i = 0; i < dim; i++) {
	    x[i] += alpha * pd[i] + omega * sd[i];
	    r[i] = s[i] - omega * t[i];
	}

	// check convergence
	err = l2norm (r);
	if (err <= tol*bnorm) break;
	if (omega == 0.0) {
	    cerr << "BiCGSTAB solver fails to converge" << endl;
	    break;
	}
	//cerr << "BICGSTAB: err = " << err/bnorm << endl;
    }
    if (toastVerbosity > 1) {
        if (err > tol*bnorm)
	    cout << "BiCGSTAB residual " << err/bnorm << " after " << niter
		 << " iterations" << endl;
	else
	    cout << "BiCGSTAB converged after " << niter
		 << " iterations (Residual " << err/bnorm << ")" << endl;
    }
    tol = err/bnorm;
    return niter;
}


// ==========================================================================
// Generalised minimal residuals (GMRES)
// ==========================================================================

template<class MT>
MATHLIB int GMRES (const TMatrix<MT> &A, const TVector<MT> &b, TVector<MT> &x,
    double &tol, TPreconditioner<MT> *precon, int restart, int maxit,
    void (*clbk)(void*))
{
    return gmres (restart, A, b, x, precon, tol, maxit, clbk);
}

template<class MT>
MATHLIB int GMRES (TVector<MT> (*Av_clbk)(const TVector<MT> &v, void *context),
    void *context, const TVector<MT> &b, TVector<MT> &x, double tol,
    TPreconditioner<MT> *precon, int restart, int maxit, int *iter, double *res)
{
    int it = gmres (restart, Av_clbk, b, x, precon, tol, maxit, context);
    if (iter) *iter = it;
    if (res) *res = tol;
    return 0;
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class MATHLIB TMatrix<double>;
template class MATHLIB TMatrix<toast::complex>;
template class MATHLIB TMatrix<int>;
#ifdef TOAST_FEATURE_SINGLEPREC
template class MATHLIB TMatrix<float>;
template class MATHLIB TMatrix<scomplex>;
#endif

template MATHLIB ostream &operator<< (ostream &os, const RMatrix &mat);
template MATHLIB ostream &operator<< (ostream &os, const CMatrix &mat);
template MATHLIB ostream &operator<< (ostream &os, const SCMatrix &mat);

template RVector  Ax (const RMatrix &A, const RVector &x);
template CVector  Ax (const CMatrix &A, const CVector &x);
template IVector  Ax (const IMatrix &A, const IVector &x);
#ifdef TOAST_FEATURE_SINGLEPREC
template FVector  Ax (const FMatrix &A, const FVector &x);
template SCVector Ax (const SCMatrix &A, const SCVector &x);
#endif

template RVector  ATx (const RMatrix &A, const RVector &x);
template CVector  ATx (const CMatrix &A, const CVector &x);
template IVector  ATx (const IMatrix &A, const IVector &x);
#ifdef TOAST_FEATURE_SINGLEPREC
template FVector  ATx (const FMatrix &A, const FVector &x);
template SCVector ATx (const SCMatrix &A, const SCVector &x);
#endif
template MATHLIB RSymMatrix  ATA (const RMatrix &A);
template MATHLIB FSymMatrix  ATA (const FMatrix &A);
template MATHLIB CSymMatrix  ATA (const CMatrix &A);
template MATHLIB SCSymMatrix ATA (const SCMatrix &A);

template MATHLIB RSymMatrix  AAT (const RMatrix &A);
template MATHLIB FSymMatrix  AAT (const FMatrix &A);
template MATHLIB CSymMatrix  AAT (const CMatrix &A);
template MATHLIB SCSymMatrix AAT (const SCMatrix &A);

template MATHLIB RVector  ATA_diag (const RMatrix &A);
template MATHLIB CVector  ATA_diag (const CMatrix &A);
#ifdef TOAST_FEATURE_SINGLEPREC
template MATHLIB FVector  ATA_diag (const FMatrix &A);
template MATHLIB SCVector ATA_diag (const SCMatrix &A);
#endif

template MATHLIB int PCG (const RMatrix &A, const RVector &b, RVector &x, double &tol,
     RPreconditioner *precon, int maxit);
template MATHLIB int PCG (RVector (*Mv_clbk)( const RVector &v,
    void *context), void * context, const RVector &b, RVector &x,
    double &tol, RPreconditioner *precon, int maxit);
template MATHLIB void PCG (const RMatrix &A, const RVector *b, RVector *x,
    int nrhs, double tol, int maxit, RPreconditioner *precon, IterativeSolverResult *res);
#ifdef TOAST_FEATURE_SINGLEPREC
template MATHLIB int PCG (const FMatrix &A, const FVector &b, FVector &x, double &tol,
     FPreconditioner *precon, int maxit);
template MATHLIB void PCG (const FMatrix &A, const FVector *b, FVector *x,
    int nrhs, double tol, int maxit, FPreconditioner *precon, IterativeSolverResult *res);
#endif

template MATHLIB int BiCGSTAB (const RMatrix &A, const RVector &b,
    RVector &x, double &tol, RPreconditioner *precon, int maxit);
template MATHLIB int BiCGSTAB (const IMatrix &A, const IVector &b,
    IVector &x, double &tol, IPreconditioner *precon, int maxit);
#ifdef TOAST_FEATURE_SINGLEPREC
template MATHLIB int BiCGSTAB (const FMatrix &A, const FVector &b,
    FVector &x, double &tol, FPreconditioner *precon, int maxit);
#endif

template MATHLIB void BiCGSTAB (const RMatrix &A, const RVector *b, RVector *x,
    int nrhs, double tol, int maxit, RPreconditioner *precon,
    IterativeSolverResult *res);
template MATHLIB void BiCGSTAB (const CMatrix &A, const CVector *b, CVector *x,
    int nrhs, double tol, int maxit, CPreconditioner *precon,
    IterativeSolverResult *res);
template MATHLIB void BiCGSTAB (const IMatrix &A, const IVector *b, IVector *x,
    int nrhs, double tol, int maxit, IPreconditioner *precon,
    IterativeSolverResult *res);
#ifdef TOAST_FEATURE_SINGLEPREC
template MATHLIB void BiCGSTAB (const FMatrix &A, const FVector *b, FVector *x,
    int nrhs, double tol, int maxit, FPreconditioner *precon,
    IterativeSolverResult *res);
template MATHLIB void BiCGSTAB (const SCMatrix &A, const SCVector *b,
    SCVector *x, int nrhs, double tol, int maxit,
    SCPreconditioner *precon, IterativeSolverResult *res);
#endif

template MATHLIB int BiCGSTAB (RVector(*Mv_clbk)(const RVector &v,
    void *context), void *context, const RVector &b, RVector &x,
    double &tol, RPreconditioner *precon, int maxit);
template MATHLIB int BiCGSTAB (CVector(*Mv_clbk)(const CVector &v,
    void *context), void *context, const CVector &b, CVector &x,
    double &tol, CPreconditioner *precon, int maxit);
template MATHLIB int BiCGSTAB (void (* MVM)(RVector &),  const RVector &b, 
    RVector &x, double &tol, RPreconditioner *precon, int maxit);

template MATHLIB int GMRES (const RMatrix &A, const RVector &b, RVector &x,
    double &tol, RPreconditioner *precon, int restart, int maxit,
    void(*clbk)(void*));
template MATHLIB int GMRES (const CMatrix &A, const CVector &b, CVector &x,
    double &tol, CPreconditioner *precon, int restart, int maxit,
    void(*clbk)(void*));
#ifdef TOAST_FEATURE_SINGLEPREC
template MATHLIB int GMRES (const FMatrix &A, const FVector &b, FVector &x,
    double &tol, FPreconditioner *precon, int restart, int maxit,
    void(*clbk)(void*));
template MATHLIB int GMRES (const SCMatrix &A, const SCVector &b, SCVector &x,
    double &tol, SCPreconditioner *precon, int restart, int maxit,
    void(*clbk)(void*));
#endif

template MATHLIB int GMRES (RVector (*Av_clbk)(const RVector &v, void *context),
    void *context, const RVector &b, RVector &x, double tol,
    RPreconditioner *precon, int restart, int maxit, int *iter, double *res);
template MATHLIB int GMRES (CVector (*Av_clbk)(const CVector &v, void *context),
    void *context, const CVector &b, CVector &x, double tol,
    CPreconditioner *precon, int restart, int maxit, int *iter, double *res);

#endif // NEED_EXPLICIT_INSTANTIATION
