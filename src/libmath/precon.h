// -*-C++-*-
// ==========================================================================
// Module mathlib
// File precon.h
// Declaration of tmplate class TPreconditioner
// Preconditioner for iterative matrix solution methods
// ==========================================================================

#ifndef __PRECON_H
#define __PRECON_H

#include "vector.h"
#include "matrix.h"
#include "dgmatrix.h"
#include "crmatrix.h"

/// \defgroup iterative linear solver preconditioner types
//@{
typedef enum {
    PRECON_NULL,                     ///< no preconditioner
    PRECON_DIAG,                     ///< diagonal (Jacobi) preconditioner
    PRECON_ICH,                      ///< incomplete Choleski
    PRECON_DILU,                     ///< diagonal incomplete LU
    PRECON_CG_MULTIGRID,             ///< multigrid
    PRECON_ILU,                      ///< incomplete LU
} PreconType;
//@}

// ==========================================================================
// class TPreconditioner

template<class MT> class TPreconditioner {
public:
    TPreconditioner() {}
    virtual ~TPreconditioner() {}

    virtual PreconType Type() = 0;

    virtual void Reset (const TMatrix<MT> *A) = 0;
    // Reset preconditioner for matrix A
    
    virtual void Apply (const TVector<MT> &r, TVector<MT> &s)=0;
    // Apply preconditioner to r and return the result in s
    // e.g. s = M^-1 r for a preconditioner matrix M
    virtual void Apply (const TDenseMatrix<MT> &r, TDenseMatrix<MT> &s) const =0;
       

    static TPreconditioner *Create (PreconType type);
    // create a new preconditioner of type 'type' and return pointer to it
};

// ==========================================================================
// class TPrecon_Null
// Dummy preconditioner (M = I)

template<class MT> class TPrecon_Null: public TPreconditioner<MT> {
public: 
    TPrecon_Null() {}
    PreconType Type() { return PRECON_NULL; }
    void Reset (const TMatrix<MT>*) {}
    void Apply (const TVector<MT> &r, TVector<MT> &s) { s = r; }
    void Apply (const TDenseMatrix<MT> &r, TDenseMatrix<MT> &s) const
    { ERROR_UNDEF; };

};

// ==========================================================================
// class TPrecon_Diag
// diagonal preconditioner (M = diag(A))

template<class MT> class TPrecon_Diag: public TPreconditioner<MT> {
public:
    TPrecon_Diag() {}
    PreconType Type() { return PRECON_DIAG; }
    void Reset (const TMatrix<MT> *A);
    void ResetFromDiagonal (const TVector<MT> &diag);
    void Apply (const TVector<MT> &r, TVector<MT> &s);
    void Apply (const TDenseMatrix<MT> &r, TDenseMatrix<MT> &s) const;

private:
    TVector<MT> idiag;
};

// ==========================================================================
// class TPrecon_IC
// Incomplete Cholesky preconditioner

template<class MT> class TPrecon_IC: public TPreconditioner<MT> {
public:
    TPrecon_IC() {}
    PreconType Type() { return PRECON_ICH; }
    void Reset (const TMatrix<MT> *A);
    void Apply (const TVector<MT> &r, TVector<MT> &s);
    void Apply (const TDenseMatrix<MT> &r, TDenseMatrix<MT> &s) const
    { ERROR_UNDEF; };

private:
    TCompRowMatrix<MT> L;
    TVector<MT>d;
};

// ==========================================================================
// class TPrecon_DILU
// D-ILU: Simple incomplete LU preconditioner where only diagonal elements
// are modified

template<class MT> class TPrecon_DILU: public TPreconditioner<MT> {
public:
    TPrecon_DILU() {}
    PreconType Type() { return PRECON_DILU; }
    void Reset (const TMatrix<MT> *);
    void Apply (const TVector<MT> &r, TVector<MT> &s);
    void Apply (const TDenseMatrix<MT> &r, TDenseMatrix<MT> &s) const
    { ERROR_UNDEF; };

private:
    int dim;                   // problem dimension
    TVector<MT>ipivot;         // inverse pivots
    const TMatrix<MT> *A;  // pointer to matrix
};


// ==========================================================================
// class TPrecon_CG_Multigrid
// Uses multigrid CG solution as preconditioner
// Suitable for SPD systems

template<class MT> class TPrecon_CG_Multigrid: public TPreconditioner<MT> {
public:
    TPrecon_CG_Multigrid (const TCompRowMatrix<MT> *AA,
        TCompRowMatrix<MT> **PP, TCompRowMatrix<MT> **RR,
        TPreconditioner<MT> **pprecon, int nlvl, int nngs = 2, int nnmg = 2) : 
      A(AA), P(PP), R(RR), precon(pprecon), maxlevel(nlvl), ngs(nngs),
      nmg(nnmg) {}

    PreconType Type() { return PRECON_CG_MULTIGRID; }

    void Reset (const TMatrix<MT> *);
    // _A points to an *array* of nlvl system matrices, from finest to
    // coarsest grid

    void Apply (const TVector<MT> &r, TVector<MT> &s);
    void Apply (const TDenseMatrix<MT> &r, TDenseMatrix<MT> &s) const
    { ERROR_UNDEF; };


private:
    const TCompRowMatrix<MT> *A; // pointer to system matrix array
    TCompRowMatrix<MT> **P;      // interpolator matrix array
    TCompRowMatrix<MT> **R;      // restrictor matrix array
    TPreconditioner<MT> **precon; // coarse-level CG preconditioner
    int maxlevel, ngs, nmg;
    double MGM (const TVector<MT> &r, TVector<MT> &x, int level) const;
    void Smoother (const TCompRowMatrix<MT> &A, const TVector<MT> &r,
        TVector<MT> &x, int itermax) const;
};

// ==========================================================================
// template typedefs

typedef TPreconditioner<double>         RPreconditioner;
typedef TPreconditioner<std::complex<double> > CPreconditioner;
typedef TPreconditioner<int>            IPreconditioner;
#ifdef TOAST_FEATURE_SINGLEPREC
typedef TPreconditioner<float>          FPreconditioner;
typedef TPreconditioner<std::complex<float> >       SCPreconditioner;
#endif

typedef TPrecon_Null<double>            RPrecon_Null;
typedef TPrecon_Null<std::complex<double> >    CPrecon_Null;
typedef TPrecon_Null<int>               IPrecon_Null;
#ifdef TOAST_FEATURE_SINGLEPREC
typedef TPrecon_Null<float>             FPrecon_Null;
typedef TPrecon_Null<std::complex<float> >          SCPrecon_Null;
#endif

typedef TPrecon_Diag<double>            RPrecon_Diag;
typedef TPrecon_Diag<std::complex<double> >    CPrecon_Diag;
typedef TPrecon_Diag<int>               IPrecon_Diag;
#ifdef TOAST_FEATURE_SINGLEPREC
typedef TPrecon_Diag<float>             FPrecon_Diag;
typedef TPrecon_Diag<std::complex<float> >          SCPrecon_Diag;
#endif

typedef TPrecon_IC<double>              RPrecon_IC;
typedef TPrecon_IC<std::complex<double> >      CPrecon_IC;
typedef TPrecon_IC<int>                 IPrecon_IC;
#ifdef TOAST_FEATURE_SINGLEPREC
typedef TPrecon_IC<float>               FPrecon_IC;
typedef TPrecon_IC<std::complex<float> >            SCPrecon_IC;
#endif

typedef TPrecon_DILU<double>            RPrecon_DILU;
typedef TPrecon_DILU<std::complex<double> >    CPrecon_DILU;
typedef TPrecon_DILU<int>               IPrecon_DILU;
#ifdef TOAST_FEATURE_SINGLEPREC
typedef TPrecon_DILU<float>             FPrecon_DILU;
typedef TPrecon_DILU<std::complex<float> >          SCPrecon_DILU;
#endif

typedef TPrecon_CG_Multigrid<double>    RPrecon_CG_Multigrid;
typedef TPrecon_CG_Multigrid<std::complex<double> > CPrecon_CG_Multigrid;
typedef TPrecon_CG_Multigrid<int>       IPrecon_CG_Multigrid;
#ifdef TOAST_FEATURE_SINGLEPREC
typedef TPrecon_CG_Multigrid<float>     FPrecon_CG_Multigrid;
typedef TPrecon_CG_Multigrid<std::complex<float> >  SCPrecon_CG_Multigrid;
#endif

#endif // !__PRECON_H
