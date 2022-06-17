// =========================================================================
// spsolver.h 
// Interface for sparse direct solvers
// =========================================================================

#ifndef __DIRECTSOLVER_H
#define __DIRECTSOLVER_H

#include "mathlib.h"
#include "stoastlib.h"

#include "cholmod.h"
#include "umfpack.h"


// =========================================================================

template<class T> class TSPDirectSolver;

template<class T> class TSPDirectSolver {

public:
    /**
     * \brief Constructor. Creates a sparse direct forward solver instance.
     */
    TSPDirectSolver ();

    /**
     * \brief Destructor. Destroys the forward solver instance.
     */
    ~TSPDirectSolver();

    /**
     * \brief Perform a symbolic factorisation of the supplied symmetric sparse matrix
     * \return True on success, else false
     */
    bool factorise_symbolic(const TCompRowMatrix<T> &F);

    /**
     * \brief Perform a numeric factorisation of the supplied symmetric sparse matrix, which
     *        must have previously undregone symolic factorisation.
     * \return True on success, else false
     */
    bool factorise_numeric(const TCompRowMatrix<T> &F);

    /**
     * \brief Solve Ax=b for a single RHS b, using the factors stored in the direct solver
     *        resulting from a previous call to the numeric factoisation function. The input
     *        solution vector x is modified in place.
     * \return True on success, else false
     */
    bool solve(const TCompRowMatrix<T> &F, TVector<T> &x, const TVector<T> &b);

private:

    // UMFPACK data for complex direct solves
    double up_c[UMFPACK_CONTROL];
    void *up_symbolic;
    void *up_numeric;

    // CHOLMOD data for real direct solves
    cholmod_common  cm_c;               // CHOLMOD common data structure
    cholmod_factor *cm_L;               // CHOLMOD factor


};



// ==========================================================================
// template typedefs

typedef TSPDirectSolver<float>          FSPDirectSolver;
typedef TSPDirectSolver<double>         RSPDirectSolver;
typedef TSPDirectSolver<std::complex<float> > SCSPDirectSolver;
typedef TSPDirectSolver<std::complex<double> > CSPDirectSolver;


#endif