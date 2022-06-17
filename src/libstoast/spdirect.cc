#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"


// =========================================================================
// Nonmember declarations

cholmod_sparse cholmod_wrap(const RCompRowMatrix &F, bool pattern = false)
{
    cholmod_sparse A;
    A.nrow = F.nRows();
    A.ncol = F.nCols();
    A.nzmax = F.nVal();
    A.p = F.rowptr;
    A.i = F.colidx;
    A.nz = 0;
    A.x = pattern ? 0 : (void *) F.ValPtr();
    A.z = 0;
    A.stype = 0;
    A.itype = CHOLMOD_INT;   
    A.xtype = pattern ? CHOLMOD_PATTERN : CHOLMOD_REAL;
    A.dtype = CHOLMOD_DOUBLE;
    A.sorted = 1;
    A.packed = 1;
    return A;
}

cholmod_sparse cholmod_wrap(const FCompRowMatrix &F,bool pattern = false)
{
    cholmod_sparse A;
    A.nrow = F.nRows();
    A.ncol = F.nCols();
    A.nzmax = F.nVal();
    A.p = F.rowptr;
    A.i = F.colidx;
    A.nz = 0;
    A.x = pattern ? 0 : (void *) F.ValPtr();
    A.z = 0;
    A.stype = 0;
    A.itype = CHOLMOD_INT;   
    A.xtype = pattern ? CHOLMOD_PATTERN : CHOLMOD_REAL;
    A.dtype = CHOLMOD_SINGLE;
    A.sorted = 1;
    A.packed = 1;
    return A;
}

cholmod_dense cholmod_wrap(const RVector &v)
{
    cholmod_dense b;
    b.nrow = v.Dim();
    b.ncol = 1;
    b.nzmax = b.nrow;
    b.d = b.nrow;
    b.x = (void *) v.data_buffer();
    b.z = 0;
    b.xtype = CHOLMOD_REAL;
    b.dtype = CHOLMOD_DOUBLE;
    return b;
}

cholmod_dense cholmod_wrap(const FVector &v)
{
    cholmod_dense b;
    b.nrow = v.Dim();
    b.ncol = 1;
    b.nzmax = b.nrow;
    b.d = b.nrow;
    b.x = (void *) v.data_buffer();
    b.z = 0;
    b.xtype = CHOLMOD_REAL;
    b.dtype = CHOLMOD_SINGLE;
    return b;
}


// =========================================================================

template<>
TSPDirectSolver<std::complex<float> >::TSPDirectSolver ()
{
    #ifndef TOAST_NO_GPL
    umfpack_zi_defaults(up_c);
    #else
    ERROR;
    #endif
}


template<>
TSPDirectSolver<std::complex<double> >::TSPDirectSolver ()
{
    #ifndef TOAST_NO_GPL
    umfpack_zi_defaults(up_c);
    #else
    ERROR;
    #endif
}

template<>
TSPDirectSolver<float>::TSPDirectSolver ()
{
    cm_L = 0;
    cholmod_start (&cm_c);
    #ifdef TOAST_NO_GPL
    cm_c.supernoal = CHOLMOD_SIMPLICIAL;
    #endif
}

template<>
TSPDirectSolver<double>::TSPDirectSolver ()
{
    cm_L = 0;
    cholmod_start (&cm_c);
    #ifdef TOAST_NO_GPL
    cm_c.supernoal = CHOLMOD_SIMPLICIAL;
    #endif
}

// =========================================================================

template<>
TSPDirectSolver<float>::~TSPDirectSolver ()
{
    // Free CHOLMOD factorisation object
    cholmod_free_factor(&cm_L, &cm_c);
    
    // Free CHOLMOD solve workspace 
    // cholmod_free_dense (& (cholmod_dense *)cm_X, &cm_c);
    // cholmod_free_dense (& (cholmod_dense *)cm_Y, &cm_c);
    // cholmod_free_dense (& (cholmod_dense *)cm_E, &cm_c);
    cholmod_finish (&cm_c) ;   
}

template<>
TSPDirectSolver<double>::~TSPDirectSolver ()
{
    // Free CHOLMOD factorisation object
    cholmod_free_factor(&cm_L, &cm_c);  

    // Free CHOLMOD solve workspace 
    // cholmod_free_dense (& (cholmod_dense *)cm_X, &cm_c);
    // cholmod_free_dense (& (cholmod_dense *)cm_Y, &cm_c);
    // cholmod_free_dense (& (cholmod_dense *)cm_E, &cm_c);
    cholmod_finish (&cm_c) ;    
}

template<>
TSPDirectSolver<std::complex<float> >::~TSPDirectSolver ()
{
    // Free UMFPACK factorisations
    // umfpack_zi_free_numeric(&up_numeric);
    // umfpack_zi_free_symbolic(&up_symbolic);
}

template<>
TSPDirectSolver<std::complex<double> >::~TSPDirectSolver ()
{
    // Free UMFPACK factorisations
    // umfpack_zi_free_numeric(&up_numeric);
    // umfpack_zi_free_symbolic(&up_symbolic);
}

// =========================================================================

template<>
bool TSPDirectSolver<float>::factorise_symbolic (const FCompRowMatrix &F) 
{
    cholmod_sparse A = cholmod_wrap(F, true);
    cholmod_free_factor(&cm_L, &cm_c);
    cm_L = cholmod_analyze (&A, &cm_c);
    xASSERT(cm_L, "System matrix analysis failed!");
    return cm_L != NULL;
}

template<>
bool TSPDirectSolver<double>::factorise_symbolic (const RCompRowMatrix &F) 
{
    cholmod_sparse A = cholmod_wrap(F, true);
    cholmod_free_factor(&cm_L, &cm_c);
    cm_L = cholmod_analyze (&A, &cm_c);
    xASSERT(cm_L, "System matrix analysis failed!");
    return cm_L != NULL;
}

template<>
bool TSPDirectSolver<std::complex<float> >::factorise_symbolic (const SCCompRowMatrix &F) 
{
    int n = F.nRows();
    CVector dcv(n);
    for(int i = 0; i < n; ++i)
    {
        dcv[i] = std::complex<double>(*(F.ValPtr() +i ));
    }     
    // umfpack_zi_free_numeric(&up_numeric);
    // umfpack_zi_free_symbolic(&up_symbolic);
    int status = umfpack_zi_symbolic(n, n, F.rowptr, F.colidx, (const double *)dcv.data_buffer(),
                                     NULL, &up_symbolic, up_c, NULL);
    xASSERT(status >= 0, "System matrix analysis failed!");
    return status == 0;
}

template<>
bool TSPDirectSolver<std::complex<double> >::factorise_symbolic (const CCompRowMatrix &F) 
{
    int n = F.nRows();
    // umfpack_zi_free_numeric(&up_numeric);
    // umfpack_zi_free_symbolic(&up_symbolic);
    int status = umfpack_zi_symbolic(n, n, F.rowptr, F.colidx, (const double *) F.ValPtr(), 
                                     NULL, &up_symbolic, up_c, NULL);
    xASSERT(status >= 0, "System matrix analysis failed!");
    return status == 0;
}


// =========================================================================

template<>
bool TSPDirectSolver<float>::factorise_numeric (const FCompRowMatrix &F) 
{
    cholmod_sparse A = cholmod_wrap(F);
    int status = cholmod_factorize(&A, cm_L, &cm_c);
    // SP TODO: Check return and the common structure for failure
    return true;
}

template<>
bool TSPDirectSolver<double>::factorise_numeric (const RCompRowMatrix &F) 
{
    cholmod_sparse A = cholmod_wrap(F);
    int status = cholmod_factorize(&A, cm_L, &cm_c);
    // SP TODO: Check return and the common structure for failure
    return true;
}

template<>
bool TSPDirectSolver<std::complex<float> >::factorise_numeric (const SCCompRowMatrix &F) 
{
    int n = F.nCols();
    CVector dcv(n);
    for (int i = 0; i < n; ++i)
    {
        dcv[i] = std::complex<double>(*(F.ValPtr() + i));
    }
    // umfpack_zi_free_numeric(&up_numeric);
    int status = umfpack_zi_numeric(F.rowptr, F.colidx, (const double *)dcv.data_buffer(), NULL,
                                    up_symbolic, &up_numeric, up_c, NULL);
    xASSERT(status >= 0, "System matrix factorisation failed");
    return status == 0;
}

template<>
bool TSPDirectSolver<std::complex<double> >::factorise_numeric (const CCompRowMatrix &F) 
{
    // umfpack_zi_free_numeric(&up_numeric);
    int status = umfpack_zi_numeric(F.rowptr, F.colidx, (const double *)F.ValPtr(), NULL,
                                    up_symbolic, &up_numeric, up_c, NULL);
    xASSERT(status >= 0, "System matrix factorisation failed");
    return status == 0;
}

// =========================================================================

template<>
bool TSPDirectSolver<float>::solve (const FCompRowMatrix &F, FVector &x, const FVector &b) 
{
    cholmod_dense cm_b = cholmod_wrap(b);
    cholmod_dense *cm_x = cholmod_solve(CHOLMOD_A, cm_L, &cm_b, (cholmod_common *)&cm_c);

    // SP TODO: Use cholmod solve 2 to avoid this
    float *xd = (float *)cm_x->x;
    for (int i = 0; i < x.Dim(); i++)
    {
        x[i] = xd[i];
    }
    cholmod_free_dense(&cm_x, (cholmod_common *)&cm_c);

    // cholmod_dense b = cholmod_wrap(qvec);
    // cholmod_solve2(CHOLMOD_A, cm_L, &b, NULL, & (cholmod_dense *)cm_X, NULL, &(cholmod_dense *)cm_Y, &(cholmod_dense *)cm_E, &(cholmod_common)cm_c);
    // double *xd = (double *) cm_X->x;
    // for(int i = 0; i < phi.Dim(); i++) {
    //     phi[i] = xd[i];
    // }

    return true;
}


template<>
bool TSPDirectSolver<double>::solve (const RCompRowMatrix &F, RVector &x, const RVector &b) 
{
    cholmod_dense cm_b = cholmod_wrap(b);
    cholmod_dense *cm_x = cholmod_solve(CHOLMOD_A, cm_L, &cm_b, (cholmod_common *)&cm_c);

    // SP TODO: Use cholmod solve 2 to avoid this
    double *xd = (double *)cm_x->x;
    for (int i = 0; i < x.Dim(); i++)
    {
        x[i] = xd[i];
    }
    cholmod_free_dense(&cm_x, (cholmod_common *)&cm_c);

    // cholmod_dense b = cholmod_wrap(qvec);
    // cholmod_solve2(CHOLMOD_A, cm_L, &b, NULL, & (cholmod_dense *)cm_X, NULL, &(cholmod_dense *)cm_Y, &(cholmod_dense *)cm_E, &(cholmod_common)cm_c);
    // double *xd = (double *) cm_X->x;
    // for(int i = 0; i < phi.Dim(); i++) {
    //     phi[i] = xd[i];
    // }

    return true;
}

template<>
bool TSPDirectSolver<std::complex<float> >::solve (const SCCompRowMatrix &F, SCVector &x, const SCVector &b) 
{
    int n = F.nRows();
    CVector ux(n);

    CVector ub(n);
    for (int i = 0; i < n; i++)
    {
        ub[i] = std::complex<double>(b[i]);
    }

    CVector uf(n);
    for (int i = 0; i < n; ++i)
    {
        uf[i] = std::complex<double>(*(F.ValPtr() + i));
    }

    int status = umfpack_zi_solve(UMFPACK_A, F.rowptr, F.colidx, (const double *)uf.data_buffer(), NULL,
                                    (double *)ux.data_buffer(), NULL,
                                    (const double *)ub.data_buffer(), NULL,
                                    up_numeric, up_c, NULL);
    xASSERT(status >= 0, "Sytem matrix solve failed");

    for (int i = 0; i < n; i++)
    {
        x[i] = std::complex<float>(ux[i]);
    }
    return status == 0;
}

template<>
bool TSPDirectSolver<std::complex<double> >::solve (const CCompRowMatrix &F, CVector &x, const CVector &b) 
{
    int status = umfpack_zi_solve(UMFPACK_A, F.rowptr, F.colidx, (const double *)F.ValPtr(), NULL,
                                      (double *)x.data_buffer(), NULL,
                                      (const double *)b.data_buffer(), NULL,
                                      up_numeric, up_c, NULL);
    xASSERT(status >= 0, "Sytem matrix solve failed");
    return status == 0;
}

// =========================================================================

#ifdef NEED_EXPLICIT_INSTANTIATION

template class STOASTLIB TSPDirectSolver<float>;
template class STOASTLIB TSPDirectSolver<double>;
template class STOASTLIB TSPDirectSolver<std::complex<double> >;
template class STOASTLIB TSPDirectSolver<std::complex<float> >;

#endif