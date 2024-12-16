// ========================================================================
// Implementation of Qvec/Mvec functions
// ========================================================================

#include "toastdef.h"
#include "calc_qmvec.h"

#ifdef TOAST_THREAD
#include "task.h"
#endif

// ============================================================================

void Integrate_Lin_Cosine (double d, double a, double x0, double x1,
    double &int_cos_u0, double &int_cos_u1)
{
    double arg1 = 2.0*a / (Pi*Pi*(x0-x1));

    int_cos_u0 = arg1 * (-2*a * cos(Pi*(d-x0)/(2*a)) +
			 2*a * cos(Pi*(d-x1)/(2*a)) +
			 Pi * (x0-x1) * sin(Pi*(d-x0)/(2*a)));
    int_cos_u1 = arg1 * (2*a * cos(Pi*(d-x0)/(2*a)) -
			 2*a * cos(Pi*(d-x1)/(2*a)) -
			 Pi * (x0-x1) * sin(Pi*(d-x1)/(2*a)));
}

// ----------------------------------------------------------------------------

CVector CompleteTrigSourceVector (const Mesh &mesh, int order)
{
    // currently only works with 2D circular mesh centered at origin
    int el, sd, nnode, *node;
    int n = mesh.nlen();
    double phi0, phi1, a, f0, f1;
    Element *pel;
    CVector qvec (n);

    for (el = 0; el < mesh.elen(); el++) {
	pel = mesh.elist[el];
	xASSERT (pel->Type() == ELID_TRI3, "Element type not supported");
	nnode = pel->nNode();
	node  = pel->Node;
	for (sd = 0; sd < pel->nSide(); sd++) {
	    if (!pel->IsBoundarySide (sd)) continue;
	    Node &nd0 = mesh.nlist[node[pel->SideNode (sd, 0)]];
	    Node &nd1 = mesh.nlist[node[pel->SideNode (sd, 1)]];
	    phi0 = atan2 (nd0[1], nd0[0]);
	    phi1 = atan2 (nd1[1], nd1[0]);

	    if (fabs (phi0-phi1) > Pi) {
		if (phi1 > phi0) phi0 += 2.0*Pi;
		else             phi1 += 2.0*Pi;
	    }
	    if (order) {
		a    = 2.0*Pi/4.0/order;
		Integrate_Lin_Cosine (0, a, phi0, phi1, f0, f1);
	    } else {
		f0 = f1 = 0.0;
	    }
	    f0 += fabs (phi1-phi0);
	    f1 += fabs (phi1-phi0);
	    qvec[node[pel->SideNode(sd,0)]] += f0;
	    qvec[node[pel->SideNode(sd,1)]] += f1;
	}
    }
    return qvec;
}

#ifdef TOAST_THREAD
struct Qvec_Threaddata {
    QMMesh *mesh;
    SourceMode qtype;
    SRC_PROFILE qprof;
    double qwidth;
    CCompRowMatrix *qvec;
};

void Qvec_engine (task_data *td)
{
    int i;
    int itask = td->proc;
    int ntask = td->np;
    Qvec_Threaddata *thdata = (Qvec_Threaddata*)td->data;
    QMMesh *mesh = thdata->mesh;
    int n = mesh->nlen();
    int nq = mesh->nQ;
    int q0 = (itask*nq)/ntask;
    int q1 = ((itask+1)*nq)/ntask;
    int dq = q1-q0;
    SourceMode qtype = thdata->qtype;
    SRC_PROFILE qprof = thdata->qprof;
    double qwidth = thdata->qwidth;
    CCompRowMatrix *qvec = thdata->qvec;

    CCompRowMatrix qvec_part(dq, n);

    for (i = q0; i < q1; i++) {
	CVector q(n);
	switch (qprof) {
	case PROF_POINT:
	    SetReal (q, QVec_Point (*mesh, mesh->Q[i], qtype));
	    break;
	case PROF_GAUSSIAN:
	    SetReal (q, QVec_Gaussian (*mesh, mesh->Q[i], qwidth, qtype));
	    break;
	case PROF_COSINE:
	    SetReal (q, QVec_Cosine (*mesh, mesh->Q[i], qwidth, qtype));
	    break;
	case PROF_COMPLETETRIG:
	    q = CompleteTrigSourceVector (*mesh, i);
	    break;
	}
	qvec_part.SetRow (i-q0, q);
    }

    Task::UserMutex_lock();
    qvec->SetRows (q0, qvec_part);
    Task::UserMutex_unlock();
}
#endif

void CalcQvec(QMMesh *mesh, SourceMode qtype,
              SRC_PROFILE qprof, double qwidth, CCompRowMatrix *qvec) {
  int n, nQ;

  n = mesh->nlen();
  nQ = mesh->nQ;

  // build the source vectors
  qvec->New(nQ, n);
  
#ifndef TOAST_THREAD
    for (int i = 0; i < nQ; i++) {
	CVector q(n);
	switch (qprof) {
	case PROF_POINT:
	    SetReal (q, QVec_Point (*mesh, mesh->Q[i], qtype));
	    break;
	case PROF_GAUSSIAN:
	    SetReal (q, QVec_Gaussian (*mesh, mesh->Q[i], qwidth, qtype));
	    break;
	case PROF_COSINE:
	    SetReal (q, QVec_Cosine (*mesh, mesh->Q[i], qwidth, qtype));
	    break;
	case PROF_COMPLETETRIG:
	    q = CompleteTrigSourceVector (*mesh, i);
	    break;
	}
	qvec.SetRow (i, q);
    }
#else
    Qvec_Threaddata thdata = {
	mesh,
	qtype,
	qprof,
	qwidth,
	qvec
    };
    Task::Multiprocess (Qvec_engine, &thdata);
#endif

}



#ifdef TOAST_THREAD
struct Mvec_Threaddata {
    QMMesh *mesh;
    SRC_PROFILE mprof;
    double mwidth;
    RVector *ref;
    CCompRowMatrix *mvec;
};

void Mvec_engine (task_data *td)
{
    const double c0 = 0.3;
    int i, j;
    int itask = td->proc;
    int ntask = td->np;
    Mvec_Threaddata *thdata = (Mvec_Threaddata*)td->data;
    QMMesh *mesh = thdata->mesh;
    RVector *ref = thdata->ref;
    int n = mesh->nlen();
    int nm = mesh->nM;
    int m0 = (itask*nm)/ntask;
    int m1 = ((itask+1)*nm)/ntask;
    int dm = m1-m0;
    SRC_PROFILE mprof = thdata->mprof;
    double mwidth = thdata->mwidth;
    CCompRowMatrix *mvec = thdata->mvec;

    CCompRowMatrix mvec_part(dm, n);

    for (i = m0; i < m1; i++) {
	CVector m(n);
	switch (mprof) {
	case PROF_GAUSSIAN:
	    SetReal (m, QVec_Gaussian (*mesh, mesh->M[i], mwidth,
				       SRCMODE_NEUMANN));
	    break;
	case PROF_COSINE:
	    SetReal (m, QVec_Cosine (*mesh, mesh->M[i], mwidth,
				     SRCMODE_NEUMANN));
	    break;
	case PROF_COMPLETETRIG:
	    m = CompleteTrigSourceVector (*mesh, i);
	    break;
	default:
		xERROR("Unknown source profile");
	}
	if (ref) {
	    RVector &rref = *ref;
	    for (j = 0; j < n; j++) m[j] *= c0/(2.0*rref[j]*A_Keijzer(rref[j]));
	}
	mvec_part.SetRow (i-m0, m);
    }

    Task::UserMutex_lock();
    mvec->SetRows (m0, mvec_part);
    Task::UserMutex_unlock();
}
#endif


void CalcMvec(QMMesh *mesh, SRC_PROFILE mprof, double mwidth,
              RVector *ref, bool apply_c2a, CCompRowMatrix *mvec) {
  int n, nM;

  n = mesh->nlen();
  nM = mesh->nM;

  // build the measurement vectors
  mvec->New(nM, n);
  
    #ifndef TOAST_THREAD
    const double c0 = 0.3;
    for (int i = 0; i < nM; i++) {
	CVector m(n);
	switch (mprof) {
	case PROF_GAUSSIAN:
	    SetReal (m, QVec_Gaussian (*mesh, mesh->M[i], mwidth,
				       SRCMODE_NEUMANN));
	    break;
	case PROF_COSINE:
	    SetReal (m, QVec_Cosine (*mesh, mesh->M[i], mwidth,
				     SRCMODE_NEUMANN));
	    break;
	case PROF_COMPLETETRIG:
	    m = CompleteTrigSourceVector (*mesh, i);
	    break;
	}
	if (apply_c2a)
	    for (int j = 0; j < n; j++) m[j] *= c0/(2.0*ref[j]*A_Keijzer(ref[j]));
	mvec.SetRow (i, m);
    }
#else
    Mvec_Threaddata thdata = {
	mesh,
	mprof,
	mwidth,
	(apply_c2a ? ref : 0),
	mvec
    };
    Task::Multiprocess (Mvec_engine, &thdata);
#endif

}

