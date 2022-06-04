// ========================================================================
// Implementation of Jacobian methods for script interfaces
// ========================================================================

#include "calc_mesh.h"

void BuildMesh(QMMesh **mesh, int nvtx, int nel, int dim, int nnd0,
               double *vtx, int *idx, int *eltp)

{
    int i, j, k;

    *mesh = new QMMesh;;
   
    // create node list
    (*mesh)->nlist.New (nvtx);
    for (i = 0; i < nvtx; i++) {
	(*mesh)->nlist[i].New(dim);
	(*mesh)->nlist[i].SetBndTp (BND_NONE); // don't know
    }
    for (j = k = 0; j < dim; j++) {
	for (i = 0; i < nvtx; i++) {
	    (*mesh)->nlist[i][j] = vtx[k++];
	}
    }

    // create element list
    Element *el, **list = new Element*[nel];
    for (i = 0; i < nel; i++) {
	switch (eltp[i]) {
	case ELID_TRI3OLD:
	    list[i] = new Triangle3old;
	    break;
	case ELID_TET4:
	    list[i] = new Tetrahedron4;
	    break;
	case ELID_WDG6:
	    list[i] = new Wedge6;
	    break;
	case ELID_VOX8:
	    list[i] = new Voxel8;
	    break;
	case ELID_TRI6:
	    list[i] = new Triangle6;
	    break;
	case ELID_TET10:
	    list[i] = new Tetrahedron10;
	    break;
	case ELID_TRI6_IP:
	    list[i] = new Triangle6_ip;
	    break;
	case ELID_TRI10:
	    list[i] = new Triangle10;
	    break;
	case ELID_TRI10_IP:
	    list[i] = new Triangle10_ip;
	    break;
	case ELID_TET10_IP:
	    list[i] = new Tetrahedron10_ip;
	    break;
	case ELID_PIX4:
	    list[i] = new Pixel4;
	    break;
	case ELID_TRI3:
	    list[i] = new Triangle3;
	    break;
	case ELID_TRI3D3:
	    list[i] = new Triangle3D3;
	    break;
	case ELID_TRI3D6:
	    list[i] = new Triangle3D6;
	    break;
	default:
	    std::cerr << "Element type not supported!" << std::endl;
	    list[i] = 0;
	    break;
	}
    }
    (*mesh)->elist.SetList (nel, list);
    delete []list;

    for (j = k = 0; j < nnd0; j++) {
	for (i = 0; i < nel; i++) {
	    if (el = (*mesh)->elist[i]) {
		if (j < el->nNode())
		    el->Node[j] = idx[k];
	    }
	    k++;
	}
    }

    // check mesh consistency
    if ((*mesh)->Shrink()) {
    	std::cerr << "warning: removed unused nodes" << std::endl;
    	nvtx = (*mesh)->nlen();
    }
    
    // set up mesh
    (*mesh)->Setup();

}


void CalcFields (QMMesh *mesh, Raster *raster,
    const CCompRowMatrix &qvec, const RVector &mua, const RVector &mus,
    const RVector &ref, double freq, char *solver, double tol,
    CVector **dphi)
{
    const double c0 = 0.3;
    int n    = mesh->nlen();
    int dim  = mesh->Dimension();
    int nQ   = qvec.nRows();

    CFwdSolver FWS (mesh, solver, tol);

    // Solution in mesh basis
    Solution msol(OT_NPARAM, n);
    msol.SetActive (OT_CMUA, true);
    msol.SetActive (OT_CKAPPA, true);

    // Set optical coefficients
    msol.SetParam (OT_CMUA, mua*c0/ref);
    msol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    RVector c2a(n);
    for (int i = 0; i < n; i++) {
	    c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));
    }
    msol.SetParam (OT_C2A, c2a);

    FWS.SetDataScaling (DATA_LOG);
    FWS.SetPrecon (PRECON_ICH);

    double omega = freq * 2.0*Pi*1e-6; // convert from MHz to rad

    // Calculate direct and adjoint fields
    FWS.Allocate ();
    FWS.Reset (msol, omega);

    // build the field vectors
    *dphi = new CVector[nQ];
    for (int i = 0; i < nQ; i++) {
        (*dphi)[i].New (n);
    }

    FWS.CalcFields (qvec, *dphi);

}


void CalcFields (QMMesh *mesh, Raster *raster,
    const RCompRowMatrix &qvec, const RVector &mua, const RVector &mus,
    const RVector &ref, double freq, char *solver, double tol,
    RVector **dphi)
{
    const double c0 = 0.3;
    int n    = mesh->nlen();
    int dim  = mesh->Dimension();
    int nQ   = qvec.nRows();

    RFwdSolver FWS (mesh, solver, tol);

    // Solution in mesh basis
    Solution msol(OT_NPARAM, n);
    msol.SetActive (OT_CMUA, true);
    msol.SetActive (OT_CKAPPA, true);

    // Set optical coefficients
    msol.SetParam (OT_CMUA, mua*c0/ref);
    msol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    RVector c2a(n);
    for (int i = 0; i < n; i++) {
	    c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));
    }
    msol.SetParam (OT_C2A, c2a);

    FWS.SetDataScaling (DATA_LOG);
    FWS.SetPrecon (PRECON_ICH);

    // Calculate direct and adjoint fields
    FWS.Allocate ();

    FWS.Reset (msol, 0);

    // build the field vectors
    *dphi = new RVector[nQ];
    for (int i = 0; i < nQ; i++) {
        (*dphi)[i].New (n);
    }

    FWS.CalcFields (qvec, *dphi);
}                                                    
