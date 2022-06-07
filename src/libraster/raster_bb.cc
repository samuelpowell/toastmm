#define RASTERLIB_IMPLEMENTATION
#include "rasterlib.h"

// ============================================================================

static double bessi0 (double x)
{
    double ax, ans, y;

    if ((ax = fabs(x)) < 3.75) {
        y = x / 3.75;
	y *= y;
	ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 +
		    y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
    } else {
        y = 3.75/ax;
	ans=(exp(ax)/sqrt(ax)) * (0.39894228 + y*(0.1328592e-1 +
	     y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2 +
             y*(-0.2057706e-1 + y*(0.2635537e-1 + y*(-0.1647633e-1 +
	     y*0.392377e-2))))))));
    }
    return ans;
}

// ============================================================================

static double bessi (int n, double x)
{
    const double ACC = 40.0;
    const double BIGNO = 1e10;
    const double BIGNI = 1e-10;

    double tox, bi, bip, bim, ans;
    int j;

    if (n < 2) {
        cerr << "Index n less than 2 in bessi" << endl;
	exit (1);
    }
    if (x == 0.0)
        return 0.0;
    else {
        tox = 2.0/fabs(x);
	bip = ans = 0.0;
	bi = 1.0;
	for (j = 2*(n+(int)sqrt(ACC*n)); j > 0; j--) {
	    bim = bip + j*tox * bi;
	    bip = bi;
	    bi = bim;
	    if (fabs(bi) > BIGNO) {
	        ans *= BIGNI;
		bi *= BIGNI;
		bip *= BIGNI;
	    }
	    if (j == n) ans = bip;
	}
	ans *= bessi0 (x) / bi;
	return x < 0.0 && (n & 1) ? -ans : ans;
    }
}


Raster_BesselBlob::Raster_BesselBlob (const IVector &_bdim,
    const IVector &_gdim, Mesh *mesh, double _alpha, double _sup,
    RDenseMatrix *bb)
    : Raster_Blob (_bdim, _gdim, mesh, _sup, bb), alpha(_alpha)
{
    a2 = sup*sup;
    scale = 1.0/bessi(2,alpha);
    scale *= ComputeBasisScale();

    ComputeNodeValues();
}

double Raster_BesselBlob::Value_nomask (const Point &p, int i, bool is_solidx)
    const
{
    int ib;
    if (is_solidx) {
	RANGE_CHECK (i >= 0 && i < slen);
	ib = GetBasisIdx(i);
    } else {
	RANGE_CHECK (i >= 0 && i < blen);
	ib = i;
    }
    int d;
    IVector b;
    GetBasisIndices(ib, b);
    RVector r(dim);
    for (d = 0; d < dim; d++) {
	r[d] = (b[d]-npad) + (bbmin[d]-p[d])*igrid[d];
    }
    double r2 = l2normsq(r);
    if (r2 >= a2) return 0.0; // p outside support of basis i
    double s = sqrt (1.0 - r2/a2);
    return scale * s*s * bessi(2,alpha*s);
}
