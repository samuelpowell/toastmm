#define RASTERLIB_IMPLEMENTATION
#include "rasterlib.h"
#include "raster_bb2.h"


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


Raster_Blob2_BB::Raster_Blob2_BB (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, double _sup, double _alpha, double diagscale,
    RDenseMatrix *bb, double _map_tol, int _npad)
: Raster_Blob2 (_bdim, _gdim, mesh, _sup, _alpha, diagscale, bb, _map_tol,
		_npad)
{
    a2 = sup*sup;
}

double Raster_Blob2_BB::RadValue (double r) const
{
    if (r >= sup) return 0.0;

    double r2 = r*r;
    double s = sqrt(1.0 - r2/a2);
    return s*s * bessi(2,sprm*s);
}
