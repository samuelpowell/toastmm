#define RASTERLIB_IMPLEMENTATION
#include "rasterlib.h"

Raster_RampBlob::Raster_RampBlob (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, double _sup, RDenseMatrix *bb)
    : Raster_Blob (_bdim, _gdim, mesh, _sup, bb)
{
    isup = 1.0/sup;
    scale = 1.0;
    scale *= ComputeBasisScale();

    ComputeNodeValues();
}

double Raster_RampBlob::Value_nomask (const Point &p, int i, bool is_solidx)
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
    GetBasisIndices(ib,b);
    RVector r(dim);
    for (d = 0; d < dim; d++) {
	r[d] = (b[d]-npad) + (bbmin[d]-p[d])*igrid[d];
    }
    double rad = l2norm(r);
    if (rad >= sup) return 0.0; // p outside support of basis i
    return scale * (1.0 - rad*isup);
}
