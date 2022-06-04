// -*- C++ -*-
#include <Python.h>
#include <numpy/arrayobject.h>

#include <fstream>
#include <iostream>

#include "objmgr.h"
#include "toastdef.h"
#include "mathlib.h"
#include "felib.h"
#include "stoastlib.h"

// Computation routines
#include "calc_gradient.h"
#include "calc_jacobian.h"
#include "calc_mesh.h"

#define TOAST_NPY_INT NPY_INT   // The numpy type for integer output to Python
typedef int nint;               // The C type for integer input from Python

static_assert (sizeof(npy_int) == sizeof(nint),  "Numpy integer size mismatch");

struct module_state {
  PyObject *error;
};

#define GETSTATE(m) ((struct module_state *)PyModule_GetState(m))

#define GETMESH(mesh, hmesh)                                    \
  {                                                             \
    if (!(mesh = (QMMesh *)g_meshmgr.Get(hmesh))) {             \
      PyErr_SetString(PyExc_ValueError, "Invalid mesh handle"); \
      return NULL;                                              \
    }                                                           \
  }

#define GETRASTER(raster, hraster)                                \
  {                                                               \
    if (!(raster = (Raster *)g_rastermgr.Get(hraster))) {         \
      PyErr_SetString(PyExc_ValueError, "Invalid raster handle"); \
      return NULL;                                                \
    }                                                             \
  }

#define GETREGUL(reg, hreg)                                               \
  {                                                                       \
    if (!(reg = (Regularisation *)g_regmgr.Get(hreg))) {                  \
      PyErr_SetString(PyExc_ValueError, "Invalid regularisation handle"); \
      return NULL;                                                        \
    }                                                                     \
  }

// ===========================================================================
// A python wrapper object that handled deallocation of C++ allocated memory
// for setting up a scipy CSR matrix

typedef struct
{
  PyObject_HEAD int *rp;
  int *ci;
  void *val;
} CSR_mem_wrapper;

// ===========================================================================
// Globals

ObjectManager<Mesh> g_meshmgr;
ObjectManager<Raster> g_rastermgr;
ObjectManager<Regularisation> g_regmgr;

// ===========================================================================
// Helper functions

// Assert a PyObject to be an array of specifed rank and type
inline bool AssertArray(PyArrayObject *arr, int ndim, int type) {

  int rank = PyArray_NDIM(arr);
  if(rank != ndim) {
    std::cerr << "Array rank mismatch, req: " << ndim << " got: " << rank << std::endl;
    return false;
  }
  
  int tp = PyArray_TYPE(arr);
  if(!PyArray_EquivTypenums(tp, type)) {
     std::cerr << "Array type mismatch, req: " << type << " got: " << tp << std::endl;
     return false;
  }

  return true;
}

// Assert a PyObect to be an array of specified rank, type, and dimensions
inline bool AssertArrayDims(PyArrayObject *arr, int type, int dim0) {
  if (!AssertArray(arr, 1, type)) {
    return false;
  }
  npy_intp *dims = PyArray_DIMS(arr);
  if (dims[0] != dim0) {
    std::cerr << "Array dimension mismatch, reqd: " << dim0 << " got: " << dims[0] << std::endl;
    return false;
  }
  return true;
}

// Assert a PyObect to be an array of specified rank, type, and dimensions
inline bool AssertArrayDims(PyArrayObject *arr, int type, int dim0, int dim1) {
  if (!AssertArray(arr, 2, type)) {
    return false;
  }
  npy_intp *dims = PyArray_DIMS(arr);
  if ((dims[0] != dim0) || (dims[1] != dim1)) {
    std::cerr << "Array dimension mismatch, reqd: [" << dim0 << ", " << dim1 << "] got: [" 
              << dims[0] << ", " << dims[1] << "]" << std::endl;
    return false;
  }
  return true;
}

// Copy a row or column vector from python to toast, assuming that this is a
// one-dimensional ndarray.
RVector CopyVector(PyObject *pyvec) {
  npy_intp *dims = PyArray_DIMS((PyArrayObject *)pyvec);
  int dim = dims[0];
  RVector v(dim);

  if (PyArray_TYPE((PyArrayObject *)pyvec) != NPY_DOUBLE) {
    std::cerr << "Attempt to form RVector from incorrect data type" << std::endl;
  } else {
    memcpy(v.data_buffer(), PyArray_DATA((PyArrayObject *)pyvec), dim * sizeof(double));
  }

  return v;
}

// Copy a toast vector to a 1D python array object, output will be a rand-1 array
// of doubles.
void CopyVector(PyObject **pyvec, const RVector &vec) {
  npy_intp dim = vec.Dim();
  *pyvec = PyArray_SimpleNew(1, &dim, NPY_DOUBLE);
  double *pydata = (double *)PyArray_DATA((PyArrayObject *)*pyvec);
  const double *data = vec.data_buffer();
  memcpy(pydata, data, dim * sizeof(double));
}

// ===========================================================================

static PyObject *mesh_read(PyObject *self, PyObject *args) {
  const char *meshname;

  if (!PyArg_ParseTuple(args, "s", &meshname)) {
    return NULL;
  }

  QMMesh *mesh = new QMMesh;
  std::ifstream ifs(meshname);
  if(!ifs) {
    PyErr_SetString(PyExc_ValueError, "Mesh file could not be opened");
    return NULL;
  }
  ifs >> *mesh;

  mesh->Setup();
  int hmesh = g_meshmgr.Add(mesh);

  return Py_BuildValue("i", hmesh);
}

// ===========================================================================

static PyObject *mesh_write(PyObject *self, PyObject *args) {
  const char *meshname;
  int hmesh;
  QMMesh *mesh;

  if (!PyArg_ParseTuple(args, "is", &hmesh, &meshname)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);
  std::ofstream ofs(meshname);
  if(!ofs) {
    PyErr_SetString(PyExc_ValueError, "Mesh file could not be opened");
    return NULL;
  }
  ofs << *mesh;

  Py_RETURN_NONE;
}

// ===========================================================================

static PyObject *mesh_clear(PyObject *self, PyObject *args) {
  int hmesh;

  if (!PyArg_ParseTuple(args, "i", &hmesh)) {
    return NULL;
  }

  g_meshmgr.Delete(hmesh);
  Py_RETURN_NONE;
}

// ===========================================================================

static PyObject *mesh_data(PyObject *self, PyObject *args) {
  int hmesh;
  QMMesh *mesh;

  if (!PyArg_ParseTuple(args, "i", &hmesh)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  int i, j;
  int nlen = mesh->nlen();
  int elen = mesh->elen();
  int dim = mesh->Dimension();
  npy_intp node_dims[2] = {nlen, dim};

  PyObject *nodelist = PyArray_SimpleNew(2, node_dims, NPY_DOUBLE);
  double *v, *vtx_data = (double *)PyArray_DATA((PyArrayObject *)nodelist);

  for (i = 0, v = vtx_data; i < nlen; i++) {
    for (j = 0; j < dim; j++) {
      *v++ = mesh->nlist[i][j];
    }
  }

  // max number of nodes per element
  int nnd = mesh->elist[0]->nNode();
  for (i = 1; i < elen; i++) {
    nnd = std::max(nnd, mesh->elist[i]->nNode());
  }
  npy_intp el_dims[2] = {elen, nnd};

  PyObject *idx = PyArray_SimpleNew(2, el_dims, TOAST_NPY_INT);
  nint *e, *el_data = (nint *)PyArray_DATA((PyArrayObject *)idx);

  // element index list
  // (0-based; value -1 indicates unused matrix entry)
  for (i = 0, e = el_data; i < elen; i++) {
    for (j = 0; j < mesh->elist[i]->nNode(); j++) {
      *e++ = mesh->elist[i]->Node[j];
    }
    for (; j < nnd; j++) {
      *e++ = -1;
    }
  }

  // element type list (see element.h)
  PyObject *eltp = PyArray_SimpleNew(1, el_dims, TOAST_NPY_INT);
  nint *et, *etp_data = (nint *)PyArray_DATA((PyArrayObject *)eltp);

  for (i = 0, et = etp_data; i < elen; i++) {
    *et++ = mesh->elist[i]->Type();
  }

  return Py_BuildValue("OOO", nodelist, idx, eltp);
}

// ===========================================================================

static PyObject *toast_surf_data(PyObject *self, PyObject *args) {
  int hmesh;
  Mesh *mesh;

  if (!PyArg_ParseTuple(args, "i", &hmesh)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  int i, j, k;
  int nlen = mesh->nlen();
  int dim = mesh->Dimension();
  int nbnd = mesh->nlist.NumberOf(BND_ANY);
  int *bndidx = new int[nlen];
  npy_intp node_dims[2] = {nbnd, dim};

  double *v, *vtx_data = new double[nbnd * dim];
  for (i = 0, v = vtx_data; i < nlen; i++) {
    if (mesh->nlist[i].isBnd()) {
      for (j = 0; j < dim; j++) {
        *v++ = mesh->nlist[i][j];
      }
    }
  }

  // vertex coordinate list
  PyObject *vtx = PyArray_SimpleNewFromData(2, node_dims, NPY_FLOAT64, vtx_data);

  for (j = k = 0; j < nlen; j++) {
    bndidx[j] = (mesh->nlist[j].isBnd() ? k++ : -1);
  }

  // boundary element index list
  // note: this currently assumes that all elements contain the
  // same number of vertices!

  int nnd = 0, nface, sd, nn, nd, bn, *bndellist, *bndsdlist;
  nface = mesh->BoundaryList(&bndellist, &bndsdlist);
  for (j = 0; j < nface; j++) {
    nnd = std::max(nnd, mesh->elist[bndellist[j]]->nSideNode(bndsdlist[j]));
  }
  npy_intp face_dims[2] = {nface, nnd};

  nint *id, *idx_data = new nint[nface * nnd];
  for (j = 0, id = idx_data; j < nface; j++) {
    Element *pel = mesh->elist[bndellist[j]];
    sd = bndsdlist[j];
    nn = pel->nSideNode(sd);
    for (i = 0; i < nnd; i++) {
      if (i < nn) {
        nd = pel->Node[pel->SideNode(sd, i)];
        bn = bndidx[nd];
      } else {
        bn = -1;
      }
      *id++ = bn;
    }
  }
  PyObject *face = PyArray_SimpleNewFromData(2, face_dims, TOAST_NPY_INT, idx_data);

  // generate nodal permutation index list
  nint *p, *perm_data = new nint[nbnd];
  for (i = 0, p = perm_data; i < nlen; i++) {
    if (bndidx[i] >= 0) {
      *p++ = i;
    }
  }
  npy_intp perm_dim = nbnd;
  PyObject *perm = PyArray_SimpleNewFromData(1, &perm_dim, TOAST_NPY_INT, perm_data);

  // cleanup
  delete[] bndidx;

  return Py_BuildValue("OOO", vtx, face, perm);
}

// ===========================================================================

static PyObject *toast_mesh_node_count(PyObject *self, PyObject *args) {
  int hmesh;
  Mesh *mesh;
  if (!PyArg_ParseTuple(args, "i", &hmesh)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);
  return Py_BuildValue("i", mesh->nlen());
}

// ===========================================================================

static PyObject *toast_mesh_element_count(PyObject *self, PyObject *args) {
  int hmesh;
  Mesh *mesh;

  if (!PyArg_ParseTuple(args, "i", &hmesh)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);
  return Py_BuildValue("i", mesh->elen());
}

// ===========================================================================

static PyObject *toast_mesh_dim(PyObject *self, PyObject *args) {
  int hmesh;
  Mesh *mesh;

  if (!PyArg_ParseTuple(args, "i", &hmesh)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);
  return Py_BuildValue("i", mesh->Dimension());
}

// ===========================================================================

static PyObject *toast_mesh_bb(PyObject *self, PyObject *args) {
  int hmesh;
  Mesh *mesh;

  if (!PyArg_ParseTuple(args, "i", &hmesh)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  int i, dim = mesh->Dimension();
  Point pmin(dim), pmax(dim);
  mesh->BoundingBox(pmin, pmax);

  double *bb_data = new double[dim * 2];
  npy_intp bb_dims[2] = {dim, 2};

  for (i = 0; i < dim; i++) {
    bb_data[i * 2] = pmin[i];
    bb_data[i * 2 + 1] = pmax[i];
  }
  PyObject *bb = PyArray_SimpleNewFromData(2, bb_dims, NPY_FLOAT64, bb_data);

  return Py_BuildValue("O", bb);
}

// ===========================================================================

static PyObject *toast_basis_pos(PyObject *self, PyObject *args) {
  int hraster;
  Raster *raster;
  RDenseMatrix pos;

  if (!PyArg_ParseTuple(args, "i", &hraster)) {
    return NULL;
  }

  GETRASTER(raster, hraster);

  raster->BasisVoxelPositions(pos);
  npy_intp pos_dims[2] = {pos.nRows(), pos.nCols()};

  PyObject *pypos = PyArray_SimpleNew(2, pos_dims, NPY_FLOAT64);
  double *data = (double *)PyArray_DATA((PyArrayObject *)pypos);
  memcpy(data, pos.ValPtr(), pos.nRows() * pos.nCols() * sizeof(double));

  return Py_BuildValue("N", pypos);
}

// ===========================================================================

static PyObject *toast_sol_pos(PyObject *self, PyObject *args) {
  int hraster;
  Raster *raster;
  RDenseMatrix pos;

  if (!PyArg_ParseTuple(args, "i", &hraster)) {
    return NULL;
  }
  GETRASTER(raster, hraster);

  raster->SolutionVoxelPositions(pos);
  npy_intp pos_dims[2] = {pos.nRows(), pos.nCols()};

  PyObject *pypos = PyArray_SimpleNew(2, pos_dims, NPY_FLOAT64);
  double *data = (double *)PyArray_DATA((PyArrayObject *)pypos);
  memcpy(data, pos.ValPtr(), pos.nRows() * pos.nCols() * sizeof(double));

  return Py_BuildValue("N", pypos);
}

// ===========================================================================

static PyObject *toast_raster_glen(PyObject *self, PyObject *args) {
  int hraster;
  Raster *raster;

  if (!PyArg_ParseTuple(args, "i", &hraster)) {
    return NULL;
  }

  GETRASTER(raster, hraster);
  return Py_BuildValue("i", raster->GLen());
}

// ===========================================================================

static PyObject *toast_raster_blen(PyObject *self, PyObject *args) {
  int hraster;
  Raster *raster;

  if (!PyArg_ParseTuple(args, "i", &hraster)) {
    return NULL;
  }

  GETRASTER(raster, hraster);
  return Py_BuildValue("i", raster->BLen());
}

// ===========================================================================

static PyObject *toast_raster_slen(PyObject *self, PyObject *args) {
  int hraster;
  Raster *raster;

  if (!PyArg_ParseTuple(args, "i", &hraster)) {
    return NULL;
  }

  GETRASTER(raster, hraster);
  return Py_BuildValue("i", raster->SLen());
}

// ===========================================================================

static PyObject *toast_raster_matrix(PyObject *self, PyObject *args) {
  int i, hraster;
  Raster *raster;
  const char *mapstr;
  char srcid, tgtid;
  const RCompRowMatrix *matrix;

  if (!PyArg_ParseTuple(args, "is", &hraster, &mapstr)) {
    return NULL;
  }

  GETRASTER(raster, hraster);
  if (strlen(mapstr) != 4) {
    std::cerr << "toast.BasisMatrix: mapping string not recognised" << std::endl;
    return NULL;
  }
  if (!strncmp(mapstr + 1, "->", 2)) {
    srcid = (mapstr[0]);
    tgtid = (mapstr[3]);
  } else if (!strncmp(mapstr + 1, "<-", 2)) {
    srcid = (mapstr[3]);
    tgtid = (mapstr[0]);
  } else {
    std::cerr << "toast.BasisMatrix: mapping string not recognised" << std::endl;
    return NULL;
  }

  switch (srcid) {
    case 'M':
      switch (tgtid) {
        case 'G':
          matrix = (const RCompRowMatrix *)&raster->Mesh2GridMatrix();  // .B
          break;
        case 'B':
          matrix = (const RCompRowMatrix *)&raster->Mesh2BasisMatrix();  // .C
          break;
        default:
          std::cerr << "toast.BasisMatrix: source and target combination not recognised" << std::endl;
          return NULL;
      }
      break;
    case 'B':
      switch (tgtid) {
        case 'G':
          matrix = (const RCompRowMatrix *)&raster->Basis2GridMatrix();  // .GI
          break;
        case 'M':
          matrix = (const RCompRowMatrix *)&raster->Basis2MeshMatrix();  // .CI
          break;
        case 'S':
          matrix = (const RCompRowMatrix *)&raster->Basis2SolMatrix();  // .D
          break;
        default:
          std::cerr << "toast.BasisMatrix: source and target combination not recognised" << std::endl;
          return NULL;
      }
      break;
    case 'G':
      switch (tgtid) {
        case 'B':
          matrix = (const RCompRowMatrix *)&raster->Grid2BasisMatrix();  // .G
          break;
        case 'M':
          matrix = (const RCompRowMatrix *)&raster->Grid2MeshMatrix();  // .BI
          break;
        default:
          std::cerr << "toast.BasisMatrix: source and target combination not recognised" << std::endl;
          return NULL;
      }
      break;
    default:
      std::cerr << "toast.BasisMatrix: source not recognised" << std::endl;
      return NULL;
  }

  const idxtype *rowptr, *colidx;
  npy_intp nnz = matrix->GetSparseStructure(&rowptr, &colidx);
  const double *Mval = matrix->ValPtr();
  npy_intp nrp = matrix->nRows() + 1;

  // Allocate the numpy arrays for the CSR matrix
  PyObject *py_rp = PyArray_SimpleNew(1, &nrp, TOAST_NPY_INT);
  PyObject *py_ci = PyArray_SimpleNew(1, &nnz, TOAST_NPY_INT);
  PyObject *py_vl = PyArray_SimpleNew(1, &nnz, NPY_DOUBLE);

  // Copy the data over
  nint *rp = (nint *)PyArray_DATA((PyArrayObject *)py_rp);
  for (i = 0; i < nrp; i++) {
    rp[i] = rowptr[i];
  }
  nint *ci = (nint *)PyArray_DATA((PyArrayObject *)py_ci);
  for (i = 0; i < nnz; i++) {
    ci[i] = colidx[i];
  }
  double *val = (double *)PyArray_DATA((PyArrayObject *)py_vl);
  for (i = 0; i < nnz; i++) {
    val[i] = Mval[i];
  }

  return Py_BuildValue("NNN", py_rp, py_ci, py_vl);
}

// ===========================================================================

static PyObject *toast_raster_sol2basis(PyObject *self, PyObject *args) {
  int i, hraster;
  Raster *raster;

  if (!PyArg_ParseTuple(args, "i", &hraster)) {
    return NULL;
  }

  GETRASTER(raster, hraster);

  npy_intp slen = raster->SLen();
  PyObject *py_s2b = PyArray_SimpleNew(1, &slen, TOAST_NPY_INT);
  nint *sol2basis = (nint *)PyArray_DATA((PyArrayObject *)py_s2b);
  for (i = 0; i < raster->SLen(); i++) {
    sol2basis[i] = raster->Sol2Basis(i);  // raster->Sol2Basis(i) or raster->GetBasisIdx(i) ??
  }
  return Py_BuildValue("N", py_s2b);
}

// ===========================================================================

static PyObject *toast_raster_basis2sol(PyObject *self, PyObject *args) {
  int i, hraster;
  Raster *raster;

  if (!PyArg_ParseTuple(args, "i", &hraster)) {
    return NULL;
  }

  GETRASTER(raster, hraster);

  npy_intp blen = raster->BLen();
  PyObject *py_b2s = PyArray_SimpleNew(1, &blen, TOAST_NPY_INT);
  nint *basis2sol = (nint *)PyArray_DATA((PyArrayObject *)py_b2s);
  for (i = 0; i < raster->SLen(); i++) {
    basis2sol[i] = raster->Basis2Sol(i);
  }
  return Py_BuildValue("N", py_b2s);
}

// ===========================================================================

static PyObject *toast_make_mesh(PyObject *self, PyObject *args) {
  PyArrayObject *py_ndlist, *py_ellist, *py_eltp;

  if (!PyArg_ParseTuple(args, "O!O!O!", 
                              &PyArray_Type, &py_ndlist, 
                              &PyArray_Type, &py_ellist,
                              &PyArray_Type, &py_eltp)) {
    return NULL;
  }

  if (!AssertArray(py_ndlist, 2, NPY_DOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "vertices must be a two-dimensional array of doubles");
    return NULL;
  }
  npy_intp *node_dims = PyArray_DIMS(py_ndlist);
  int nvtx = node_dims[0];
  int dim = node_dims[1];
  double *vtx = (double *)PyArray_DATA(py_ndlist);

  if (!AssertArray(py_ellist, 2, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "elements must be a two-dimensional array of integers");
    return NULL;
  }
  npy_intp *ell_dims = PyArray_DIMS(py_ellist);
  int nel = ell_dims[0];
  int nnd0 = ell_dims[1];
  nint *idx = (nint *)PyArray_DATA(py_ellist);

  if (!AssertArrayDims(py_eltp, TOAST_NPY_INT, nel)) {
    PyErr_SetString(PyExc_ValueError, "element types must be a one-dimensional array of integers");
    return NULL;
  }
  nint *etp = (nint *)PyArray_DATA(py_eltp);


  QMMesh *mesh;
  BuildMesh(&mesh, nvtx, nel, dim, nnd0, vtx, idx, etp);

  int hmesh = g_meshmgr.Add(mesh);
  return Py_BuildValue("i", hmesh);
}

// ===========================================================================

static PyObject *toast_make_raster(PyObject *self, PyObject *args) {
  int hmesh, hraster;
  QMMesh *mesh;
  RDenseMatrix *bb = 0;
  PyArrayObject *py_size, *py_size_intm;

  if (!PyArg_ParseTuple(args, "iO!O!", &hmesh, 
                                       &PyArray_Type, &py_size, 
                                       &PyArray_Type, &py_size_intm)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);
  int mdim = mesh->Dimension();
  
  if (!AssertArrayDims(py_size, TOAST_NPY_INT, mdim)) {
    PyErr_SetString(PyExc_ValueError, "raster dimensions must be integeral and match mesh dimensionality");
    return NULL;
  }
  npy_intp *dims = PyArray_DIMS(py_size);
  int dim = dims[0];
  nint *size = (nint *)PyArray_DATA(py_size);
  IVector bdim(dim, size);

  if (!AssertArrayDims(py_size_intm, TOAST_NPY_INT, mdim)) {
    PyErr_SetString(PyExc_ValueError, "raster intermediate dimensions must be integeral and match mesh dimensionality");
    return NULL;
  }
  npy_intp *dims_intm = PyArray_DIMS(py_size_intm);
  int dim_intm = dims_intm[0];
  nint *size_intm = (nint *)PyArray_DATA(py_size_intm);
  IVector gdim(dim_intm, size_intm);

  Raster *raster;
  raster = new Raster_Pixel(bdim, gdim, mesh, bb);
  hraster = g_rastermgr.Add(raster);

  return Py_BuildValue("i", hraster);
}

// ===========================================================================

static PyObject *toast_clear_raster(PyObject *self, PyObject *args) {
  int hraster;

  if (!PyArg_ParseTuple(args, "i", &hraster)) {
    return NULL;
  }

  g_rastermgr.Delete(hraster);
  Py_RETURN_NONE;
}

// ===========================================================================

template <typename T>
void MapBasis(const Raster *raster, char srcid, char tgtid,
              const TVector<T> &src, TVector<T> &tgt) {
  switch (srcid) {
    case 'M':
      switch (tgtid) {
        case 'G':
          raster->Map_MeshToGrid(src, tgt);
          break;
        case 'B':
          raster->Map_MeshToBasis(src, tgt);
          break;
        case 'S':
          raster->Map_MeshToSol(src, tgt);
          break;
        default:
          std::cerr << "toast.MapBasis: target id not recognised" << std::endl;
      }
      break;
    case 'G':
      switch (tgtid) {
        case 'M':
          raster->Map_GridToMesh(src, tgt);
          break;
        case 'B':
          raster->Map_GridToBasis(src, tgt);
          break;
        case 'S':
          raster->Map_GridToSol(src, tgt);
          break;
        default:
          std::cerr << "toast.MapBasis: target id not recognised" << std::endl;
      }
      break;
    case 'B':
      switch (tgtid) {
        case 'G':
          raster->Map_BasisToGrid(src, tgt);
          break;
        case 'S':
          raster->Map_BasisToSol(src, tgt);
          break;
        case 'M':
          raster->Map_BasisToMesh(src, tgt);
          break;
        default:
          std::cerr << "toast.MapBasis: target id not recognised" << std::endl;
      }
      break;
    case 'S':
      switch (tgtid) {
        case 'B':
          raster->Map_SolToBasis(src, tgt);
          break;
        case 'G':
          raster->Map_SolToGrid(src, tgt);
          break;
        case 'M':
          raster->Map_SolToMesh(src, tgt);
          break;
        default:
          std::cerr << "toast.MapBasis: target id not recognised"
                    << std::endl;
      }
      break;
    default:
      std::cerr << "toast.MapBasis: source id not recognised" << std::endl;
  }
}

static PyObject *toast_map_basis(PyObject *self, PyObject *args) {
  int hraster;
  Raster *raster;
  const char *mapstr;
  char srcid, tgtid;
  PyArrayObject *py_srcvec;

  if (!PyArg_ParseTuple(args, "isO!", &hraster, &mapstr, 
                                      &PyArray_Type, &py_srcvec)) {
    return NULL;
  }

  GETRASTER(raster, hraster);

  if (strlen(mapstr) != 4) {
    std::cerr << "toast.MapBasis: mapping string not recognised" << std::endl;
    return NULL;
  }
  if (!strncmp(mapstr + 1, "->", 2)) {
    srcid = (mapstr[0]);
    tgtid = (mapstr[3]);
  } else if (!strncmp(mapstr + 1, "<-", 2)) {
    srcid = (mapstr[3]);
    tgtid = (mapstr[0]);
  } else {
    std::cerr << "toast.MapBasis: mapping string not recognised" << std::endl;
    return NULL;
  }

  npy_intp *dims = PyArray_DIMS(py_srcvec);
  int dtype = PyArray_TYPE(py_srcvec);
  int nsrc = dims[0];
  int ntgt;
  switch (tgtid) {
    case 'M':
      ntgt = raster->mesh().nlen();
      break;
    case 'G':
      ntgt = raster->GLen();
      break;
    case 'B':
      ntgt = raster->BLen();
      break;
    case 'S':
      ntgt = raster->SLen();
      break;
    default:
      std::cerr << "toast.MapBasis: target id not recognised" << std::endl;
      return NULL;
  }
  npy_intp py_ntgt = ntgt;
  PyObject *py_tgtvec = PyArray_SimpleNew(1, &py_ntgt, dtype);

  switch (dtype) {
    case NPY_DOUBLE: {

      if (!AssertArray(py_srcvec, 1, NPY_DOUBLE)) {
        PyErr_SetString(PyExc_ValueError, "data to be mapped must be one-dimensional, double or complex double");
        return NULL;
      }

      double *src_data = (double *)PyArray_DATA(py_srcvec);
      RVector tgt, src(nsrc, src_data);
      MapBasis(raster, srcid, tgtid, src, tgt);
      double *tgt_data = (double *)PyArray_DATA((PyArrayObject *) py_tgtvec);
      memcpy(tgt_data, tgt.data_buffer(), ntgt * sizeof(double));
    } break;
    case NPY_CDOUBLE: {

      if (!AssertArray(py_srcvec, 1, NPY_CDOUBLE)) {
        PyErr_SetString(PyExc_ValueError, "data to be mapped must be one-dimensional, double or complex double");
        return NULL;
      }

      std::complex<double> *src_data = (std::complex<double> *)PyArray_DATA(py_srcvec);
      CVector tgt, src(nsrc, src_data);
      MapBasis(raster, srcid, tgtid, src, tgt);
      std::complex<double> *tgt_data = (std::complex<double> *)PyArray_DATA((PyArrayObject *) py_tgtvec);
      memcpy(tgt_data, tgt.data_buffer(), ntgt * sizeof(std::complex<double>));
    } break;
#ifdef UNDEF
    case NPY_FLOAT: {
      float *src_data = (float *)PyArray_DATA(py_srcvec);
      FVector tgt, src(nsrc, src_data);
      MapBasis(raster, srcid, tgtid, src, tgt);
      float *tgt_data = (float *)PyArray_DATA(py_tgtvec);
      memcpy(tgt_data, tgt.data_buffer(), ntgt * sizeof(float));
    } break;
    case NPY_CFLOAT: {
      scomplex *src_data = (scomplex *)PyArray_DATA(py_srcvec);
      SCVector tgt, src(nsrc, src_data);
      MapBasis(raster, srcid, tgtid, src, tgt);
      scomplex *tgt_data = (scomplex *)PyArray_DATA(py_tgtvec);
      memcpy(tgt_data, tgt.data_buffer(), ntgt * sizeof(scomplex));
    } break;
#endif
    default:
      std::cerr << "toast.MapBasis: vector type not recognised" << std::endl;
      return NULL;
  }

  return Py_BuildValue("O", py_tgtvec);
}

// ===========================================================================

static PyObject *toast_readqm(PyObject *self, PyObject *args) {
  const char *qmname;
  int hmesh;
  QMMesh *mesh;

  if (!PyArg_ParseTuple(args, "is", &hmesh, &qmname)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  std::ifstream ifs(qmname);
  if(!ifs) {
    PyErr_SetString(PyExc_ValueError, "QM file could not be opened");
    return NULL;
  }
  mesh->LoadQM(ifs);

  std::cout << "QM: " << mesh->nQ << " sources, " << mesh->nM
            << " detectors, " << mesh->nQM << " measurements" << std::endl;

  Py_RETURN_NONE;
}

// ===========================================================================

static PyObject *toast_read_nim(PyObject *self, PyObject *args) {
  const char *nimname;
  char cbuf[256];
  int i, j = 0, idx, imgsize = 0;

  if (!PyArg_ParseTuple(args, "si", &nimname, &idx)) {
    return NULL;
  }

  std::ifstream ifs(nimname);
  if (!ifs.getline(cbuf, 256)) {
    return NULL;
  }
  if (strcmp(cbuf, "NIM") && strcmp(cbuf, "RIM")) {
    return NULL;
  }
  do {
    ifs.getline(cbuf, 256);
    if (!strncasecmp(cbuf, "ImageSize", 9))
      sscanf(cbuf + 11, "%d", &imgsize);
  } while (strcasecmp(cbuf, "EndHeader"));

  if (!imgsize) {
    return NULL;
  }

  double *img = new double[imgsize];
  for (;;) {
    do {
      ifs.getline(cbuf, 256);
    } while (ifs.good() && strncasecmp(cbuf, "Image", 5));

    if (!ifs.good()) {
      break;
    }

    for (i = 0; i < imgsize; i++) {
      ifs >> img[i];
    }

    if (j++ == idx) {
      break;
    }
  }

  npy_intp dims = imgsize;
  PyObject *nim = PyArray_SimpleNewFromData(1, &dims, NPY_DOUBLE, img);

  return Py_BuildValue("O", nim);
}

// ===========================================================================

static PyObject *toast_write_nim(PyObject *self, PyObject *args) {
  const char *nimname, *meshname;
  PyArrayObject *py_nim;

  if (!PyArg_ParseTuple(args, "ssO!", &nimname, &meshname, 
                                      &PyArray_Type, &py_nim)) {
    return NULL;
  }

  std::ofstream ofs(nimname);
  ofs << "NIM" << std::endl;
  ofs << "Mesh = " << meshname << std::endl;
  ofs << "SolutionType = N/A" << std::endl;

  npy_intp nd = PyArray_NDIM((PyArrayObject *)py_nim);
  npy_intp n = 1;
  for (int i = 0; i < nd; i++) {
    n *= PyArray_DIM((PyArrayObject *)py_nim, i);
  }
  ofs << "ImageSize = " << n << std::endl;

  ofs << "EndHeader" << std::endl;

  ofs << "Image 0" << std::endl;

  ofs.precision(12);
  ofs.setf(std::ios::scientific);
  for (int i = 0; i < n; i++) {
    ofs << (*(double *)PyArray_GETPTR1((PyArrayObject *)py_nim, i)) << ' ';
  }
  ofs << std::endl;

  Py_RETURN_NONE;
}

// ===========================================================================

static PyObject *toast_sysmat_cw(PyObject *self, PyObject *args) {
  bool elbasis = false;  // for now

  int i, hmesh;
  QMMesh *mesh;
  PyArrayObject *py_mua, *py_mus, *py_ref;

  if (!PyArg_ParseTuple(args, "iO!O!O!", &hmesh, 
                                         &PyArray_Type, &py_mua, 
                                         &PyArray_Type, &py_mus, 
                                         &PyArray_Type, &py_ref)) {
    return NULL;
  }
  GETMESH(mesh, hmesh);

  int nlen = mesh->nlen();

  if (!AssertArrayDims(py_mua, NPY_DOUBLE, nlen)) {
    PyErr_SetString(PyExc_ValueError, "mua must be a one-dimensional vector of doubles");
    return NULL;
  }
  if (!AssertArrayDims(py_mus, NPY_DOUBLE, nlen)) {
    PyErr_SetString(PyExc_ValueError, "mus must be a one-dimensional vector of doubles");
    return NULL;
  }
  if (!AssertArrayDims(py_ref, NPY_DOUBLE, nlen)) {
    PyErr_SetString(PyExc_ValueError, "ref must be a one-dimensional vector of doubles");
    return NULL;
  }
  double *mua = (double *)PyArray_DATA(py_mua);
  double *mus = (double *)PyArray_DATA(py_mus);
  double *ref = (double *)PyArray_DATA(py_ref);

  RVector prm(nlen);

  // Set optical coefficients
  Solution sol(OT_NPARAM, nlen);
  for (i = 0; i < nlen; i++) {
    prm[i] = mua[i] * c0 / ref[i];
  }
  sol.SetParam(OT_CMUA, prm);
  for (i = 0; i < nlen; i++) {
    prm[i] = c0 / (3.0 * ref[i] * (mua[i] + mus[i]));
  }
  sol.SetParam(OT_CKAPPA, prm);
  for (int i = 0; i < nlen; i++) {
    prm[i] = c0 / (2.0 * ref[i] * A_Keijzer(ref[i]));
  }
  sol.SetParam(OT_C2A, prm);

  // Create forward solver to initialise system matrix
  RFwdSolver FWS(mesh, LSOLVER_ITERATIVE, 1e-10);
  FWS.SetDataScaling(DATA_LOG);

  FWS.Allocate();
  FWS.AssembleSystemMatrix(sol, 0, elbasis);

  const idxtype *rowptr, *colidx;
  npy_intp nnz = FWS.F->GetSparseStructure(&rowptr, &colidx);
  const double *Fval = FWS.F->ValPtr();
  npy_intp nrp = nlen + 1;

  // Allocate the numpy arrays for the CSR matrix
  PyObject *py_rp = PyArray_SimpleNew(1, &nrp, TOAST_NPY_INT);
  PyObject *py_ci = PyArray_SimpleNew(1, &nnz, TOAST_NPY_INT);
  PyObject *py_vl = PyArray_SimpleNew(1, &nnz, NPY_DOUBLE);

  // Copy the data over
  nint *rp = (nint *)PyArray_DATA((PyArrayObject *)py_rp);
  for (i = 0; i < nrp; i++) {
    rp[i] = (int) rowptr[i];
  }
  nint *ci = (nint *)PyArray_DATA((PyArrayObject *)py_ci);
  for (i = 0; i < nnz; i++) {
    ci[i] = (int) colidx[i];
  }
  double *val = (double *)PyArray_DATA((PyArrayObject *)py_vl);
  for (i = 0; i < nnz; i++) {
    val[i] = Fval[i];
  }

  return Py_BuildValue("NNN", py_rp, py_ci, py_vl);
}

// ===========================================================================

static PyObject *toast_sysmat(PyObject *self, PyObject *args) {
  bool elbasis = false;  // for now

  int i, hmesh;
  double freq;
  QMMesh *mesh;
  PyArrayObject *py_mua, *py_mus, *py_ref;

  if (!PyArg_ParseTuple(args, "iO!O!O!d", &hmesh, 
                                          &PyArray_Type, &py_mua, 
                                          &PyArray_Type, &py_mus, 
                                          &PyArray_Type, &py_ref, 
                                          &freq)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  double omega = freq * 2.0 * Pi * 1e-6;
  int nlen = mesh->nlen();

  if (!AssertArrayDims(py_mua, NPY_DOUBLE, nlen)) {
    PyErr_SetString(PyExc_ValueError, "mua must be a one-dimensional vector of doubles");
    return NULL;
  }
  if (!AssertArrayDims(py_mus, NPY_DOUBLE, nlen)) {
    PyErr_SetString(PyExc_ValueError, "mus must be a one-dimensional vector of doubles");
    return NULL;
  }
  if (!AssertArrayDims(py_ref, NPY_DOUBLE, nlen)) {
    PyErr_SetString(PyExc_ValueError, "ref must be a one-dimensional vector of doubles");
    return NULL;
  }
  double *mua = (double *)PyArray_DATA(py_mua);
  double *mus = (double *)PyArray_DATA(py_mus);
  double *ref = (double *)PyArray_DATA(py_ref);

  RVector prm(nlen);

  // Set optical coefficients
  Solution sol(OT_NPARAM, nlen);
  for (i = 0; i < nlen; i++) {
    prm[i] = mua[i] * c0 / ref[i];
  }
  sol.SetParam(OT_CMUA, prm);
  for (i = 0; i < nlen; i++) {
    prm[i] = c0 / (3.0 * ref[i] * (mua[i] + mus[i]));
  }
  sol.SetParam(OT_CKAPPA, prm);
  for (int i = 0; i < nlen; i++) {
    prm[i] = c0 / (2.0 * ref[i] * A_Keijzer(ref[i]));
  }
  sol.SetParam(OT_C2A, prm);

  // Create forward solver to initialise system matrix
  CFwdSolver FWS(mesh, LSOLVER_ITERATIVE, 1e-10);
  FWS.SetDataScaling(DATA_LOG);

  FWS.Allocate();
  FWS.AssembleSystemMatrix(sol, omega, elbasis);

  const idxtype *rowptr, *colidx;
  npy_intp nnz = FWS.F->GetSparseStructure(&rowptr, &colidx);
  const std::complex<double> *Fval = FWS.F->ValPtr();
  npy_intp nrp = nlen + 1;

  // Allocate the numpy arrays for the CSR matrix
  PyObject *py_rp = PyArray_SimpleNew(1, &nrp, TOAST_NPY_INT);
  PyObject *py_ci = PyArray_SimpleNew(1, &nnz, TOAST_NPY_INT);
  PyObject *py_vl = PyArray_SimpleNew(1, &nnz, NPY_CDOUBLE);

  // Copy the data over
  nint *rp = (nint *)PyArray_DATA((PyArrayObject *)py_rp);
  for (i = 0; i < nrp; i++) {
    rp[i] = rowptr[i];
  }
  nint *ci = (nint *)PyArray_DATA((PyArrayObject *)py_ci);
  for (i = 0; i < nnz; i++) {
    ci[i] = colidx[i];
  }
  std::complex<double> *val = (std::complex<double> *)PyArray_DATA((PyArrayObject *)py_vl);
  for (i = 0; i < nnz; i++) {
    val[i] = Fval[i];
  }

  return Py_BuildValue("NNN", py_rp, py_ci, py_vl);
}

// ===========================================================================

void CalcQvec(const QMMesh *mesh, SourceMode qtype,
              SRC_PROFILE qprof, double qwidth, CCompRowMatrix *qvec) {
  int i, n, nQ;

  n = mesh->nlen();
  nQ = mesh->nQ;

  // build the source vectors
  qvec->New(nQ, n);

  for (i = 0; i < nQ; i++) {
    CVector q(n);
    switch (qprof) {
      case PROF_POINT:
        SetReal(q, QVec_Point(*mesh, mesh->Q[i], qtype));
        break;
      case PROF_GAUSSIAN:
        SetReal(q, QVec_Gaussian(*mesh, mesh->Q[i], qwidth, qtype));
        break;
      case PROF_COSINE:
        SetReal(q, QVec_Cosine(*mesh, mesh->Q[i], qwidth, qtype));
        break;
      case PROF_COMPLETETRIG:
        std::cerr << "Not implemented" << std::endl;
        // q = CompleteTrigSourceVector (*mesh, i);
        break;
    }
    qvec->SetRow(i, q);
  }
}

static PyObject *toast_qvec(PyObject *self, PyObject *args, PyObject *keywds) {
  int i, hmesh;
  const char *typestr = "Neumann";
  const char *profstr = "Gaussian";
  double qwidth = 1.0;
  QMMesh *mesh;

  static const char *kwlist[] = {"mesh", "type", "shape", "width", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "i|ssd", (char **)kwlist, &hmesh, &typestr, &profstr, &qwidth)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  SourceMode qtype = SRCMODE_NEUMANN;
  if (!strcasecmp(typestr, "Neumann")) {
    qtype = SRCMODE_NEUMANN;
  } else if (!strcasecmp(typestr, "Isotropic")) {
    qtype = SRCMODE_ISOTROPIC;
  } else {
    std::cerr << "toast.Qvec: Invalid source type" << std::endl;
  }

  SRC_PROFILE qprof = PROF_POINT;
  if (!strcasecmp(profstr, "Point")) {
    qprof = PROF_POINT;
  } else if (!strcasecmp(profstr, "Gaussian")) {
    qprof = PROF_GAUSSIAN;
  } else if (!strcasecmp(profstr, "Cosine")) {
    qprof = PROF_COSINE;
  } else if (!strcasecmp(profstr, "TrigBasis")) {
    qprof = PROF_COMPLETETRIG;
  } else {
    std::cerr << "toast.Qvec: Invalid source profile" << std::endl;
  }

  if (qprof != PROF_POINT) {
    if (qwidth <= 0.0) {
      std::cerr << "toast.Qvec: Invalid source width" << std::endl;
    }
  }

  CCompRowMatrix qvec;
  CalcQvec(mesh, qtype, qprof, qwidth, &qvec);

  const idxtype *rowptr, *colidx;
  npy_intp nnz = qvec.GetSparseStructure(&rowptr, &colidx);
  const std::complex<double> *qval = qvec.ValPtr();
  npy_intp nrp = qvec.nRows() + 1;

  // Allocate the numpy arrays for the CSR matrix
  PyObject *py_rp = PyArray_SimpleNew(1, &nrp, TOAST_NPY_INT);
  PyObject *py_ci = PyArray_SimpleNew(1, &nnz, TOAST_NPY_INT);
  PyObject *py_vl = PyArray_SimpleNew(1, &nnz, NPY_CDOUBLE);

  // Copy the data over
  nint *rp = (nint *)PyArray_DATA((PyArrayObject *)py_rp);
  for (i = 0; i < nrp; i++) {
    rp[i] = rowptr[i];
  }
  nint *ci = (nint *)PyArray_DATA((PyArrayObject *)py_ci);
  for (i = 0; i < nnz; i++) {
    ci[i] = colidx[i];
  }
  std::complex<double> *val = (std::complex<double> *)PyArray_DATA((PyArrayObject *)py_vl);
  for (i = 0; i < nnz; i++) {
    val[i] = qval[i];
  }

  return Py_BuildValue("NNN", py_rp, py_ci, py_vl);
}

// ===========================================================================

void CalcMvec(const QMMesh *mesh, SRC_PROFILE mprof, double mwidth,
              RVector *ref, CCompRowMatrix *mvec) {
  int n, nM;
  int i, j;

  n = mesh->nlen();
  nM = mesh->nM;

  // build the measurement vectors
  mvec->New(nM, n);
  for (i = 0; i < nM; i++) {
    CVector m(n);
    switch (mprof) {
      case PROF_GAUSSIAN:
        SetReal(m, QVec_Gaussian(*mesh, mesh->M[i], mwidth, SRCMODE_NEUMANN));
        break;
      case PROF_COSINE:
        SetReal(m, QVec_Cosine(*mesh, mesh->M[i], mwidth, SRCMODE_NEUMANN));
        break;
      case PROF_COMPLETETRIG:
        std::cerr << "Not implemented" << std::endl;
        // m = CompleteTrigSourceVector (*mesh, i);
        break;
      default:
        std::cerr << "Not implemented" << std::endl;
        break;
    }
    for (j = 0; j < n; j++) {
      m[j] *= c0 / (2.0 * (*ref)[j] * A_Keijzer((*ref)[j]));
    }
    mvec->SetRow(i, m);
  }
}

static PyObject *toast_mvec(PyObject *self, PyObject *args, PyObject *keywds) {
  int i, hmesh;
  const char *profstr = "Gaussian";
  double mwidth = 1.0;
  QMMesh *mesh;
  double refind = 1.0;

  static const char *kwlist[] = {"mesh", "shape", "width", "ref", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "isdd", (char **)kwlist, &hmesh, &profstr, &mwidth, &refind)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  SRC_PROFILE mprof = PROF_POINT;
  if (!strcasecmp(profstr, "Point")) {
    mprof = PROF_POINT;
  } else if (!strcasecmp(profstr, "Gaussian")) {
    mprof = PROF_GAUSSIAN;
  } else if (!strcasecmp(profstr, "Cosine")) {
    mprof = PROF_COSINE;
  } else if (!strcasecmp(profstr, "TrigBasis")) {
    mprof = PROF_COMPLETETRIG;
  } else {
    std::cerr << "toast.Qvec: Invalid source profile" << std::endl;
  }

#ifdef UNDEF
  npy_intp *ref_dims = PyArray_DIMS((PyArrayObject *)py_ref);
  int reflen = ref_dims[0] * ref_dims[1];
  RVector ref(mesh->nlen());
  if (reflen == 1) {
    ref = *(double *)PyArray_DATA((PyArrayObject *)py_ref);
  } else {
    if (reflen != mesh->nlen())
      return NULL;
    memcpy(ref.data_buffer(), PyArray_DATA((PyArrayObject *)py_ref), reflen * sizeof(double));
  }
#endif
  RVector ref(mesh->nlen());
  ref = refind;

  if (mprof != PROF_POINT) {
    if (mwidth <= 0.0) {
      std::cerr << "toast.Mvec: Invalid detector width" << std::endl;
    }
  }

  CCompRowMatrix mvec;
  CalcMvec(mesh, mprof, mwidth, &ref, &mvec);

  const idxtype *rowptr, *colidx;
  npy_intp nnz = mvec.GetSparseStructure(&rowptr, &colidx);
  const std::complex<double> *mval = mvec.ValPtr();
  npy_intp nrp = mvec.nRows() + 1;

  // Allocate the numpy arrays for the CSR matrix
  PyObject *py_rp = PyArray_SimpleNew(1, &nrp, TOAST_NPY_INT);
  PyObject *py_ci = PyArray_SimpleNew(1, &nnz, TOAST_NPY_INT);
  PyObject *py_vl = PyArray_SimpleNew(1, &nnz, NPY_CDOUBLE);

  // Copy the data over
  nint *rp = (nint *)PyArray_DATA((PyArrayObject *)py_rp);
  for (i = 0; i < nrp; i++) {
    rp[i] = rowptr[i];
  }
  nint *ci = (nint *)PyArray_DATA((PyArrayObject *)py_ci);
  for (i = 0; i < nnz; i++) {
    ci[i] = colidx[i];
  }
  std::complex<double> *val = (std::complex<double> *)PyArray_DATA((PyArrayObject *)py_vl);
  for (i = 0; i < nnz; i++) {
    val[i] = mval[i];
  }

  return Py_BuildValue("NNN", py_rp, py_ci, py_vl);
}

// ===========================================================================

static PyObject *toast_fields(PyObject *self, PyObject *args) {
  int hmesh, hraster, nQ;
  double freq;
  QMMesh *mesh;
  Raster *raster;

  PyArrayObject *py_qvec_vl, *py_qvec_rp, *py_qvec_ci;
  PyArrayObject *py_mua, *py_mus, *py_ref;

  if (!PyArg_ParseTuple(args, "iiiO!O!O!O!O!O!d",
                        &hmesh, &hraster, &nQ, 
                        &PyArray_Type, &py_qvec_vl, 
                        &PyArray_Type, &py_qvec_rp, 
                        &PyArray_Type, &py_qvec_ci,
                        &PyArray_Type, &py_mua, 
                        &PyArray_Type, &py_mus, 
                        &PyArray_Type, &py_ref,
                        &freq)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);
  if (hraster < 0) {
    raster = NULL;
  } else {
    GETRASTER(raster, hraster);
  }

  int n = mesh->nlen();

  if (!AssertArray(py_qvec_rp, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "qvec row pointer must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_qvec_ci, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "qvec column indices must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_qvec_vl, 1, NPY_CDOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "qvec vals must be a one-dimensional vector of complex doubles");
    return NULL;
  }
  nint *qrowptr = (nint *)PyArray_DATA(py_qvec_rp);
  nint *qcolidx = (nint *)PyArray_DATA(py_qvec_ci);
  std::complex<double> *qval = (std::complex<double> *)PyArray_DATA(py_qvec_vl);
  CCompRowMatrix qvec(nQ, n, qrowptr, qcolidx, qval);

  if (!AssertArrayDims(py_mua, NPY_DOUBLE, n)) {
    PyErr_SetString(PyExc_ValueError, "mua must be a one-dimensional vector of doubles of length equal to mesh nodes");
    return NULL;
  }
  if (!AssertArrayDims(py_mus, NPY_DOUBLE, n)) {
    PyErr_SetString(PyExc_ValueError, "mus must be a one-dimensional vector of doubles of length equal to mesh nodes");
    return NULL;
  }
  if (!AssertArrayDims(py_ref, NPY_DOUBLE, n)) {
    PyErr_SetString(PyExc_ValueError, "ref must be a one-dimensional vector of doubles of length equal to mesh nodes");
    return NULL;
  }
  double *mua_ptr = (double *)PyArray_DATA(py_mua);
  RVector mua(n, mua_ptr, SHALLOW_COPY);

  double *mus_ptr = (double *)PyArray_DATA(py_mus);
  RVector mus(n, mus_ptr, SHALLOW_COPY);

  double *ref_ptr = (double *)PyArray_DATA(py_ref);
  RVector ref(n, ref_ptr, SHALLOW_COPY);

  int slen = (raster ? raster->SLen() : n);

  CVector *dphi;
  CVector sphi(slen);
 
  CalcFields(mesh, raster, qvec, mua, mus, ref, freq, "DIRECT", 1e-12, &dphi);

  // Map or copy to output
  npy_intp dmx_dims[2] = {slen, nQ};
  PyObject *dmx = PyArray_SimpleNew(2, dmx_dims, NPY_CDOUBLE);
  std::complex<double> *dmx_ptr = (std::complex<double> *)PyArray_DATA((PyArrayObject *)dmx);
  for (int i = 0; i < nQ; i++) {
    std::complex<double> *dp = dmx_ptr + i;
    if (raster) {
      raster->Map_MeshToSol(dphi[i], sphi);
    } else {
      sphi = dphi[i];
    }
    for (int j = 0; j < slen; j++) {
      *dp = sphi[j];
      dp += nQ;
    }
  }
  delete[] dphi;

  PyObject *ret = Py_BuildValue("O", dmx);
  Py_DECREF(dmx);

  if (ret) {
    return ret;
  } else {
    Py_RETURN_NONE;
  }
}

// ===========================================================================

// Calculate Jacobian from given direct and adjoint fields and boundary
// projection data

static PyObject *toast_jacobian(PyObject *self, PyObject *args) {
  int hmesh, hraster, i, j;
  QMMesh *mesh;
  Raster *raster;
  PyArrayObject *py_dphi, *py_aphi, *py_proj;

  if (!PyArg_ParseTuple(args, "iiO!O!O!", &hmesh, &hraster, 
                                          &PyArray_Type, &py_dphi,
                                          &PyArray_Type, &py_aphi, 
                                          &PyArray_Type, &py_proj)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);
  GETRASTER(raster, hraster);

  int n = mesh->nlen();
  int nq = mesh->nQ;
  int nm = mesh->nM;
  int nqm = mesh->nQM;

  // set up direct fields
  if (!AssertArrayDims(py_dphi, NPY_CDOUBLE, n, nq)) {
    PyErr_SetString(PyExc_ValueError, "dphi must be a two-dimensional array of dimension nodes x sources");
    return NULL;
  }

  std::complex<double> *dphi_ptr = (std::complex<double> *)PyArray_DATA(py_dphi);
  CVector *dphi = new CVector[nq];
  for (i = 0; i < nq; i++) {
    dphi[i].New(n);
    for (j = 0; j < n; j++) {
      dphi[i][j] = dphi_ptr[i + j * nq];
    }
  }

  // set up adjoint fields
  if (!AssertArrayDims(py_aphi, NPY_CDOUBLE, n, nq)) {
    PyErr_SetString(PyExc_ValueError, "aphi must be a two-dimensional array of dimension nodes x detectors");
    return NULL;
  }

  std::complex<double> *aphi_ptr = (std::complex<double> *)PyArray_DATA(py_aphi);
  CVector *aphi = new CVector[nm];
  for (i = 0; i < nm; i++) {
    aphi[i].New(n);
    for (j = 0; j < n; j++) {
      aphi[i][j] = aphi_ptr[i + j * nm];
    }
  }

  // copy projections
  if (!AssertArrayDims(py_proj, NPY_CDOUBLE, nqm)) {
    PyErr_SetString(PyExc_ValueError, "projection data must be a one-dimensional array dimension sources x detectors");
    return NULL;
  }

  std::complex<double> *proj_ptr = (std::complex<double> *)PyArray_DATA(py_proj);
  CVector proj(nqm, proj_ptr, SHALLOW_COPY);

  RDenseMatrix J;
  CalcJacobian(mesh, raster, dphi, aphi, &proj, DATA_LOG, J);

  npy_intp J_dims[2] = {J.nRows(), J.nCols()};
  PyObject *pyJ = PyArray_SimpleNew(2, J_dims, NPY_DOUBLE);
  double *pyJ_data = (double *)PyArray_DATA((PyArrayObject *)pyJ);
  const double *J_data = J.ValPtr();
  memcpy(pyJ_data, J_data, J.nRows() * J.nCols() * sizeof(double));

  // cleanup
  delete[] dphi;
  delete[] aphi;

  PyObject *res = Py_BuildValue("O", pyJ);
  Py_DECREF(pyJ);
  return res;
}


static PyObject *toast_jacobian_optical(PyObject *self, PyObject *args) {
  int hmesh, hraster;
  QMMesh *mesh;
  Raster *raster;
  char *solver;
  double freq, tol;
  PyArrayObject *py_qvec_vl, *py_qvec_rp, *py_qvec_ci;
  PyArrayObject *py_mvec_vl, *py_mvec_rp, *py_mvec_ci;
  PyArrayObject *py_mua, *py_mus, *py_ref;

  if (!PyArg_ParseTuple(args, "iiO!O!O!O!O!O!O!O!O!dsd", &hmesh, &hraster,
                        &PyArray_Type, &py_qvec_vl,
                        &PyArray_Type, &py_qvec_rp,
                        &PyArray_Type, &py_qvec_ci,
                        &PyArray_Type, &py_mvec_vl,
                        &PyArray_Type, &py_mvec_rp,
                        &PyArray_Type, &py_mvec_ci,
                        &PyArray_Type, &py_mua,
                        &PyArray_Type, &py_mus,
                        &PyArray_Type, &py_ref,
                        &freq, &solver, &tol))
    return NULL;

  GETMESH(mesh, hmesh);
  GETRASTER(raster, hraster);

  int n = mesh->nlen();
  int nQ = mesh->nQ;
  int nM = mesh->nM;

  if (!AssertArray(py_qvec_rp, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "qvec row pointer must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_qvec_ci, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "qvec column indices must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_qvec_vl, 1, NPY_CDOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "qvec vals must be a one-dimensional vector of complex doubles");
    return NULL;
  }
  nint *qrowptr = (nint *)PyArray_DATA(py_qvec_rp);
  nint *qcolidx = (nint *)PyArray_DATA(py_qvec_ci);
  std::complex<double> *qval = (std::complex<double> *)PyArray_DATA(py_qvec_vl);
  const CCompRowMatrix qvec(nQ, n, qrowptr, qcolidx, qval /*, SHALLOW_COPY*/);

  if (!AssertArray(py_mvec_rp, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "mvec row pointer must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_mvec_ci, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "mvec column indices must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_mvec_vl, 1, NPY_CDOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "mvec vals must be a one-dimensional vector of complex doubles");
    return NULL;
  }
  nint *mrowptr = (nint *)PyArray_DATA(py_mvec_rp);
  nint *mcolidx = (nint *)PyArray_DATA(py_mvec_ci);
  std::complex<double> *mval = (std::complex<double> *)PyArray_DATA(py_mvec_vl);
  const CCompRowMatrix mvec(nM, n, mrowptr, mcolidx, mval /*, SHALLOW_COPY*/);

  if (!AssertArrayDims(py_mua, NPY_DOUBLE, n)) {
    PyErr_SetString(PyExc_ValueError, "mua must be a one-dimensional vector of doubles of length equal to mesh nodes");
    return NULL;
  }
  if (!AssertArrayDims(py_mus, NPY_DOUBLE, n)) {
    PyErr_SetString(PyExc_ValueError, "mua must be a one-dimensional vector of doubles of length equal to mesh nodes");
    return NULL;
  }
  if (!AssertArrayDims(py_ref, NPY_DOUBLE, n)) {
    PyErr_SetString(PyExc_ValueError, "mua must be a one-dimensional vector of doubles of length equal to mesh nodes");
    return NULL;
  }
  double *mua_ptr = (double *)PyArray_DATA(py_mua);
  const RVector mua(n, mua_ptr, SHALLOW_COPY);

  double *mus_ptr = (double *)PyArray_DATA(py_mus);
  const RVector mus(n, mus_ptr, SHALLOW_COPY);

  double *ref_ptr = (double *)PyArray_DATA(py_ref);
  const RVector ref(n, ref_ptr, SHALLOW_COPY);

  RDenseMatrix J;
  CalcJacobian(mesh, raster, qvec, mvec, mua, mus, ref, freq, solver, tol, J);
    
  npy_intp J_dims[2] = {J.nRows(), J.nCols()};
  PyObject *pyJ = PyArray_SimpleNew(2, J_dims, NPY_DOUBLE);
  double *pyJ_data = (double *)PyArray_DATA((PyArrayObject *)pyJ);
  const double *J_data = J.ValPtr();
  memcpy(pyJ_data, J_data, J.nRows() * J.nCols() * sizeof(double));
  
  PyObject *res = Py_BuildValue("O", pyJ);
  Py_DECREF(pyJ);
  return res;
}


static PyObject *toast_jacobianCW(PyObject *self, PyObject *args) {
  int hmesh, hraster;
  QMMesh *mesh;
  Raster *raster;
  char *solver;
  double tol;
  PyArrayObject *py_qvec_vl, *py_qvec_rp, *py_qvec_ci;
  PyArrayObject *py_mvec_vl, *py_mvec_rp, *py_mvec_ci;
  PyArrayObject *py_mua, *py_mus, *py_ref;

  if (!PyArg_ParseTuple(args, "iiO!O!O!O!O!O!O!O!O!sd", &hmesh, &hraster,
                        &PyArray_Type, &py_qvec_vl,
                        &PyArray_Type, &py_qvec_rp,
                        &PyArray_Type, &py_qvec_ci,
                        &PyArray_Type, &py_mvec_vl,
                        &PyArray_Type, &py_mvec_rp,
                        &PyArray_Type, &py_mvec_ci,
                        &PyArray_Type, &py_mua,
                        &PyArray_Type, &py_mus,
                        &PyArray_Type, &py_ref,
                        &solver, &tol)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);
  GETRASTER(raster, hraster);

  int n = mesh->nlen();
  int nQ = mesh->nQ;
  int nM = mesh->nM;

  if (!AssertArray(py_qvec_rp, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "qvec row pointer must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_qvec_ci, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "qvec column indices must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_qvec_vl, 1, NPY_DOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "qvec vals must be a one-dimensional vector of doubles");
    return NULL;
  }
  nint *qrowptr = (nint *)PyArray_DATA(py_qvec_rp);
  nint *qcolidx = (nint *)PyArray_DATA(py_qvec_ci);
  double *qval = (double *)PyArray_DATA(py_qvec_vl);
  const RCompRowMatrix qvec(nQ, n, qrowptr, qcolidx, qval /*, SHALLOW_COPY*/);

  if (!AssertArray(py_mvec_rp, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "mvec row pointer must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_mvec_ci, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "mvec column indices must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_mvec_vl, 1, NPY_DOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "mvec vals must be a one-dimensional vector of doubles");
    return NULL;
  }
  nint *mrowptr = (nint *)PyArray_DATA(py_mvec_rp);
  nint *mcolidx = (nint *)PyArray_DATA(py_mvec_ci);
  double *mval = (double *)PyArray_DATA(py_mvec_vl);
  const RCompRowMatrix mvec(nM, n, mrowptr, mcolidx, mval /*, SHALLOW_COPY*/);

  if (!AssertArrayDims(py_mua, NPY_DOUBLE, n)) {
    PyErr_SetString(PyExc_ValueError, "mua must be a one-dimensional vector of doubles of length equal to mesh nodes");
    return NULL;
  }
  if (!AssertArrayDims(py_mus, NPY_DOUBLE, n)) {
    PyErr_SetString(PyExc_ValueError, "mua must be a one-dimensional vector of doubles of length equal to mesh nodes");
    return NULL;
  }
  if (!AssertArrayDims(py_ref, NPY_DOUBLE, n)) {
    PyErr_SetString(PyExc_ValueError, "mua must be a one-dimensional vector of doubles of length equal to mesh nodes");
    return NULL;
  }
  double *mua_ptr = (double *)PyArray_DATA(py_mua);
  const RVector mua(n, mua_ptr, SHALLOW_COPY);

  double *mus_ptr = (double *)PyArray_DATA(py_mus);
  const RVector mus(n, mus_ptr, SHALLOW_COPY);

  double *ref_ptr = (double *)PyArray_DATA(py_ref);
  const RVector ref(n, ref_ptr, SHALLOW_COPY);

  RDenseMatrix J;
  CalcJacobianCW(mesh, raster, qvec, mvec, mua, mus, ref, solver, tol, J);

  npy_intp J_dims[2] = {J.nRows(), J.nCols()};
  PyObject *pyJ = PyArray_SimpleNew(2, J_dims, NPY_DOUBLE);
  double *pyJ_data = (double *)PyArray_DATA((PyArrayObject *)pyJ);
  const double *J_data = J.ValPtr();
  memcpy(pyJ_data, J_data, J.nRows() * J.nCols() * sizeof(double));
  
  PyObject *res = Py_BuildValue("O", pyJ);
  Py_DECREF(pyJ);
  return res;
}


static PyObject *toast_gradient(PyObject *self, PyObject *args) {
  int hmesh, hraster;
  double freq;
  QMMesh *mesh;
  Raster *raster;
  PyArrayObject *py_qvec_vl, *py_qvec_rp, *py_qvec_ci;
  PyArrayObject *py_mvec_vl, *py_mvec_rp, *py_mvec_ci;
  PyArrayObject *py_mua, *py_mus, *py_ref;
  PyArrayObject *py_data, *py_sd;

  const char *solver = "direct";  // for now
  const double tol = 1e-12;       // for now

  if (!PyArg_ParseTuple(args, "iiO!O!O!O!O!O!O!O!O!dO!O!", &hmesh, &hraster,
                        &PyArray_Type, &py_qvec_vl, 
                        &PyArray_Type, &py_qvec_rp, 
                        &PyArray_Type, &py_qvec_ci,
                        &PyArray_Type, &py_mvec_vl, 
                        &PyArray_Type, &py_mvec_rp,
                        &PyArray_Type, &py_mvec_ci,
                        &PyArray_Type, &py_mua, 
                        &PyArray_Type, &py_mus,
                        &PyArray_Type, &py_ref,
                        &freq, 
                        &PyArray_Type, &py_data, 
                        &PyArray_Type, &py_sd)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);
  GETRASTER(raster, hraster);

  int n = mesh->nlen();
  int nQ = mesh->nQ;
  int nM = mesh->nM;
  int nQM = mesh->nQM;

  if (!AssertArray(py_qvec_rp, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "qvec row pointer must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_qvec_ci, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "qvec column indices must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_qvec_vl, 1, NPY_CDOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "qvec vals must be a one-dimensional vector of complex doubles");
    return NULL;
  }
  nint *qrowptr = (nint *)PyArray_DATA(py_qvec_rp);
  nint *qcolidx = (nint *)PyArray_DATA(py_qvec_ci);
  std::complex<double> *qval = (std::complex<double> *)PyArray_DATA(py_qvec_vl);
  CCompRowMatrix qvec(nQ, n, qrowptr, qcolidx, qval /*, SHALLOW_COPY*/);

  if (!AssertArray(py_mvec_rp, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "mvec row pointer must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_mvec_ci, 1, TOAST_NPY_INT)) {
    PyErr_SetString(PyExc_ValueError, "mvec column indices must be a one-dimensional vector of integers");
    return NULL;
  }
  if (!AssertArray(py_mvec_vl, 1, NPY_CDOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "mvec vals must be a one-dimensional vector of complex doubles");
    return NULL;
  }
  nint *mrowptr = (nint *)PyArray_DATA(py_mvec_rp);
  nint *mcolidx = (nint *)PyArray_DATA(py_mvec_ci);
  std::complex<double> *mval = (std::complex<double> *)PyArray_DATA(py_mvec_vl);
  CCompRowMatrix mvec(nM, n, mrowptr, mcolidx, mval /*, SHALLOW_COPY*/);

  if (!AssertArrayDims(py_mua, NPY_DOUBLE, n)) {
    PyErr_SetString(PyExc_ValueError, "mua must be a one-dimensional vector of doubles of length equal to mesh nodes");
    return NULL;
  }
  if (!AssertArrayDims(py_mus, NPY_DOUBLE, n)) {
    PyErr_SetString(PyExc_ValueError, "mua must be a one-dimensional vector of doubles of length equal to mesh nodes");
    return NULL;
  }
  if (!AssertArrayDims(py_ref, NPY_DOUBLE, n)) {
    PyErr_SetString(PyExc_ValueError, "mua must be a one-dimensional vector of doubles of length equal to mesh nodes");
    return NULL;
  }
  int twonqm = nQM * 2;
  if (!AssertArrayDims(py_data, NPY_DOUBLE, twonqm)) {
    PyErr_SetString(PyExc_ValueError, "data must be a one-dimensional vector of doubles of length 2 x source x detector");
    return NULL;
  }
  if (!AssertArrayDims(py_sd, NPY_DOUBLE, twonqm)) {
    PyErr_SetString(PyExc_ValueError, "SD must be a one-dimensional vector of doubles of length 2 x source x detector");
    return NULL;
  }
  double *mua_ptr = (double *)PyArray_DATA(py_mua);
  RVector mua(n, mua_ptr, SHALLOW_COPY);

  double *mus_ptr = (double *)PyArray_DATA(py_mus);
  RVector mus(n, mus_ptr, SHALLOW_COPY);

  double *ref_ptr = (double *)PyArray_DATA(py_ref);
  RVector ref(n, ref_ptr, SHALLOW_COPY);

  double *data_ptr = (double *)PyArray_DATA(py_data);
  RVector data(nQM * 2, data_ptr, SHALLOW_COPY);

  double *sd_ptr = (double *)PyArray_DATA(py_sd);
  RVector sd(nQM * 2, sd_ptr, SHALLOW_COPY);

	int nth;
  #ifdef TOAST_THREAD
  nth = Task::GetThreadCount();
  #else
  nth = 1;
  #endif

  CFwdSolver FWS (mesh, solver, tol, nth);
  FWS.SetDataScaling (DATA_LOG);
  
  npy_intp grad_dim = raster->SLen() * 2;
  PyObject *py_grad = PyArray_SimpleNew(1, &grad_dim, NPY_DOUBLE);
  double *py_grad_data = (double *)PyArray_DATA((PyArrayObject *)py_grad);
  memset(py_grad_data, 0, grad_dim * sizeof(double));
  RVector grad(grad_dim, py_grad_data, SHALLOW_COPY);

  GetGradientCplx (mesh, raster, FWS, mua, mus, ref, freq, data, sd, qvec, mvec, grad);

  PyObject *res = Py_BuildValue("O", py_grad);
  Py_DECREF(py_grad);

  return res;
}

static PyObject *toast_krylov(PyObject *self, PyObject *args) {
  PyObject *py_x, *py_J;

  if (!PyArg_ParseTuple(args, "OO", &py_x, &py_J)) {
    return NULL;
  }

  npy_intp *J_dims = PyArray_DIMS((PyArrayObject *)py_J);
  std::cerr << "Jm=" << J_dims[0] << ", Jn=" << J_dims[1] << std::endl;

  Py_RETURN_NONE;
}

// ===========================================================================
// ===========================================================================
// Regularisation methods

static PyObject *toast_regul(PyObject *self, PyObject *args, PyObject *keywds) {
  Regularisation *reg = 0;
  const char *regtype;
  int hraster, hreg;
  double tau, beta = 1;
  void *kapref = 0;
  bool istensor = false;  // reference diffusivity in tensor format?
  Raster *raster;
  PyArrayObject *py_x;

  static const char *kwlist[] = {"regtype", "raster", "x", "tau", "beta", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "siO!d|d", (char **)kwlist, &regtype, &hraster, 
                                                            &PyArray_Type, &py_x, 
                                                            &tau, &beta)) {
    return NULL;
  }
  GETRASTER(raster, hraster);

  if (!AssertArray(py_x, 1, NPY_DOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "x must be a one-dimensional vector of doubles");
    return NULL;
  }
  npy_intp *x_dims = PyArray_DIMS(py_x);
  int xlen = x_dims[0];

  RVector x0(xlen);
  memcpy(x0.data_buffer(), PyArray_DATA(py_x), xlen * sizeof(double));

  if (!strcasecmp(regtype, "TK0")) {
    RVector xs(x0.Dim());
    xs = 1;
    reg = new Tikhonov0(tau, &x0, &xs);
  } else if (!strcasecmp(regtype, "TK1")) {
    reg = new TK1(tau, &x0, raster, kapref, istensor);
  } else if (!strcasecmp(regtype, "TV")) {
    reg = new TV(tau, beta, &x0, raster, 0, false);
  }

  hreg = g_regmgr.Add(reg);
  return Py_BuildValue("i", hreg);
}

// ===========================================================================

static PyObject *toast_regul_value(PyObject *self, PyObject *args) {
  int hreg;
  Regularisation *reg;
  PyArrayObject *py_x;

  if (!PyArg_ParseTuple(args, "iO!", &hreg, &PyArray_Type, &py_x)) {
    return NULL;
  }
  GETREGUL(reg, hreg);

  if (!AssertArray(py_x, 1, NPY_DOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "x must be a one-dimensional vector of doubles");
    return NULL;
  }
  npy_intp *x_dims = PyArray_DIMS(py_x);
  int xlen = x_dims[0];

  RVector x(xlen);
  memcpy(x.data_buffer(), PyArray_DATA(py_x), xlen * sizeof(double));

  double rval = reg->GetValue(x);
  return Py_BuildValue("d", rval);
}

// ===========================================================================

static PyObject *toast_regul_grad(PyObject *self, PyObject *args) {
  int hreg;
  Regularisation *reg;
  PyArrayObject *py_x;

  if (!PyArg_ParseTuple(args, "iO!", &hreg, &PyArray_Type, &py_x)) {
    return NULL;
  }

  GETREGUL(reg, hreg);

  if (!AssertArray(py_x, 1, NPY_DOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "x must be a one-dimensional vector of doubles");
    return NULL;
  }
  npy_intp *x_dims = PyArray_DIMS(py_x);
  int xlen = x_dims[0];

  RVector x(xlen);
  memcpy(x.data_buffer(), PyArray_DATA(py_x), xlen * sizeof(double));

  RVector grad(xlen);
  grad = reg->GetGradient(x);

  npy_intp grad_dim = grad.Dim();
  PyObject *py_grad = PyArray_SimpleNew(1, &grad_dim, NPY_DOUBLE);
  double *py_grad_data = (double *)PyArray_DATA((PyArrayObject *)py_grad);
  const double *grad_data = grad.data_buffer();
  memcpy(py_grad_data, grad_data, grad.Dim() * sizeof(double));

  PyObject *res = Py_BuildValue("O", py_grad);
  Py_DECREF(py_grad);
  return res;
}

// ===========================================================================

static PyObject *toast_regul_hdiag(PyObject *self, PyObject *args) {
  int hreg;
  Regularisation *reg;
  PyArrayObject *py_x;
  PyObject *py_hdiag;

  if (!PyArg_ParseTuple(args, "iO!", &hreg, &PyArray_Type, &py_x)) {
    return NULL;
  }

  GETREGUL(reg, hreg);

  if (!AssertArray(py_x, 1, NPY_DOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "x must be a (one-dimensional) vector of doubles");
    return NULL;
  }
  RVector x = CopyVector((PyObject *) py_x);
  RVector diag(x.Dim());
  diag = reg->GetHessianDiag(x);
  CopyVector(&py_hdiag, diag);

  PyObject *res = Py_BuildValue("O", py_hdiag);
  Py_DECREF(py_hdiag);
  return res;
}

// ===========================================================================

static PyObject *toast_element_dof(PyObject *self, PyObject *args) {
  int hmesh, elid;
  QMMesh *mesh;

  if (!PyArg_ParseTuple(args, "ii", &hmesh, &elid)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);
  if (elid < 0 || elid >= mesh->nlen()) {
    return NULL;
  }

  Element *pel = mesh->elist[elid];
  npy_intp nnd = (npy_intp)pel->nNode();

  PyObject *pydof = PyArray_SimpleNew(1, &nnd, TOAST_NPY_INT);
  nint *data = (nint *)PyArray_DATA((PyArrayObject *)pydof);

  for(int i = 0; i < nnd; nnd++)
  {
    data[i] = pel->Node[i];
  }
  //memcpy(data, pel->Node, nnd * sizeof(int));

  return Py_BuildValue("N", pydof);
}

// ===========================================================================

static PyObject *toast_element_size(PyObject *self, PyObject *args) {
  int hmesh, elid;
  QMMesh *mesh;

  if (!PyArg_ParseTuple(args, "ii", &hmesh, &elid)) {
    return NULL;
  }
  GETMESH(mesh, hmesh);

  if (elid < 0) {
    // run over entire mesh
    npy_intp nel = (npy_intp)mesh->elen();
    PyObject *pyelsize = PyArray_SimpleNew(1, &nel, NPY_DOUBLE);
    double *data = (double *)PyArray_DATA((PyArrayObject *)pyelsize);
    for (int i = 0; i < mesh->elen(); i++) {
      data[i] = mesh->ElSize(i);
    }
    return Py_BuildValue("N", pyelsize);
  } else {
    // single element
    if (elid >= mesh->nlen()) {
      return NULL;
    }
    return Py_BuildValue("d", mesh->ElSize(elid));
  }
}

// ===========================================================================

static PyObject *toast_element_region(PyObject *self, PyObject *args) {
  int hmesh, elid;
  QMMesh *mesh;

  if (!PyArg_ParseTuple(args, "ii", &hmesh, &elid)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  if (elid < 0) {
    // run over entire mesh
    npy_intp nel = (npy_intp)mesh->elen();
    PyObject *pyelreg = PyArray_SimpleNew(1, &nel, TOAST_NPY_INT);
    nint *data = (nint *)PyArray_DATA((PyArrayObject *)pyelreg);
    for (int i = 0; i < mesh->elen(); i++) {
      data[i] = mesh->elist[i]->Region();
    }
    return Py_BuildValue("N", pyelreg);
  } else {  // single element
    if (elid >= mesh->elen()) {
      return NULL;
    }
    return Py_BuildValue("i", (int) mesh->elist[elid]->Region());
  }
}

// ===========================================================================

static PyObject *toast_element_setregion(PyObject *self, PyObject *args) {
  int hmesh, elid, reg;
  QMMesh *mesh;

  if (!PyArg_ParseTuple(args, "iii", &hmesh, &elid, &reg)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  if (elid >= 0 && elid < mesh->elen()) {
    mesh->elist[elid]->SetRegion(reg);
  }

  Py_RETURN_NONE;
}

// ===========================================================================

static PyObject *toast_mesh_setregion(PyObject *self, PyObject *args) {
  int hmesh;
  QMMesh *mesh;
  PyArrayObject *py_reglist;

  if (!PyArg_ParseTuple(args, "iO!", &hmesh, &PyArray_Type, &py_reglist)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  std::cout << "Warning, dimensionality not checked" << std::endl;

  npy_intp *dims = PyArray_DIMS(py_reglist);
  int *reg = (int *)PyArray_DATA(py_reglist);

  if (dims[0] * dims[1] < mesh->elen()) {
    return NULL;
  }

  for (int i = 0; i < mesh->elen(); i++) {
    mesh->elist[i]->SetRegion(reg[i]);
  }

  Py_RETURN_NONE;
}

// ===========================================================================

static PyObject *toast_element_data(PyObject *self, PyObject *args) {
  int hmesh, elid, i, j;
  QMMesh *mesh;

  if (!PyArg_ParseTuple(args, "ii", &hmesh, &elid)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  std::cout << "Warning, dimensionality not checked" << std::endl;

  Element *pel = mesh->elist[elid];
  int dim = pel->Dimension();
  int nnd = pel->nNode();
  int nsd = pel->nSide();
  int nsn = pel->nSideNode(0);
  for (i = 1; i < nsd; i++) {
    nsn = ::max(nsn, pel->nSideNode(i));
  }

  npy_intp node_dims[2] = {nnd, dim};
  PyObject *nodelist = PyArray_SimpleNew(2, node_dims, NPY_DOUBLE);
  double *v, *vtx_data = (double *)PyArray_DATA((PyArrayObject *)nodelist);
  for (i = 0, v = vtx_data; i < nnd; i++) {
    for (j = 0; j < dim; j++) {
      *v++ = mesh->nlist[pel->Node[i]][j];
    }
  }

  npy_intp el_dims = nnd;
  PyObject *idx = PyArray_SimpleNew(1, &el_dims, TOAST_NPY_INT);
  nint *e, *el_data = (nint *)PyArray_DATA((PyArrayObject *)idx);
  for (i = 0, e = el_data; i < nnd; i++) {
    *e++ = pel->Node[i];
  }

  return Py_BuildValue("OOi", nodelist, idx, pel->Type());
}

// ===========================================================================

static PyObject *toast_element_mat(PyObject *self, PyObject *args) {
  int hmesh, elid, sideidx, i, j, k, l;
  QMMesh *mesh;
  const char *intstr;
  PyArrayObject *py_prm;
  PyObject *elmat = NULL;

  if (!PyArg_ParseTuple(args, "iisO!i", &hmesh, &elid, &intstr, &PyArray_Type, &py_prm, &sideidx)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  Element *pel = mesh->elist[elid];
  int dim = pel->Dimension();
  int nnd = pel->nNode();

  if (!strcmp(intstr, "F")) {
    npy_intp dims = nnd;
    elmat = PyArray_SimpleNew(1, &dims, NPY_DOUBLE);
    double *data = (double *)PyArray_DATA((PyArrayObject *)elmat);
    for (i = 0; i < nnd; i++) {
      data[i] = pel->IntF(i);
    }
  } else if (!strcmp(intstr, "FF")) {
    npy_intp dims[2] = {nnd, nnd};
    elmat = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double *data = (double *)PyArray_DATA((PyArrayObject *)elmat);
    for (i = 0; i < nnd; i++) {
      data[i * nnd + i] = pel->IntFF(i, i);
      for (j = 0; j < i; j++) {
        data[i * nnd + j] = data[j * nnd + i] = pel->IntFF(i, j);
      }
    }
  } else if (!strcmp(intstr, "FFF")) {
    npy_intp dims[3] = {nnd, nnd, nnd};
    elmat = PyArray_SimpleNew(3, dims, NPY_DOUBLE);
    double *data = (double *)PyArray_DATA((PyArrayObject *)elmat);
    for (i = 0; i < nnd; i++) {
      for (j = 0; j < nnd; j++) {
        for (k = 0; k < nnd; k++) {
          data[(i * nnd + j) * nnd + k] = pel->IntFFF(i, j, k);
        }
      }
    }
  } else if (!strcmp(intstr, "DD")) {
    npy_intp dims[2] = {nnd, nnd};
    elmat = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double *data = (double *)PyArray_DATA((PyArrayObject *)elmat);
    for (i = 0; i < nnd; i++) {
      data[i * nnd + i] = pel->IntDD(i, i);
      for (j = 0; j < i; j++) {
        data[i * nnd + j] = data[j * nnd + i] = pel->IntDD(i, j);
      }
    }
  } else if (!strcmp(intstr, "FD")) {
    npy_intp dims[3] = {nnd, nnd, dim};
    elmat = PyArray_SimpleNew(3, dims, NPY_DOUBLE);
    double *data = (double *)PyArray_DATA((PyArrayObject *)elmat);
    for (i = 0; i < nnd; i++) {
      for (j = 0; j < nnd; j++) {
        RVector fd = pel->IntFD(i, j);
        if (fd.Dim() == 0) {
          return NULL;
        }
        for (k = 0; k < dim; k++) {
          data[i * nnd * dim + j * dim + k] = fd[k];
        }
      }
    }
  } else if (!strcmp(intstr, "FDD")) {
    npy_intp dims[3] = {nnd, nnd, nnd};
    elmat = PyArray_SimpleNew(3, dims, NPY_DOUBLE);
    double *data = (double *)PyArray_DATA((PyArrayObject *)elmat);
    for (i = 0; i < nnd; i++) {
      for (j = 0; j < nnd; j++) {
        for (k = 0; k < nnd; k++) {
          data[i * nnd * nnd + j * nnd + k] = pel->IntFDD(i, j, k);
        }
      }
    }
  } else if (!strcmp(intstr, "dd")) {
    npy_intp dims[4] = {nnd, dim, nnd, dim};
    elmat = PyArray_SimpleNew(4, dims, NPY_DOUBLE);
    double *v, *data = (double *)PyArray_DATA((PyArrayObject *)elmat);
    RSymMatrix intdd = pel->Intdd();
    for (i = 0, v = data; i < nnd; i++) {
      for (j = 0; j < dim; j++) {
        for (k = 0; k < nnd; k++) {
          for (l = 0; l < dim; l++) {
            *v++ = intdd(i * dim + j, k * dim + l);
          }
        }
      }
    }
  } else if (!strcmp(intstr, "PFF")) {
    if ((PyObject *)py_prm == Py_None) {
      return NULL;
    }
    RVector prm = RVector(mesh->nlen(), (double *)PyArray_DATA(py_prm), SHALLOW_COPY);
    RSymMatrix intPFF = pel->IntPFF(prm);
    npy_intp dims[2] = {nnd, nnd};
    elmat = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double *data = (double *)PyArray_DATA((PyArrayObject *) elmat);
    for (i = 0; i < nnd; i++) {
      data[i * nnd + i] = intPFF(i, i);
      for (j = 0; j < i; j++) {
        data[i * nnd + j] = data[j * nnd + i] = intPFF(i, j);
      }
    }
  } else if (!strcmp(intstr, "PDD")) {
    if ((PyObject *) py_prm == Py_None) {
      return NULL;
    }
    RVector prm = RVector(mesh->nlen(), (double *)PyArray_DATA(py_prm), SHALLOW_COPY);
    RSymMatrix intPDD = pel->IntPDD(prm);
    npy_intp dims[2] = {nnd, nnd};
    elmat = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double *data = (double *)PyArray_DATA((PyArrayObject *)elmat);
    for (i = 0; i < nnd; i++) {
      data[i * nnd + i] = intPDD(i, i);
      for (j = 0; j < i; j++) {
        data[i * nnd + j] = data[j * nnd + i] = intPDD(i, j);
      }
    }
  } else if (!strcmp(intstr, "BndF")) {
    npy_intp dims = nnd;
    elmat = PyArray_SimpleNew(1, &dims, NPY_DOUBLE);
    double *data = (double *)PyArray_DATA((PyArrayObject *)elmat);
    if (sideidx >= 0) {
      for (i = 0; i < nnd; i++) {
        data[i] = pel->SurfIntF(i, sideidx);
      }
    } else {
      RVector bndintf = pel->BndIntF();
      for (i = 0; i < nnd; i++) {
        data[i] = bndintf[i];
      }
    }
  } else if (!strcmp(intstr, "BndFF")) {
    int ii, jj;
    npy_intp dims[2] = {nnd, nnd};
    elmat = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double *data = (double *)PyArray_DATA((PyArrayObject *)elmat);
    memset(data, 0, nnd * nnd * sizeof(double));
    if (sideidx >= 0) {
      for (ii = 0; ii < pel->nSideNode(sideidx); ii++) {
        i = pel->SideNode(sideidx, ii);
        data[i * nnd + i] = pel->SurfIntFF(i, i, sideidx);
        for (jj = 0; jj < ii; jj++) {
          j = pel->SideNode(sideidx, jj);
          data[i * nnd + j] = data[j * nnd + i] = pel->SurfIntFF(i, j, sideidx);
        }
      }
    } else {
      for (sideidx = 0; sideidx < pel->nSide(); sideidx++) {
        if (!pel->IsBoundarySide(sideidx)) {
          continue;
        }
        for (ii = 0; ii < pel->nSideNode(sideidx); ii++) {
          i = pel->SideNode(sideidx, ii);
          data[i * nnd + i] += pel->SurfIntFF(i, i, sideidx);
          for (jj = 0; jj < ii; jj++) {
            j = pel->SideNode(sideidx, jj);
            data[i * nnd + j] = data[j * nnd + i] += pel->SurfIntFF(i, j, sideidx);
          }
        }
      }
    }
  } else if (!strcmp(intstr, "BndPFF")) {
    if ((PyObject *) py_prm == Py_None) {
      return NULL;
    }
    RVector prm = RVector(mesh->nlen(), (double *)PyArray_DATA(py_prm),
                          SHALLOW_COPY);
    npy_intp dims[2] = {nnd, nnd};
    elmat = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double *data = (double *)PyArray_DATA((PyArrayObject *)elmat);
    if (sideidx >= 0) {
      PyErr_SetString(PyExc_ValueError, "Not implemented yet");
      return NULL;
    } else {
      for (i = 0; i < nnd; i++) {
        for (j = 0; j < nnd; j++) {
          data[i + j * nnd] = pel->BndIntPFF(i, j, prm);
        }
      }
    }
  }

  if (!elmat) {
    return NULL;
  }

  return Py_BuildValue("O", elmat);
}

// ===========================================================================

static PyObject *toast_element_shapef(PyObject *self, PyObject *args) {
  int i, j, hmesh, elid, global = 0;
  QMMesh *mesh;
  PyArrayObject *py_pos;
  if (!PyArg_ParseTuple(args, "iiO!|i", &hmesh, &elid, &PyArray_Type, &py_pos, &global)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  if (elid < 0 || elid >= mesh->elen()) {
    return NULL;
  }

  if ((PyObject *)py_pos == Py_None) {
    return NULL;
  }

  int dim = mesh->Dimension();

  Element *pel = mesh->elist[elid];
  int nn = pel->nNode();

  double *pos = (double *)PyArray_DATA(py_pos);
  int nd = PyArray_NDIM(py_pos);
  npy_intp *dims = PyArray_DIMS(py_pos);
  if (dims[0] != dim) {
    return NULL;
  }
  int npoint = (nd == 1 ? 1 : dims[1]);

  npy_intp outd[2] = {nn, npoint};
  PyObject *py_shapef = PyArray_SimpleNew(2, outd, NPY_DOUBLE);
  double *shapef = (double *)PyArray_DATA((PyArrayObject *)py_shapef);
  Point pt(dim);
  for (j = 0; j < npoint; j++) {
    for (i = 0; i < dim; i++) {
      pt[i] = pos[i * npoint + j];
    }
    RVector fun = (global ? pel->GlobalShapeF(mesh->nlist, pt) : pel->LocalShapeF(pt));
    for (i = 0; i < nn; i++) {
      shapef[i * npoint + j] = fun[i];
    }
  }

  return Py_BuildValue("O", py_shapef);
}

// ===========================================================================

static PyObject *toast_element_shaped(PyObject *self, PyObject *args) {
  int i, j, k, hmesh, elid, global = 0;
  QMMesh *mesh;
  PyArrayObject *py_pos;
  if (!PyArg_ParseTuple(args, "iiO!|i", &hmesh, &elid, &PyArray_Type, &py_pos, &global)) {
    return NULL;
  }

  GETMESH(mesh, hmesh);

  if (elid < 0 || elid >= mesh->elen()) {
    return NULL;
  }

  int dim = mesh->Dimension();

  Element *pel = mesh->elist[elid];
  int nn = pel->nNode();

  double *pos = (double *)PyArray_DATA(py_pos);
  int nd = PyArray_NDIM(py_pos);
  npy_intp *dims = PyArray_DIMS(py_pos);
  if (dims[0] != dim) {
    return NULL;
  }
  int npoint = (nd == 1 ? 1 : dims[1]);

  npy_intp outd[3] = {nn, dim, npoint};
  PyObject *py_shaped = PyArray_SimpleNew(npoint == 1 ? 2 : 3, outd, NPY_DOUBLE);
  double *shaped = (double *)PyArray_DATA((PyArrayObject *)py_shaped);
  Point pt(dim);
  for (j = 0; j < npoint; j++) {
    for (i = 0; i < dim; i++) {
      pt[i] = pos[i * npoint + j];
    }
    RDenseMatrix fgrad = (global ? pel->GlobalShapeD(mesh->nlist, pt) : pel->LocalShapeD(pt));
    double *pg = fgrad.ValPtr();
    for (i = 0; i < nn; i++) {
      for (k = 0; k < dim; k++) {
        shaped[j * nn * dim + i * dim + k] = pg[i + k * nn];
      }
    }
  }

  return Py_BuildValue("O", py_shaped);
}

// ===========================================================================

static PyObject *toast_test(PyObject *self, PyObject *args) {
  npy_intp dmx_dims[2] = {100, 10000};
  PyObject *dmx = PyArray_SimpleNew(2, dmx_dims, NPY_CDOUBLE);
  PyObject *ret = Py_BuildValue("O", dmx);
  Py_DECREF(dmx);
  return ret;
}

// ===========================================================================

static PyMethodDef ToastMethods[] = {
    {"ReadMesh", mesh_read, METH_VARARGS, "Read a Toast mesh from file"},
    {"WriteMesh", mesh_write, METH_VARARGS, "Write a Toast mesh to file"},
    {"ClearMesh", mesh_clear, METH_VARARGS, "Delete a mesh from memory"},
    {"MeshData", mesh_data, METH_VARARGS, "Extract node and element data from a mesh"},
    {"SurfData", toast_surf_data, METH_VARARGS, "Extract surface node and face data from a mesh"},
    {"meshSetRegion", toast_mesh_setregion, METH_VARARGS, "Set region values for all mesh elements"},
    {"ElementData", toast_element_data, METH_VARARGS, "Extract node data for a mesh element"},
    {"MeshNodeCount", toast_mesh_node_count, METH_VARARGS, "Return the number of mesh nodes"},
    {"MeshElementCount", toast_mesh_element_count, METH_VARARGS, "Return the number of mesh elements"},
    {"MeshDim", toast_mesh_dim, METH_VARARGS, "Return the mesh dimension (2 or 3)"},
    {"MeshBB", toast_mesh_bb, METH_VARARGS, "Return the corner coordinates of the mesh bounding box"},
    {"RasterBasisPoints", toast_basis_pos, METH_VARARGS, "Return the positions of the basis points in a matrix"},
    {"RasterSolutionPoints", toast_sol_pos, METH_VARARGS, "Return the positions of the solution points in a matrix"},
    {"MakeMesh", toast_make_mesh, METH_VARARGS, "Create a Toast mesh from node and element data"},
    {"MakeRaster", toast_make_raster, METH_VARARGS, "Create a mesh-to-raster mapper"},
    {"ClearRaster", toast_clear_raster, METH_VARARGS, "Delete a raster object from memory"},
    {"RasterGLen", toast_raster_glen, METH_VARARGS, "Get a raster glen"},
    {"RasterBLen", toast_raster_blen, METH_VARARGS, "Get a raster blen"},
    {"RasterSLen", toast_raster_slen, METH_VARARGS, "Get a raster slen"},
    {"RasterMatrix", toast_raster_matrix, METH_VARARGS, "Get the sparse matrix parameters for a raster transformation."},
    {"RasterBasis2Sol", toast_raster_basis2sol, METH_VARARGS, "A vector of BLen() containing indices of the corresponding element in the solution vector."},
    {"RasterSol2Basis", toast_raster_sol2basis, METH_VARARGS, "A vector of SLen() containing indices of the corresponding element in the basis vector."},
    {"MapBasis", toast_map_basis, METH_VARARGS, "Map a field from one basis to another"},
    {"ReadQM", toast_readqm, METH_VARARGS, "Load a QM description file into a mesh"},
    {"ReadNim", toast_read_nim, METH_VARARGS, "Read a nodal image file"},
    {"WriteNim", toast_write_nim, METH_VARARGS, "Write a nodal image to file"},
    {"Sysmat", toast_sysmat, METH_VARARGS, "Create a complex sparse system matrix from optical parameters"},
    {"Sysmat_CW", toast_sysmat_cw, METH_VARARGS, "Create a real sparse system matrix from optical parameters"},
    {"Qvec", (PyCFunction)toast_qvec, METH_VARARGS | METH_KEYWORDS, "Construct an array of source vectors from QM information"},
    {"Mvec", (PyCFunction)toast_mvec, METH_VARARGS | METH_KEYWORDS, "Construct an array of measurement vectors from QM information"},
    {"Fields", toast_fields, METH_VARARGS, "Calculate direct and adjoint fields"},
    {"Jacobian", toast_jacobian, METH_VARARGS, "Calculate the Jacobian matrix of the DOT forward operator"},
    {"JacobianOptical", toast_jacobian_optical, METH_VARARGS, "Calculate the Jacobian matrix of the DOT forward operator from optical parameters."},
    {"JacobianCW", toast_jacobianCW, METH_VARARGS, "Calculate the Jacobian matrix of the CW intensity data (real case)"},
    {"Gradient", toast_gradient, METH_VARARGS, "Calculate the gradient of the DOT forward operator from the optical parameters"},
    {"Krylov", toast_krylov, METH_VARARGS, "Solve Hx = y with implicit Hessian H = J^T J"},
    {"Regul", (PyCFunction)toast_regul, METH_VARARGS | METH_KEYWORDS, "Return a regularisation object"},
    {"RegValue", toast_regul_value, METH_VARARGS, "Returns the regularisation value R(x) for a parameter vector x"},
    {"RegGradient", toast_regul_grad, METH_VARARGS, "Returns the regularisation gradient for a parameter vector x"},
    {"RegHDiag", toast_regul_hdiag, METH_VARARGS, "Returns the diagonal of the Hessian of the regularisation operator"},
    {"elementDof", toast_element_dof, METH_VARARGS, "Returns a permutation array for the global degrees of freedom of the element"},
    {"elementSize", toast_element_size, METH_VARARGS, "Returns the element size"},
    {"elementRegion", toast_element_region, METH_VARARGS, "Returns the element region index"},
    {"elementSetRegion", toast_element_setregion, METH_VARARGS, "Sets the element region index"},
    {"elementData", toast_element_data, METH_VARARGS, "Returns the element geometry"},
    {"elementMat", toast_element_mat, METH_VARARGS, "Returns an integral over the element or element surface"},
    {"elementShapeF", toast_element_shapef, METH_VARARGS, "Returns shape function values for points in the element"},
    {"elementShapeD", toast_element_shaped, METH_VARARGS, "Returns shape function derivatives for points in the element"},

    {"Test", toast_test, METH_VARARGS, "A dummy test function"},
    {NULL, NULL, 0, NULL}};

// ===========================================================================

static int toastmod_traverse(PyObject *m, visitproc visit, void *arg) {
  Py_VISIT(GETSTATE(m)->error);
  return 0;
}

static int toastmod_clear(PyObject *m) {
  Py_CLEAR(GETSTATE(m)->error);
  return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "toastmod",
    NULL,
    sizeof(struct module_state),
    ToastMethods,
    NULL,
    toastmod_traverse,
    toastmod_clear,
    NULL};

PyMODINIT_FUNC PyInit_toastmod(void) {
  PyObject *module = PyModule_Create(&moduledef);
  import_array();

#ifdef TOAST_THREAD
	Task_Init (0);
#endif

  if (module == NULL)
    return NULL;
  struct module_state *st = GETSTATE(module);

  st->error = PyErr_NewException("toastmod.Error", NULL, NULL);
  if (st->error == NULL) {
    Py_DECREF(module);
    return NULL;
  }

  return module;
}