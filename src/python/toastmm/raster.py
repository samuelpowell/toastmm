import numpy as np
from toastmm import toastmod
from scipy import sparse
from types import *


class Raster:
  """Basis mapping object.

  Syntax basis = toast.raster.Raster(mesh,grd)

  Parameters:
    mesh:   toast.mesh.Mesh object containing a valid FEM mesh.
    grd:    integer array of length 2 or 3 (corresponding to mesh dimension) containing the
            grid size for the regular basis.
    intgrd: integer array of length 2 or 3 (corresponding to mesh dimension) containing the
            grid size for the intermediate basis.

  Notes:
    A raster object allows to map between an unstructured FEM nodal basis, and a regular
    grid basis. Typically this is used when the forward problem is defined as a FEM solver,
    and the inverse problem is solved on a regular grid.
  """

  def __init__(self, mesh, grd, intgrd=None):
    self.handle = None
    self.make(mesh, grd, intgrd=intgrd)

  def __del__(self):
    self.clear()

  # def handle(self):
  #   """Returns the internal raster handle.
  #
  #   Syntax: handle = raster.handle()
  #   """
  #   return self.handle

  def make(self, mesh, grd, intgrd=None):
    """Initialise the mapper object by assigning a mesh and grid.

    Syntax: raster.make(mesh, grd, intgrd)

    Parameters:
      mesh:   toast.mesh.Mesh object
      grd:    integer array of length 2 or 3 (corresponding to mesh dimension) containing 
              the grid size for the regular basis.
      intgrd: integer array of length 2 or 3 (corresponding to mesh dimension) containing
              the grid size for the intermediate basis.
    """
    self.mesh = mesh
    grd = np.array(grd, dtype=np.intc)
    self.grd = grd

    if intgrd is None:
      intgrd = np.copy(grd)
    else:
      intgrd = np.array(intgrd, dtype=np.intc)

    self.clear()
    self.handle = toastmod.MakeRaster(mesh.handle, grd, intgrd)

  def clear(self):
    if self.handle is not None:
      toastmod.ClearRaster(self.handle)
      self.handle = None

  def map(self, mapstr, srcvec):
    """Map a scalar field from one basis to another.

    Syntax: tgt_coef = raster.map(mapstr, srcvec)

    Parameters:
      mapstr: a string of the form 'S->T' defining the source basis (S) to map from, and
              target basis (T) to map to. "S" and "T" are placeholders for one of:
                M: mesh basis
                B: raster basis (fully populated bounding box)
                S: raster basis (sparse; omitting voxels with no mesh support)
      srcvec: array of basis coefficients in source basis

    Return values:
      tgt_coef: array of basis coefficients in target basis
    """
    if srcvec.size > max(srcvec.shape):
      raise ValueError("Only vectors are supported.")
    return toastmod.MapBasis(self.handle, mapstr, srcvec.flatten())

  def basis_points(self):
    return toastmod.RasterBasisPoints(self.handle)

  def solution_points(self):
    return toastmod.RasterSolutionPoints(self.handle)

  def glen(self):
    return toastmod.RasterGLen(self.handle)

  def blen(self):
    return toastmod.RasterBLen(self.handle)

  def slen(self):
    return toastmod.RasterSLen(self.handle)

  def sol_to_basis(self):
    """Returns a vector of length Slen() (solution dimension). Each entry contains the index
    of the corresponding element in the basis vector (range 0 .. Blen()-1). Used to map
    between basis representation (full bounding box) and solution representation 
    (supported|non-masked voxels only).
    """
    return toastmod.RasterSol2Basis(self.handle)

  def basis_to_sol(self):
    """ Returns a vector of length Blen() (basis dimension). Each entry contains the index 
    of the corresponding element in the solution vector (range 0 .. Slen()-1), or -1 if the
    element is outside the support of the domain.
    """
    return toastmod.RasterBasis2Sol(self.handle)

  def get_matrix(self, mapstr):
    src_key = mapstr[0]
    trg_key = mapstr[-1]
    mat_sizes = {
        'S': self.slen(),
        'G': self.glen(),
        'B': self.blen(),
        'M': self.mesh.node_count()
    }

    # Note: If basis and grid are identical, then matrices G and GI are null. Return identity mat.
    if mat_sizes[src_key] == mat_sizes[trg_key] and src_key in ['G', 'B'
                                                               ] and trg_key in ['G', 'B']:
      return sparse.identity(self.GLen())

    rp, ci, vl = toastmod.RasterMatrix(self.handle, mapstr)
    return sparse.csr_matrix((vl, ci, rp), shape=(mat_sizes[trg_key], mat_sizes[src_key]))
