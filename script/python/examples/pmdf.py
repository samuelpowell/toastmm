# This pytoast example generates the Jacobian matrix
# of the discrete forward operator for a 2D circular problem

# Import various modules
import os
import numpy as np
from numpy import matrix
from scipy import sparse
from scipy.sparse import linalg
from numpy.random import rand
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Import the toast modules
import toastmm as ts

# Set the file paths
meshdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "meshes", "2D")
meshfile = os.path.join(meshdir, "ellips_tri10.msh")
qmfile = os.path.join(meshdir, "circle25_32x32.qm")
muafile = os.path.join(meshdir, "tgt_mua_ellips_tri10.nim")
musfile = os.path.join(meshdir, "tgt_mus_ellips_tri10.nim")

# Load the mesh and source/detector specs
mesh = ts.Mesh(meshfile)
mesh.read_qm(qmfile)

# Extract mesh geometry
nlist, elist, eltp = mesh.data()
nlen = nlist.shape[0]

grd = np.array([64, 64])
basis = ts.Raster(mesh, grd)

# Set up the linear system
qvec = mesh.qvec()
mvec = mesh.mvec()
nq = qvec.shape[0]

# Homogeneous parameter distributions
for bkg in range(2):
  if bkg == 0:
    mua = np.ones(nlen) * 0.025
    mus = np.ones(nlen) * 2.0
  else:
    mua = mesh.read_nim(muafile)
    mus = mesh.read_nim(musfile)

  ref = np.ones(nlen) * 1.4
  freq = 100

  # Calculate fields and projections
  dphi = mesh.fields(None, qvec, mua, mus, ref, freq)
  aphi = mesh.fields(None, mvec, mua, mus, ref, freq)
  proj = np.reshape(mvec.T * dphi, (-1), 'F')

  # Calculate Jacobian matrix
  J = mesh.jacobian(basis, dphi, aphi, proj)

  # Extract sensitivity regions for a single source-detector pair
  slen = int(J.shape[1] / 2)
  nqm = int(J.shape[0] / 2)
  J8_lnamp = J[10, :]
  J8_phase = J[10 + nqm]
  J8_lnamp_mua = J8_lnamp[0:slen - 1]
  J8_lnamp_kap = J8_lnamp[slen:slen * 2 - 1]
  J8_phase_mua = J8_phase[0:slen - 1]
  J8_phase_kap = J8_phase[slen:slen * 2 - 1]

  bJ8_lnamp_mua = basis.map('S->B', J8_lnamp_mua)
  bJ8_lnamp_mua = np.reshape(bJ8_lnamp_mua, grd)
  bJ8_lnamp_kap = basis.map('S->B', J8_lnamp_kap)
  bJ8_lnamp_kap = np.reshape(bJ8_lnamp_kap, grd)

  bJ8_phase_mua = basis.map('S->B', J8_phase_mua)
  bJ8_phase_mua = np.reshape(bJ8_phase_mua, grd)
  bJ8_phase_kap = basis.map('S->B', J8_phase_kap)
  bJ8_phase_kap = np.reshape(bJ8_phase_kap, grd)

  # Display sensitivity regions as images
  title = ["Homogeneous background", "Inhomogeneous background"]
  plt.figure("PMDF: " + title[bkg])
  plt.subplot(2, 2, 1)
  im = plt.imshow(bJ8_lnamp_mua)
  plt.title('mua, lnamp')
  plt.colorbar()

  plt.subplot(2, 2, 2)
  im = plt.imshow(bJ8_lnamp_kap)
  plt.title('kappa, lnamp')
  plt.colorbar()

  plt.subplot(2, 2, 3)
  im = plt.imshow(bJ8_phase_mua)
  plt.title('mua, phase')
  plt.colorbar()

  plt.subplot(2, 2, 4)
  im = plt.imshow(bJ8_phase_kap)
  plt.title('kappa, phase')
  plt.colorbar()

plt.show()
