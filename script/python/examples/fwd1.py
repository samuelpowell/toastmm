# This pytoast example solves the forward problem
# for a homogeneous 2D circular problem

# Import various modules
import os
import sys
import numpy as np
from numpy import matrix
from scipy import sparse
from scipy.sparse import linalg
from numpy.random import rand
import matplotlib.pyplot as plt

import toastmm as ts

# Set the file paths
meshdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "meshes", "2D")
meshfile = os.path.join(meshdir, "circle25_32.msh")
qmfile = os.path.join(meshdir, "circle25_32x32.qm")

# Load the mesh and source/detector specs
mesh = ts.Mesh(meshfile)
mesh.read_qm(qmfile)
nlen = mesh.node_count()

# Homogeneous parameter distributions
mua = np.ones(nlen) * 0.025
mus = np.ones(nlen) * 2.0
ref = np.ones(nlen) * 1.4
freq = 100

# Set up the linear system
smat = mesh.sysmat(mua, mus, ref, freq)
qvec = mesh.qvec()
mvec = mesh.mvec()

# Solve the linear system
nq = qvec.shape[1]
phi = np.empty(qvec.shape, dtype='complex128')
for q in range(nq):
  qq = qvec[:, q].todense()
  res = linalg.gmres(smat, qq, tol=1e-12)
  phi[:, q] = res[0]

# Project to boundary
y = mvec.transpose() * phi
logy = np.log(y)

# Display as sinogram
plt.figure(1)
im = plt.imshow(logy.real, interpolation='none')
plt.title('log amplitude')
plt.xlabel('detector index')
plt.ylabel('source index')
plt.colorbar()
plt.draw()
#plt.show()

plt.figure(2)
im = plt.imshow(logy.imag, interpolation='none')
plt.title('phase')
plt.xlabel('detector index')
plt.ylabel('source index')
plt.colorbar()
plt.show()
