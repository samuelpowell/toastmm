# This pytoast example solves the forward problem
# for a homogeneous 2D circular problem

# Import various modules
import os
import numpy as np
from numpy import matrix
from scipy import sparse
from scipy.sparse import linalg
from numpy.random import rand
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import toastmm as ts

# Set the file paths
meshdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "meshes", "2D")
meshfile = os.path.join(meshdir, "circle25_32.msh")
qmfile = os.path.join(meshdir, "circle25_32x32.qm")

# Load the mesh and source/detector specs
mesh = ts.Mesh(meshfile)
mesh.read_qm(qmfile)
nlen = mesh.node_count()

# Extract mesh geometry
nlist, elist, eltp = mesh.data()

grd = np.array([64, 64])
basis = ts.Raster(mesh, grd)

# Homogeneous parameter distributions
refind = 1.4
mua = np.ones(nlen) * 0.025
mus = np.ones(nlen) * 2.0
ref = np.ones(nlen) * refind
freq = 100

# Set up the linear system
qvec = mesh.qvec(type='Neumann', shape='Gaussian', width=2)
mvec = mesh.mvec(shape='Gaussian', width=2, ref=refind)
nq = qvec.shape[1]

phi = mesh.fields(basis, qvec, mua, mus, ref, freq)

fig1 = plt.figure()
ims = []
for q in range(nq):
  lphi = np.log(phi[:, q])
  bphi = basis.map('S->B', lphi)
  bphi = np.reshape(bphi, grd)
  ims.append((plt.imshow(bphi.real),))

im_ani = animation.ArtistAnimation(fig1, ims, interval=50, repeat=False, blit=True)
plt.show()

# Display fields
fig2 = plt.figure()
plt.subplot(1, 2, 1)
im = plt.imshow(bphi.real)
plt.title('log amplitude')
plt.colorbar()
plt.subplot(1, 2, 2)
im = plt.imshow(bphi.imag)
plt.title('phase')
plt.colorbar()
plt.show()
