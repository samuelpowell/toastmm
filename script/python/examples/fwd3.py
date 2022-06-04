# Solve a 3D problem, computing the fields

# Import various modules
import os
import numpy as np
import toast

# Set the file paths
meshdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "meshes", "3D")
meshfile = os.path.join(meshdir, "cyl3.msh")
qmfile = os.path.join(meshdir, "cyl_3ring.qm")

# Load the mesh and source/detector specs
mesh = toast.Mesh(meshfile)
mesh.ReadQM(qmfile)
nlen = mesh.NodeCount()

# Extract mesh geometry
nlist,elist,eltp = mesh.Data()

# Homogeneous parameter distributions
refind = 1.4
mua = np.ones (nlen) * 0.025
mus = np.ones (nlen) * 2.0
ref = np.ones (nlen) * refind
freq = 100

# Build the source vectors
qvec = mesh.Qvec(type='Neumann', shape='Gaussian', width=2)
mvec = mesh.Mvec(shape='Gaussian', width=2, ref=refind)
nq = qvec.shape[1]

# Compute the forward fields
phi = mesh.Fields(None, qvec, mua, mus, ref, freq)



