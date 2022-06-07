% perf_1
%
% Measure performance of TOAST library in building a mesh and constructing a Jacobian for
% a representatively large 3D head mesh, which has undergone a dirty affine transform.
%

clear all
close all

load head.mat

% Some ancillary parameters
nnd = size(node,1);
nel = size(elem,1);
eltp = 3*ones(nel,1);
ref = 1.4*ones(nnd,1);
mua = 0.01 *ones(nnd,1);
mus = 1.00 *ones(nnd,1);

tol = 1e-10;

% Build the mesh
hmesh = toastMesh(node, elem, eltp);

if any(hmesh.ElementSize() < 0)
    warning('Negative volume elements');
end

% Build the source and measure vectors
hmesh.SetQM(qpos, mpos);
qvec = hmesh.Qvec ('Neumann', 'Gaussian', 2);
mvec = hmesh.Mvec ('Gaussian', 2, ref);

% Build a basis
hbasis = toastBasis([20 20 20]);

% Matrix of perf
%
% - Solver
% - Implementation
% - Basis/mesh

% Continuous wave

phi_cw_direct_lib = toastFields(hmesh,0,qvec,mua,mus,ref,0,'direct',tol,'toast');
phi_cw_direct_mat = toastFields(hmesh,0,qvec,mua,mus,ref,0,'direct',tol,'matlab');
phi_cw_cg_lib = toastFields(hmesh,0,qvec,mua,mus,ref,0,'cg',tol,'toast');
phi_cw_cg_mat = toastFields(hmesh,0,qvec,mua,mus,ref,0,'cg',tol,'matlab');

% Frequency domain
freq = 100;

phi_fd_direct_lib = toastFields(hmesh,0,qvec,mua,mus,ref,freq,'direct',tol,'toast');
phi_fd_direct_mat = toastFields(hmesh,0,qvec,mua,mus,ref,freq,'direct',tol,'matlab');
phi_fd_bc_lib = toastFields(hmesh,0,qvec,mua,mus,ref,freq,'bicgstab',tol,'toast');
phi_fd_bc_mat = toastFields(hmesh,0,qvec,mua,mus,ref,freq,'bicgstab',tol,'matlab');
phi_fd_gm_lib = toastFields(hmesh,0,qvec,mua,mus,ref,freq,'gmres',tol,'toast');
phi_fd_gm_mat = toastFields(hmesh,0,qvec,mua,mus,ref,freq,'gmres',tol,'matlab');
