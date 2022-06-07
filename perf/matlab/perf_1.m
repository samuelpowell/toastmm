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

% The matrix of performance:
%
% Library / MATLAB
%  - DIRECT / CG
%  - Basis / Mesh

% Mesh Jacobian all styles
tic
JmDM = toastJacobianCW(hmesh, 0, qvec, mvec, mua, mus, ref, 'direct', tol, 'matlab');
t_JmDM = toc

tic
JmIM = toastJacobianCW(hmesh, 0, qvec, mvec, mua, mus, ref, 'cg', tol, 'matlab');
t_JmIM = toc

% tic
% JmDT = toastJacobianCW(hmesh, 0, real(qvec), real(mvec), mua, mus, ref, 'direct', tol, 'toast');
% t_JmDT = toc
% 
tic
JmCT = toastJacobianCW(hmesh, 0, qvec, mvec, mua, mus, ref, 'cg', tol, 'toast');
t_JmCT = toc

% Basis Jacobian all styles
hbasis = toastBasis(hmesh, [20 20 20]);

tic
JbDM = toastJacobianCW(hmesh, hbasis, qvec, mvec, mua, mus, ref, 'direct', tol, 'matlab');
t_JmDM = toc

tic
JbIM = toastJacobianCW(hmesh, hbasis, qvec, mvec, mua, mus, ref, 'cg', tol, 'matlab');
t_JmIM = toc

tic
JbDT = toastJacobianCW(hmesh, hbasis, real(qvec), real(mvec), mua, mus, ref, 'direct', tol, 'toast');
t_JmDT = toc

tic
JbCT = toastJacobianCW(hmesh, hbasis, real(qvec), real(mvec), mua, mus, ref, 'cg', tol, 'toast');
t_JmCT = toc



