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
hbasis = toastBasis(hmesh, [20 20 20]);

% Continuous wave

% Compute the fields to exclude this from the computation timing
phi_cw = toastFields(hmesh,0,[qvec mvec],mua,mus,ref,0,'cg',tol);
dphi_cw = phi_cw(:, 1:size(qvec,2));
aphi_cw = phi_cw(:, (size(qvec,2)+1) : end);
        
% Build projection data, reduce by linklist
proj_cw = reshape (mvec.' * dphi_cw, [], 1);
proj_cw = proj_cw(hmesh.DataLinkList());

tic
J_mCW = toastJacobianCW(hmesh, 0, dphi_cw, aphi_cw, proj_cw);
t_mCW = toc

tic
J_bCW = toastJacobianCW(hmesh, hbasis, dphi_cw, aphi_cw, proj_cw);
t_bCW = toc

% Frequency domain

% Compute the fields to exclude this from the computation timing
phi_fd = toastFields(hmesh,100,[qvec mvec],mua,mus,ref,0,'bicgstab',tol);
dphi_fd = phi_fd(:, 1:size(qvec,2));
aphi_fd = phi_fd(:, (size(qvec,2)+1) : end);
        
% Build projection data, reduce by linklist
proj_fd = reshape (mvec.' * dphi_fd, [], 1);
proj_fd = proj_fd(hmesh.DataLinkList());

tic
J_mFD = toastJacobianCW(hmesh, 0, dphi_fd, aphi_fd, proj_fd);
t_mFD = toc

tic
J_bFD = toastJacobianCW(hmesh, hbasis, dphi_fd, aphi_fd, proj_fd);
t_bFD = toc



