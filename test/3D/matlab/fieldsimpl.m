% fieldsimpl.m
%
% Compare performance and cross-check of toastFields using MATLAB and TOAST solver
% implementations.

close all 
clear all

% Build a system with quite a few RHS
meshdir = '../meshes/';
hmesh = toastMesh([meshdir 'cyl3.msh']);
hmesh.ReadQM([meshdir 'cyl_5ring.qm']);
nnd = hmesh.NodeCount;
mua = 0.01 * ones(nnd,1);
mus = 1.00 * ones(nnd,1);
ref = 1.4 * ones(nnd,1);
qvec = hmesh.Qvec('Neumann', 'Gaussian', 1);
mvec = hmesh.Mvec('Gaussian', 1, ref);

% CW
tic; x1  = toastFields(hmesh,0,qvec,mua,mus,ref,0,'direct',1e-10,'matlab'); t_ml_cw_direct = toc
tic; x2  = toastFields(hmesh,0,qvec,mua,mus,ref,0,'direct',1e-10,'toast'); t_st_cw_direct = toc
tic; x3  = toastFields(hmesh,0,qvec,mua,mus,ref,0,'direct',1e-10,'auto'); t_au_cw_direct = toc
tic; x4  = toastFields(hmesh,0,qvec,mua,mus,ref,0,'cg',1e-10,'matlab'); t_ml_cw_cg = toc
tic; x5  = toastFields(hmesh,0,qvec,mua,mus,ref,0,'cg',1e-10,'toast'); t_st_cw_cg = toc
tic; x6  = toastFields(hmesh,0,qvec,mua,mus,ref,0,'cg',1e-10,'auto'); t_au_cw_cg = toc

assert(norm(x1 - x2) < 1e-8);
assert(norm(x1 - x3) < 1e-8);
assert(norm(x1 - x4) < 1e-8);
assert(norm(x1 - x5) < 1e-8);
assert(norm(x1 - x6) < 1e-8);

% FD
tic; x7  = toastFields(hmesh,0,qvec,mua,mus,ref,80,'direct',1e-10,'matlab'); t_ml_fd_direct = toc
tic; x8  = toastFields(hmesh,0,qvec,mua,mus,ref,80,'direct',1e-10,'toast'); t_st_fd_direct = toc
tic; x9  = toastFields(hmesh,0,qvec,mua,mus,ref,80,'direct',1e-10,'auto'); t_au_fd_direct = toc
tic; x10 = toastFields(hmesh,0,qvec,mua,mus,ref,80,'bicgstab',1e-10,'matlab'); t_ml_fd_cg = toc
tic; x11 = toastFields(hmesh,0,qvec,mua,mus,ref,80,'bicgstab',1e-10,'toast'); t_st_fd_cg = toc
tic; x12 = toastFields(hmesh,0,qvec,mua,mus,ref,80,'bicgstab',1e-10,'auto'); t_au_fd_cg = toc

assert(norm(x7 - x8) < 1e-8);
assert(norm(x7 - x9) < 1e-8);
assert(norm(x7 - x10) < 1e-8);
assert(norm(x7 - x11) < 1e-8);
assert(norm(x7 - x12) < 1e-8);



