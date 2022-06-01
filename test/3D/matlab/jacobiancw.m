% fieldsimpl.m
%
% Compare performance and cross-check of toastFields using MATLAB and TOAST solver
% implementations.

close all 
clear all

% Build a system with quite a few RHS
meshdir = '../meshes/';
hmesh = toastMesh([meshdir 'cyl3.msh']);
hmesh.ReadQM([meshdir 'cyl_3ring.qm']);
nnd = hmesh.NodeCount;
mua = 0.01 * ones(nnd,1);
mus = 1.00 * ones(nnd,1);
ref = 1.4 * ones(nnd,1);
qvec = hmesh.Qvec('Neumann', 'Gaussian', 1);
mvec = hmesh.Mvec('Gaussian', 1, ref);
hbasis = toastBasis(hmesh, [100 100 100]);

tic; J1 = toastJacobianCW(hmesh, 0, qvec, mvec, mua, mus, ref); t_J1 = toc;
tic; J2 = toastJacobianCW(hmesh, 0, qvec, mvec, mua, mus, ref, 'cg'); t_J2 = toc;
tic; J6 = toastJacobianCW(hmesh, hbasis, qvec, mvec, mua, mus, ref); t_J1 = toc;
tic; J7 = toastJacobianCW(hmesh, hbasis, qvec, mvec, mua, mus, ref, 'cg'); t_J2 = toc;
