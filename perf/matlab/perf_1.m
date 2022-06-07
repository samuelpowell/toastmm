% perf_1
%
% Measure performance of TOAST library in building a mesh and constructing a Jacobian for
% a representatively large 3D head mesh, which has undergone a dirty affine transform.
%

clear all
close all

load head.mat

% Build the mesh
nnd = size(node,1);
eltp = 3*ones(nnd,1);
elem = elem(:,[4 1:3]);
hmesh = toastMesh(node, elem, eltp);