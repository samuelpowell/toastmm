function proj = toastProject(hMesh,mua,mus,ref,omega,qvec,mvec,solver,tol)
%toastProject         - Forward operator to calulate boundary data.
%
% Synopsis: proj = toastProject (hMesh,mua,mus,ref,omega,qvec,mvec, ...
%                                solver,tol)
%    hMesh: mesh handle
%    mua:   vector of nodal absorption values [1/mm]
%    mus:   vector of nodal scatter values [1/mm]
%    ref:   vector of nodal refractive index values
%    omega: modulation frequency [MHz]
%    qvec:  column matrix of nodal source distributions
%    mvec:  column matrix of nodal measurement distributions
%    solver: linear solver [DIRECT|CG|BICG|BICGSTAB|GMRES]
%    tol:   linear solver tolerance (ignored for solver DIRECT)
%    proj:  vector of boundary measurements.
%
% Performs forward solutions for all sources, and generates boundary data
% at all detector sites.
% The resulting projection vector is one-dimensional. It consists of two
% parts: the log of the modulation amplitude data, and the phase shift data.
% Each block contains sub-blocks for each source, and each sub-block
% consists of the data for each measurement site.
%
% The returned data vector is real of length 2*nqm, where nqm is the
% length of the permutation vector returned by toastDataLinkList, that
% is, unused source-detector combinations are removed from the result.

% Subsititute missing parameters with defaults
if nargin < 8
    solver = 'direct';
elseif nargin < 9
    tol = 1e-10;
end

% Calculate system matrix
phi = toastFields(hMesh, 0,  qvec, mua, mus, ref, omega, solver, tol, 'auto');

% Perform projection and remove unused measurements
lgamma = reshape (log(mvec.' * phi), [], 1);
lgamma = lgamma(hMesh.DataLinkList());

% Rearrange data in terms of log amplitude and phase shift blocks
if omega > 0
    proj = [real(lgamma);imag(lgamma)];
else
    proj = lgamma;
end
