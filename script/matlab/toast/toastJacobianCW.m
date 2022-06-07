function J = toastJacobianCW(mesh,basis,varargin)
% Jacobian for DOT continuous wave problem for mua parameter
%
% Syntax: [J, proj] = toastJacobianCW (mesh, basis, qvec, mvec, mua, mus, ref, solver, tol, impl)
%         [J, proj] = toastJacobianCW (mesh, basis, dphi, aphi, proj)
%
% Parameters:
%         mesh [toastMesh object]:
%             mesh object
%         basis [toastBasis object]:
%             basis mapper object (or 0 to return Jacobian in mesh basis)
%         qvec [nlen x nq sparse matrix]:
%             array of source column vectors
%         mvec [nlen x nm sparse matrix]:
%             array of measurement operator column vectors
%         mua [nlen vector]:
%             nodal absorption coefficients
%         mus [nlen vector]:
%             nodal scattering coefficients
%         ref [nlen vector]:
%             nodal refractive index coefficients
%         solver [string]:
%             linear solver (DIRECT|CG|BICG|BICGSTAB|GMRES)
%         tol [scalar]:
%             linear solver tolerance (optional, iterative solvers only)
%         impl (string):
%             solver provider (AUTO|TOAST|MATLAB)
%         dphi: [nlen x nq real matrix]
%             real direct fields (nodal basis)
%         aphi: [nlen x nm real matrix]
%             real adjoint fields (nodal basis)
%         proj: [nqm vector]
%             real (linear) boundary projection data
%
% Return values:
%         J: [nqm x slen dense real matrix]:
%             Jacobian matrix
%
%         proj: [nqm dense real matrix]:
%             loggarithm of boundary projection data
%
% Notes:  Calculates the derivative of the logarithm of the CW amplitude
%         data with respect to the absorption coefficients of the forward
%         operator.
%
%         The returned matrix is of size nqm x slen, where nqm is the number
%         of measurements, and slen is the dimension of the reconstruction
%         basis.

if isobject(basis)
    hb = basis.handle;
else
    hb = 0;
end


if nargin==5
    % Compute from fields
    dphi = varargin{1};
    aphi = varargin{2};
    proj = varargin{3};
    J = toastmex(uint32(54),mesh.handle,hb,dphi,aphi,proj);
else
    % Compute fields
    qvec = varargin{1};
    mvec = varargin{2};
    mua = varargin{3};
    mus = varargin{4};
    ref = varargin{5};
    solver = 'direct';
    tol = 1e-8;
    impl = 'auto';
    if nargin >= 8
        solver = varargin{6};
        if nargin >= 9
            tol = varargin{7};
            if nargin >= 10
                impl = varargin{8};
            end
        end
    end

    % Compute fields in mesh basis
    phi = toastFields(mesh,0,[qvec mvec],mua,mus,ref,0,solver,tol,impl);
    dphi = phi(:, 1:size(qvec,2));
    aphi = phi(:, (size(qvec,2)+1) : end);
        
    % Build projection data, reduce by linklist
    proj = reshape (mvec.' * dphi, [], 1);
    proj = proj(mesh.DataLinkList());

    % Compute Jacobian
    J = toastmex(uint32(54),mesh.handle,hb,dphi,aphi,proj);

end

proj = log(proj)

end
