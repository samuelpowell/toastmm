% toastFields  - Calculate complex photon density fields
%
% Syntax: phi = toastFields(mesh,basis,qvec,mua,mus,ref,freq,method,tol,impl)
%
% Parameters:
%         mesh (toastMesh instance):
%             mesh object
%         basis (toastBasis instance):
%             basis object (set to 0 to return fields in mesh basis)
%         qvec (complex sparse matrix n x nq):
%             matrix of nq source vectors
%         mua (real array n):
%             nodal absorption coefficients [1/mm]
%         mus (real array n):
%             nodal scattering coefficients [1/mm]
%         ref (real array n):
%             nodal refractive index values
%         freq (scalar):
%             modulation frequency (set 0 for a CW solution) [MHz]
%         method (string):
%             solver method (DIRECT|CG|BICGSTAB|GMRES)
%         tol (scalar):
%             solver tolerance
%         impl (string):
%             solver provider (AUTO|TOAST|MATLAB)
%
% Return values:
%         phi (real or complex matrix slen x nq or n x nq):
%             photon density fields for all nq sources
%
% Notes:  If the basis parameter is set to 0, the fields are returned in
%         the mesh basis. Otherwise, they are returned in the solution
%         basis. The returned matrix contains the field for source i in
%         column i.
%
%         To compute adjoint fields, pass the matrix of measurement
%         vectors (mvec) instead of qvec.

function [phi] = toastFields(mesh,basis,qvec,mua,mus,ref,freq,method,tol,impl)

    if nargin < 10
        impl = auto;
        if nargin < 9
            tol = 1e-10;
            if nargin < 8
                method = 'direct';
            end
        end
    end
    
    method = lower(method);
    impl = lower(impl);
    
    % If implementation is automatic, we use MATLAB direct solvers, and TOAST iterative
    % solvers to exploit multithreading in each case
    switch impl

    case 'auto'
        if strcmp(method, 'direct')
            lib = false;
        else
            lib = true;
        end
    case 'toast'
        lib = true;
    case 'matlab'
        lib = false;
    otherwise
        error('Implementation option unknown');
    end

    %
    % Compute using library
    %
    if lib
        if basis == 0
            bhandle = 0;
        else
            bhandle = basis.handle;
        end   
        phi = toastmex(uint32(51),mesh.handle,bhandle,qvec,mua,mus,ref,freq,method,tol);
        return;
    end
    
    %
    % Compute using MATLAB solvers
    %
    nQ   = size(qvec, 2);
    n    = mesh.NodeCount;
    maxit = 100;
    gflag = 0;
    
    % Get the system matrix
    if freq == 0
        cw = true;
        qvec = real(qvec);
        mvec = real(mvec);
        S = dotSysmat (mesh,mua,mus,ref);
    else
        cw = false;
        S = dotSysmat (mesh,mua,mus,ref,'freq',freq);
    end
    
    if strcmp(method,'direct')
        % Direct solve, let mldivide decide
        dphi = S\qvec;
    else    
        
        % Preallocate for solution and create preconditioner
        if cw
            dphi = zeros(n, nQ);
            L = ichol(S);
        else
            dphi = complex(zeros(n, nQ));
            [L,U] = ilu(S);
        end
        
        switch lower(method)
    
        case 'cg'
            if ~cw
                error('Unable to solve complex system with conjugate gradient method');
            end
            for i = 1:nQ
                [dphi(:,i), flag] = pcg(S,full(qvec(:,i)),tol,maxit,L,L');
                gflag = gflag + flag;
            end
            
        case 'bicgstab'
            if cw
                for i = 1:nQ
                    [dphi(:,i), flag] = bicgstab(S,full(qvec(:,i)),tol,maxit,L');
                    gflag = gflag + flag;
                end
            else
                for i = 1:nQ
                    [dphi(:,i), flag] = bicgstab(S,full(qvec(:,i)),tol,maxit,L,U);
                    gflag = gflag + flag;
                end
            end
    
        case 'gmres'        
            if cw
                for i = 1:nQ
                    [dphi(:,i), flag] = gmres(S,full(qvec(:,i)),[],tol,maxit,L');
                    gflag = gflag + flag;
                end
            else
                for i = 1:nQ
                    [dphi(:,i), flag] = gmres(S,full(qvec(:,i)),[],tol,maxit,L,U);
                    gflag = gflag + flag;
                end
            end
    
        otherwise
            error('Unknown solver type');
        end
    end

    if gflag ~= 0
        error('Solver did not converge');
    end
    
    if basis ~= 0
        if cw
            phi = zeros(raster.slen, nQ);
        else
            phi = complex(zeros(raster.slen, nQ));
        end
    
        for i = 1:nQ
            phi(:,i) = basis.Map('M->S', sphi(:,i));
        end
    else
        phi = dphi;
    end
    