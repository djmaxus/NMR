function [params] = gen_params(Nx,genK,calcG,genQ,genBC)
%GEN_PARAMS Summary of this function goes here
%   Detailed explanation goes here

%% spatial derivatives

Nxe = Nx + 2; % inner and boundary points number
Nw  = Nx + 1; % links number

Dxe          =              - speye(Nw, Nxe);
Dxe(:,2:end) = Dxe(:,2:end) + speye(Nw);      % (du/dx) * dx

%% compilation
params = struct(...
... space
    'dx2', (1/Nx)^2                                , ...
    'Dxe', Dxe                                     , ... (du/dx) * dx
    'Dw' , Dxe(2:end,2:end)                        , ... (dw/dx)*dx
    'Dx' , Dxe(:,2:end-1)                          , ... (dw/dx)*dx over inner points
    'Xup', movsum((1:Nxe)',2,'Endpoints','discard'), ... indexing for upstream approximations
... time
    'Dt' , [-1;1]                                  , ... (du/dt)*dt
    'Jt' , speye(Nx)                               , ...
... problem
    'K'  , sparse(diag(genK(Nx)))                  , ...
    'G'  , calcG                                   , ...
    'q'  , genQ(Nx)                                , ...
    'U0' , genBC(Nx)                                 ... boundary and initial conditions
    );
end

