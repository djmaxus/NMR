%% grid and operators
Nx = 5;
dx = 1/Nx;
dx2 = dx^2;

Nxe = Nx + 2;
Nw = Nx + 1;
Xw = (1:Nw)';

Dxe = -speye(Nw, Nxe); 
Dxe(:,2:end) = Dxe(:,2:end) + speye(Nw); 
%Dxe = Dxe./dx;

Dw = Dxe(2:end,2:end);

Dx = Dxe(:,2:end-1);

Xe = (1:Nxe)';

msXe = movsum(Xe,2,'Endpoints','discard'); % for upstream approximations

Nt = 1;
dt= 1/Nt;

dtdx2 = dt * dx2;

Odt = [-1;1];
Jdt = speye(Nx);

%% model

K = ones(Nw,1);   % non-uniform
K = sparse(diag(K));
Gf = @(u) 1+u.^2; % non-linear
Gdf = @(u) 2*u;
q = zeros(Nx,1);  % enforced

%% conditions
U = rand(Nxe,1) % initial/boundary conditions
U = repmat(U,[1,2]);

%% calculation

Udx = Dxe * U(:,2); % dU/dx
sdU = sign(Udx); % upstream direction
sX = round(0.5*(msXe+sdU)); % upstream index

G = Gf(U(:,2)); % nonlinear measurement

GG = sparse(diag(G(sX)));

W = -K*GG*Udx; % nonlinear flow

Udt = U(2:end-1,:) * Odt; % time derivative

r = dx2*Udt + dt * Dw*W - dtdx2*q % inner resisual

dG = Gdf(U(:,2)); % nonlinear measurement derivative
dG = sparse(diag(dG));

Jg = dG(sX,2:end-1); % upstream approximation

Jw = -K * (GG*Dx + Jg.*Udx); % flow jacobian

J = dx2*Jdt + dt*Dw*Jw; % inner Jacobian

full(J)



