%% grid and operators
Nx = 5; % inner points number
dx = 1/Nx; % spatial resolution
dx2 = dx^2;

Nxe = Nx + 2; % inner and boundary points number
Nw = Nx + 1; % links number

Dxe = -speye(Nw, Nxe); 
Dxe(:,2:end) = Dxe(:,2:end) + speye(Nw); % (du/dx) * dx

Dw = Dxe(2:end,2:end); % (dw/dx)*dx

Dx = Dxe(:,2:end-1); % (dw/dx)*dx over inner points

msXe = movsum((1:Nxe)',2,'Endpoints','discard'); % indexing for upstream approximations

Nt = 1;
dt= 1/Nt;

dtdx2 = dt * dx2;

Odt = [-1;1]; % (du/dt)*dt
Jdt = speye(Nx);

%% model

K = ones(Nw,1);   % non-uniformity
K = sparse(diag(K));

Gf = @(u) 1+u.^2; % non-linearity
Gdf = @(u) 2*u;

q = zeros(Nx,1);  % actuation

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



