%% grid and operators
Nx = 100; % inner points number
dx = 1/Nx; % spatial resolution
dx2 = dx^2;

Nxe = Nx + 2; % inner and boundary points number
Nw = Nx + 1; % links number

Dxe = -speye(Nw, Nxe); 
Dxe(:,2:end) = Dxe(:,2:end) + speye(Nw); % (du/dx) * dx

Dw = Dxe(2:end,2:end); % (dw/dx)*dx

Dx = Dxe(:,2:end-1); % (dw/dx)*dx over inner points

msXe = movsum((1:Nxe)',2,'Endpoints','discard'); % indexing for upstream approximations

Odt = [-1;1]; % (du/dt)*dt
Jdt = speye(Nx);

%% model

K = ones(Nw,1);   % non-uniformity
K = sparse(diag(K));

Gf = @(u) 1+u.^2; % non-linearity
Gdf = @(u) 2*u;

q = zeros(Nx,1);  % actuation

%% conditions
U0 = rand(Nxe,1); % initial/boundary conditions
U = repmat(U0,[1,2]);

%% dynamic timesteps
Nt = 1;
dt= 1/Nt;

%% calculation

Udx = Dxe * U(:,2); % dU/dx
sdU = sign(Udx); % upstream direction
sX = round(0.5*(msXe+sdU)); % upstream index

G = Gf(U(:,2)); % nonlinear measurement
GG = sparse(diag(G(sX)));

W = -K*GG*Udx; % nonlinear flow

Udt = U(2:end-1,:) * Odt; % time derivative

r = dx2*Udt + dt * Dw*W - dt*dx2*q; % inner resisual

dG = Gdf(U(:,2)); % nonlinear measurement derivative
dG = sparse(diag(dG));

Jg = dG(sX,2:end-1); % upstream approximation

Jw = -K * (GG*Dx + Jg.*Udx); % flow jacobian

J = dx2*Jdt + dt*Dw*Jw; % inner Jacobian

%% params
genK = @(Nx) speye(Nx+1);
calcG = @(u)[1+u.^2,2*u];
genQ = @(Nx) zeros(Nx,1);  % actuation
genBC = @(Nx) U0;

params = gen_params(Nx,genK,calcG,genQ,genBC);

%% calc
[res,jac] = gen_both(repmat(params.U0,[1,2]),1,params);



