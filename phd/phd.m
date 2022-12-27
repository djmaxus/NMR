clear;

%% Problem & Discretization
x_domain = {0, 1}; % domain
x_domain_length = x_domain{2} - x_domain{1};

Nf = 2000;
dx = x_domain_length / Nf;

Nu_all = Nf + 1; % grid size for u
u0 = linspace(0,1,Nu_all)'; % initial condition
ul = 1; % left boundary condition
ur = 0; % right boundary condition
u0(1) = ul;
u0(end) = ur;

kappa = linspace(1,1,Nf)';

%%
calc_linearization = @(u0,dt,u) pde_linearization(u0, dt, u, dx, ul, ur, kappa);

[U,timesteps] = ode_solve(calc_linearization,u0,100,true);
