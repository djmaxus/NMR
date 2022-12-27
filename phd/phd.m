clear;

x_domain = {0, 1}; % domain
x_domain_length = x_domain{2} - x_domain{1};

Nf = 1000;
dx = x_domain_length / Nf;

Nu_all = Nf + 1; % grid size for u
u0 = linspace(0,1,Nu_all)'; % initial condition
ul = 1; % left boundary condition
ur = 0; % right boundary condition
u0(1) = ul;
u0(end) = ur;

u = u0;

dt = 1;

kappa = linspace(1,1,Nf)';

%%

calc_linearization = @(u0,dt,u) pde_linearization(u0, dt, u, dx, ul, ur, kappa);

function [u] = ode_step(calc_linearization,u0,dt)
calc_linear_for_newton = @(u) calc_linearization(u0,dt,u);
u = newton_iterations(calc_linear_for_newton,u0,1e-3,10);
end

function [sol] = newton_iterations(calc_linear,initial,tolerance,max_iters)

% set initial approximation
current = initial;

% launch loop
for i = 0:max_iters
    [residual, jacobian] = calc_linear(current);
    
    % calculate errors
    err_abs = norm(residual) / sqrt(length(residual));
    
    % memorize initials
    if i == 0
        err_initial = err_abs;
    end
    
    % calculate relative errors
    err_rel = err_abs / err_initial;
    
    % Newton convergence
    is_conv_abs = err_abs <= tolerance;
    is_conv_rel = err_rel <= tolerance;
    
    is_conv_any = is_conv_abs || is_conv_rel;
    
    % is_conv_any
    if is_conv_any
        break;
    end
    
    % solve linear system
    delta = jacobian \ residual;
    current = current - delta;
end

if is_conv_any
    sol = current;
else
    sol = [];
end

end
