clear;

x_domain = {0, 1}; % domain
x_domain_length = x_domain{2} - x_domain{1};

Nf = 10;
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

function [sol] = ode_solve(calc_linearization,u0,dt)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
    %% function constants
    timestep_inc = 2;
    timestep_dec = 0.25;

    %% Parameter save
    sol = problem;
    
    %% timestep evaluation
    t = 0;
    dt = T/Nt;
    
    %%  dynamic timestep
    if length(varargin) == 1
        dt_lim = varargin{1};
        dt_min = dt_lim(1);
        dt_max = dt_lim(2);
    else
        dt_min = dt;
        dt_max = dt;
    end
    %% initial condition
    x0 = problem.x0;
    X = [x0,x0];

    %% dataset initialization
    sol.t = zeros(1,round(T/dt_min)+1);
    sol.X = zeros(length(problem.x0),length(sol.t));
    sol.X(:,1) = x0;

    %% solution
    wb_msg = ['Solving ODE: ',problem.name];
    wb = waitbar(t, wb_msg, 'WindowStyle', 'modal');
    nstep = 1;
    while t < T
        dt = min([dt, T-t, dt_max]);
        if dt < 1e-16
            break;
        end
        
        X(:,2) = X(:,1); % initial approximation
        
        Fun = @(x) problem.residual([X(:,1),x],dt);
        Jac = @(x) problem.jacobian([X(:,1),x],dt);
        
        % Nonlinear solver
        init_approx = @()X(:,2);
        [argout_nlin] = nlin_solver(init_approx(),Fun,Jac);
        sol.nlin(nstep) = argout_nlin;
        
        if argout_nlin.is_conv % if converged
            % solution update
            X(:,1) = argout_nlin.result;
            X(:,2) = argout_nlin.result;
            
            % time update
            nstep = nstep + 1;
            t = t+dt;
            
            % timestep increase
            if argout_nlin.iters <= argout_nlin.param.iters_min
                dt =  dt * timestep_inc;
            end
            
            % step type update
%             if and(is_step_fom,dt >= 1/Nt)
%                 is_step_fom = false;
%             end
            
            % log
            sol.X(:,nstep) = X(:,2);
            sol.t(nstep) = t;

            % waitbar
            waitbar(t/T,wb,wb_msg);
        else % if not converged
%             if adaptive_fom
%                 if (~is_step_fom)
%                     is_step_fom = true;
%                     continue;
%                 end
%             end
            if (dt < dt_min)
                warning('Not converged');
                argout_nlin.err
                break;
            end
            dt = dt * timestep_dec;
        end
    end

    sol.t(nstep+1:end) = [];
    sol.X(:,nstep+1:end) = [];
    delete(wb);
end

function [u] = ode_step(calc_linearization,u0,dt)
calc_linear_for_newton = @(u) calc_linearization(u0,dt,u);
u = newton_iterations(calc_linear_for_newton,u0,1e-3,10);
end
