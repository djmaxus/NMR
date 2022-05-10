function [sol] = solveDif(T,Nt,problem,nlin_solver,varargin)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
    %% function constants
    timestep_inc = 2;
    timestep_dec = 0.25;
    %% ROM mode processing
%     adaptive_fom = false;
%     is_step_fom = false;
%     if length(varargin) >= 2
%         while length(varargin) >=2
%             varname = varargin{1};
%             if strcmp(varname,'adaptive_fom')
%                 adaptive_fom = varargin{2};
%                 is_step_fom = true;
%             end
%             for i = 1:2
%                 varargin(1) = [];
%             end
%         end
%     end

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
