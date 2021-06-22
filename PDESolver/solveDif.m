function [sol] = solveDif(T,Nt,problem,nlin_solver,varargin)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
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
    sol.problem = problem;
    
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
    sol.X = zeros(length(problem.x0),round(T/dt_min)+1);
    sol.X(:,1) = x0;

    %% solution process
    wb_msg = ['Solving ODE: ',problem.name];
    wb = waitbar(t,wb_msg,'WindowStyle','modal');
    nstep = 1;
    while t < T
        dt = min(dt, T-t);
        if dt < 1e-16
            break;
        end
        % initial approximation
        X(:,2) = X(:,1);
        
        Fun = @(x) problem.residual([X(:,1),x],dt);
        Jac = @(x) problem.jacobian([X(:,1),x],dt);
        
        % Nonlinear solver
        [argout_nlin] = nlin_solver(X(:,2),Fun,Jac);
        sol.nlin(nstep) = argout_nlin;
        

        %% Convergency check
        
        if argout_nlin.is_conv
            % solution update
            X(:,1) = argout_nlin.result;
            X(:,2) = argout_nlin.result;
            
            % time update
            nstep = nstep + 1;
            t = t+dt;
            
            % timestep increase
            if argout_nlin.iters <= argout_nlin.param.iters_min
                dt =  dt * 2;
                dt = min(dt,dt_max);
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
        elseif argout_nlin.iters >= argout_nlin.param.iters
%             if adaptive_fom
%                 if (~is_step_fom)
%                     is_step_fom = true;
%                     continue;
%                 end
%             end
            if (dt < dt_min)
                warning(['Not converged',newline(),...
                    'Relative residual: ', num2str(argout_nlin.err.rel_g)]);
                break;
            end
            dt = dt/4;
        end
    end

    sol.t(nstep+1:end) = [];
    sol.X(:,nstep+1:end) = [];
    delete(wb);
end

% function [out] = msg_add(is_fom)
%     if is_fom
%         out = 'FOM';
%     else
%         out = 'DEF';
%     end
% end
