function [sol] = solveDif(T,Nt,dt_lim,problem,nlin_solver,varargin)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
    % Parameter save
    sol.problem = problem;
    
    %% ROM mode processing
    rom_stage = 0;
    adaptive_fom = false;
    is_step_fom = false;
    POD = [];
    Ydd = [];
    if length(varargin) >= 2
        while length(varargin) >=2
            varname = varargin{1};
            if strcmp(varname,'rom_stage')
                rom_stage = varargin{2};
            end
            if strcmp(varname,'encoder')
                POD = varargin{2};
            end
            if strcmp(varname,'dd_approx')
                Ydd = varargin{2};
            end
            if strcmp(varname,'adaptive_fom')
                adaptive_fom = varargin{2};
                is_step_fom = true;
            end
            for i = 1:2
                varargin(1) = [];
            end
        end
    end
    %% stage name
    switch rom_stage
        case 0
            sol.name = 'FOM';
        case 1
            sol.name = 'ROM_{I}';
        case 2
            sol.name = 'ROM_{II}';
        otherwise
            sol.name = 'ROM';
    end


    %% encoding status
    [is_encoded,is_gappy] = setRomStage(rom_stage);
    if is_encoded
        if isempty(POD)
            error('Variable subspace POD.X required');
        end
        if isempty(POD.X)
            error('Variable subspace POD.X required');
        end
    end
    %% residual POD
    if is_gappy
        PODr = POD.Rx;
        PODj = POD.Jnz;
    else
        PODr = [];
        PODj = [];
    end
    dd_app = and(rom_stage > 0,~isempty(Ydd));

    %
    sol.is_encoded = is_encoded;
    sol.is_gappy = is_gappy;
    sol.dd_app = dd_app;

    %% problem statement
    % time settings
    t = 0;
    dt = T/Nt;
    dt_min = dt_lim(1);
    dt_max = dt_lim(2);
    %% initial condition
    x0 = problem.x0;
    X = [x0,x0];

    %% dataset initialization
    sol.t = zeros(1,round(1/dt_min)+1);
    sol.X = zeros(length(problem.x0),round(1/dt_min)+1);
    sol.X(:,1) = x0;

    %% solution process
    wb_msg = [sol.name,': step '];
    wb = waitbar(t,wb_msg);
    nstep = 1;
    while t < T
        dt = min(dt, T-t);
        
        % initial approximation
        X(:,2) = X(:,1);

        if and(is_gappy,~is_step_fom)
            PODr = POD.Rx;
            PODj = POD.Jnz;
        else
            PODr = [];
            PODj = [];
        end
        
        Fun = @(x) problem.residual([X(:,1),x],dt,PODr);
        Jac = @(x) problem.jacobian([X(:,1),x],dt,PODj);
        
        % Nonlinear solver
        [argout_nlin] = nlin_solver.eval(X(:,2),Fun,Jac);
        sol.nlin{nstep} = argout_nlin;
        

        %% Convergency check
        
        if argout_nlin.is_conv
            % solution update
            X(:,1) = argout_nlin.result;
            X(:,2) = argout_nlin.result;
            
            % time update
            nstep = nstep + 1;
            t = t+dt;
            
            % timestep increase
            if argout_nlin.iters <= nlin_solver.kmin
                dt =  dt * 2;
                dt = min(dt,dt_max);
            end
            
            % step type update
            if and(is_step_fom,dt >= 1/Nt)
                is_step_fom = false;
            end
            
            % log
            sol.X(:,nstep) = X(:,2);
            sol.t(nstep) = t;

            % waitbar
            waitbar(t/1,wb,[wb_msg,msg_add(is_step_fom)]);
        elseif argout_nlin.iters >= kmax
            if adaptive_fom
                if (~is_step_fom)
                    is_step_fom = true;
                    continue;
                end
            end
            if (dt < dt_min)
                warning(['Not converged',newline(),...
                    'Relative residual: ', num2str(argout_nlin.error_gauss_relative)]);
                break;
            end
            dt = dt/4;
        end
    end

    sol.t(nstep+1:end) = [];
    sol.X(:,nstep+1:end) = [];

    delete(wb);
end

function [is_encoded,is_gappy] = setRomStage(rom_stage)
    is_encoded = false;
    is_gappy = false;

    if rom_stage >= 1
        is_encoded = true;
    end
    if rom_stage >= 2
        is_gappy = true;
    end

end

function [out] = msg_add(is_fom)
    if is_fom
        out = 'FOM';
    else
        out = 'DEF';
    end
end
