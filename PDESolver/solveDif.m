function [sol] = solveDif(Nx,Nt,param,varargin)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here

    % log
%     sol.time.solve = 0;
%     sol.time.proj = 0;
%     sol.time.build = 0;
% 
%     sol.nlin = cell(1,1);
%     sol.iters = 0;
    
    % Parameter save
    sol.param = param;
    
    %% ROM mode processing
    rom_stage = 0;
    train_mode = false;
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
            if strcmp(varname,'train_mode')
                train_mode = varargin{2};
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
    dt = 1/Nt;

    if train_mode
        dt_min = dt;
        dt_max = dt;
    else
        dt_min = dt;
        dt_max = dt;
    end
    %% initial condition
    x0 = param.x0;
    X = [x0,x0];
    if is_encoded
        Y = POD.X.encode(X);
    end

    %% convergence parameters
    kmax = 20;
    kmin = 5;
    crit_rel = 1e-3;
    crit_abs = 1e-3;
    nlin_solver = @(x,Res,Jac,Proj)(newton_common(x,Res,Jac,kmax,crit_abs,crit_rel,Proj));

    %% dataset initialization
    sol.t = zeros(1,round(1/dt_min)+1);
    sol.X = zeros(Nx,round(1/dt_min)+1);
    sol.X(:,1) = x0;
    if is_encoded
        sol.Y = zeros(POD.X.npc(),round(1/dt_min)+1);
    end

    %% solution process
    wb_msg = [sol.name,': step '];
    wb = waitbar(t,wb_msg);
    R0 = 1;
    nstep = 1;
    while t < 1
        dt = min(dt, 1-t);
        if dt < 1e-10
            break;
        end

        if and(is_encoded,~is_step_fom)
            Y(:,2) = Y(:,1);
            if dd_app
                Yi = Ydd(:,nstep+1);
                alpha = rand;
                Y(:,2) = alpha * Yi + (1 - alpha)*Y(:,1);
            end
        else
            X(:,2) = X(:,1);
        end

        sol.t(nstep) = t;
        sol.X(:,nstep) = X(:,1);
        if is_encoded
            sol.Y(:,nstep) = Y(:,1);
        end
        if and(is_gappy,~is_step_fom)
            PODr = POD.Rx;
            PODj = POD.Jnz;
        else
            PODr = [];
            PODj = [];
        end

        if is_encoded
            Arg = @(y)([X(:,1),POD.X.decode(y)]);
            u1 = Y(:,2);
            Proj = @(x)(POD.X.JD(x));
        else
            Arg = @(x)([X(:,1),x]);
            u1 = X(:,2);
            Proj = @(x)(speye(length(x)));
        end
        Fun = @(x) buildResidual(Arg(x),dt,param.D,param,PODr);
        Jac = @(x) buildJacobian(Arg(x),dt,param.D,param,PODj);

        [argout_nlin] = nlin_solver(u1,Fun,Jac,Proj);
        
        sol.nlin{nstep} = argout_nlin;
        
%         sol.time.build = sol.time.build + argout_nlin.build;
%         sol.time.solve = sol.time.solve + argout_nlin.solve;
%         sol.time.proj  = sol.time.proj + argout_nlin.proj;

        %% Convergency check
        
        if argout_nlin.is_conv
            X = Arg(argout_nlin.result);
            
            if (is_encoded)
                Y(:,2) = argout_nlin.result;
            end
            % time update
            nstep = nstep + 1;
            t = t+dt;
            % solution update
            if and(is_encoded,~is_step_fom)
                Y(:,1) = Y(:,2);
                X(:,1) = POD.X.decode(Y(:,1));
            else
                X(:,1) = X(:,2);
            end
            if and(is_encoded,is_step_fom)
                Y = POD.X.encode(X);
            end
            % timestep increase
            if argout_nlin.iters <= kmin
                dt =  dt * 2;
                dt = min(dt,dt_max);
            end
            % step type update
            if and(is_step_fom,dt >= 1/Nt)
                is_step_fom = false;
            end
            
            % log
%             sol.iters = sol.iters + k;
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
    if is_encoded
        sol.Y(:,nstep) = Y(:,1);
        sol.Y(:,nstep+1:end) = [];
    end

%     sol.time.lin = sol.time.solve + sol.time.proj;
%     sol.time.all = sol.time.lin + sol.time.build;

%     sol.iters_avg = sol.iters / nstep;
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
