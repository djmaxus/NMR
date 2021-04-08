function [sol] = solveDif(Nx,Nt,param,varargin)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here

    % timelog
    sol.timelog.linsol = 0;
    sol.timelog.proj = 0;
    sol.timelog.pod = 0;
    sol.timelog.resbld = 0;
    sol.timelog.jacbld = 0;
    sol.timelog.kn = zeros(1,Nt);

    % Parameter save
    sol.D = param.D;
    sol.bc = param.bc;

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

    % paramter distribution

    % space discretization
    D.val = param.D.val*Nx^2;
    D.model = param.D.model;

    %% convergence parameters
    kmax = 12;
    kmin = 4;
    crit_rel = 1e-3;
    crit_abs = 1e-3;

    %% dataset initialization
    sol.t = zeros(1,round(1/dt_min)+1);
    sol.X = zeros(Nx,round(1/dt_min)+1);
    if is_encoded
        sol.Y = zeros(POD.X.npc(),round(1/dt_min)+1);
    end
    if train_mode
        sol.Rx = {};
        sol.Jnz = {};
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
        if train_mode
            sol.Rx = [sol.Rx,cell(1)];
            sol.Rx{end} = zeros(Nx+1,kmax);
            sol.Jnz = [sol.Jnz,cell(1)];
            sol.Jnz{end} = zeros(Nx+1,2*kmax);
        end
        for k = 1:kmax

            if and(is_encoded,~is_step_fom)
                % decoding
                t_pod = tic;
                X(:,2) = POD.X.decode(Y(:,2));
                t_pod = toc(t_pod);
                sol.timelog.pod = sol.timelog.pod + t_pod;
            end

            % Residual

            t_res_build = tic;
            [Rt,W] = buildResidual(X,dt,D,param.bc,PODr);
            t_res_build = toc(t_res_build);
            sol.timelog.resbld = sol.timelog.resbld + t_res_build;
            Rx = W(2:Nx+1)-W(1:Nx);
            R = Rt + Rx;
            % source term
            Rs = param.q.val;
            R = R + Rs;


            % Jacobian
            t_jac_build = tic;
            [Jt,Jf] = buildFlow(X,dt,D,param.bc,PODj);
            t_jac_build = toc(t_jac_build);
            sol.timelog.jacbld = sol.timelog.jacbld + t_jac_build;
            [Jx] = buildJacobian(Jf,Nx);
            J = Jx + Jt;


            if train_mode
                % collect residual snapshots
                sol.Rx{end}(:,k) = W;
                sol.Jnz{end}(:,2*k-1:2*k) = Jf.*[1,-1];
            end

            % Decoder Jacobian
            if and(is_encoded,~is_step_fom)
                t_projection = tic;
                J = J * POD.X.JD(Y(:,2));
                t_projection = toc(t_projection);
                sol.timelog.proj = sol.timelog.proj + t_projection;
            end

            % convergence
            delta = norm(J'*R);

            if k == 1
                R0 = delta;
            end
            is_converged_abs = delta*dt <= crit_abs;
            is_converged_rel = delta <= crit_rel * R0;
            if is_converged_abs
                is_converged = is_converged_abs;
            else
                is_converged = is_converged_rel;
            end

            if is_converged
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
                if k <= kmin
                    dt =  dt * 2;
                    dt = min(dt,dt_max);
                end
                % log
                sol.timelog.kn(nstep) = k;
                if train_mode
                    sol.Rx{end}(:,k+1:end) = [];
                    sol.Jnz{end}(:,2*(k+1)-1:end) = [];
                end
                break;
            elseif k == kmax
                %                 dt = dt / 4;
                break;
            end

            R = J'*R;
            J = J' * J;
            % linsolve
            t_lin_solve = tic;
            dX = (J) \ (-R);
            t_lin_solve = toc(t_lin_solve);
            sol.timelog.linsol = sol.timelog.linsol + t_lin_solve;

            is_weak = sum(isnan(dX)) + sum(isinf(dX));
            if is_weak > 0
                dX(:) = 0;
            end

            % update
            if and(is_encoded,~is_step_fom)
                Y(:,2) = Y(:,2) + dX;
            else
                X(:,2) = X(:,2) + dX;
            end
        end

        if is_converged
            waitbar(t/1,wb,[wb_msg,msg_add(is_step_fom)]);
            if and(is_step_fom,dt >= 1/Nt)
                is_step_fom = false;
            end
        end

        if and(~is_converged,k >= kmax)
            if train_mode
                sol.Rx(end) = [];
                sol.Jnz(end) = [];
            end
            if adaptive_fom
                if (~is_step_fom)
                    is_step_fom = true;
                    continue;
                end
            end
            if (dt < dt_min)
                warning(['Not converged',newline(),...
                    'Relative residual: ', num2str(delta/R0)]);
                break;
            end
            dt = dt/4;
        end

    end

    sol.t(nstep) = t;
    sol.t(nstep+1:end) = [];

    sol.X(:,nstep) = X(:,1);
    sol.X(:,nstep+1:end) = [];
    if is_encoded
        sol.Y(:,nstep) = Y(:,1);
        sol.Y(:,nstep+1:end) = [];
    end

    sol.timelog.bld = sol.timelog.resbld + sol.timelog.jacbld + sol.timelog.pod;
    sol.timelog.lin = sol.timelog.linsol + sol.timelog.proj;
    sol.timelog.overall = sol.timelog.lin + sol.timelog.bld;
    sol.timelog.iters.bld = sum(sol.timelog.kn,'all');
    sol.timelog.iters.lin = sum(sol.timelog.kn-1,'all');
    sol.timelog.lin_avg = sol.timelog.lin/sol.timelog.iters.lin;
    sol.timelog.bld_avg = sol.timelog.bld/sol.timelog.iters.bld;
    sol.timelog.iters.lin_avg = sol.timelog.iters.lin / length(sol.timelog.kn);
    sol.timelog.iters.bld_avg = sol.timelog.iters.bld / (length(sol.timelog.kn)-1);
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
