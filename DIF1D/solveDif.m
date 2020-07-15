function [sol] = solveDif(Nx,Nt,bc,D,varargin)
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
    sol.D = D;
    sol.bc = bc;
    
    %% ROM mode processing
    rom_stage = 0;
    train_mode = false;
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
            for i = 1:2
                varargin(1) = [];
            end
        end
    end
    is_step_fom = false;
    %% stage name
    switch rom_stage
        case 0
            sol.name = 'FOM';
        case 1
            sol.name = 'ROM_I';
        case 2
            sol.name = 'ROM_{II}';
        otherwise
            sol.name = 'ROM';
    end
    
    
    %% encoding status
    [is_encoded,is_gappy] = setRomStage(rom_stage,is_step_fom);
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
    dt_min = dt*1e-1;
    dt_max = dt;
    %     if train_mode
    %         dt_min = min(dt,1/Nx);
    %         dt_min = dt;
    %     else
    %         dt_min = dt;
    %         dt_max = max(dt,10/Nx);
    %     end
    %% initial condition
    a = 1.0;
    p = 2;
    x0 = a * sin(p*2*pi*(1:Nx)'./Nx);
    X = [x0,x0];
    if is_encoded
        Y = zeros(POD.X.npc(),2);
        X = POD.X.decode(Y);
        GP = POD.X.u(:,1:POD.X.npc());
    end
    
    % paramter distribution
    
    % space discretization
    D = D*Nx^2;
    
    % convergence parameters
    kmax = 8;
    kmin = 4;
    crit_conv = 1e-3;
    
    % dataset initialization
    sol.t = zeros(1,round(1/dt_min)+1);
    sol.X = zeros(Nx,round(1/dt_min)+1);
    if is_encoded
        sol.Y = zeros(POD.X.npc(),round(1/dt_min)+1);
    end
    if train_mode
        sol.Rx = [];
        sol.Jnz = [];
    end
    
    % solution process
    wb = waitbar(t,'Forward problem');
    R0 = 1;
    nstep = 1;
    while t < 1
        dt = min(dt, 1-t);
        
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
            [Rt,W] = buildResidual(X,dt,D,bc,PODr);
            t_res_build = toc(t_res_build);
            sol.timelog.resbld = sol.timelog.resbld + t_res_build;
            Rx = W(2:Nx+1)-W(1:Nx);
            R = Rt + Rx;
            
            % Jacobian
            t_jac_build = tic;
            [Jt,Jf] = buildFlow(X,dt,D,bc,PODj);
            t_jac_build = toc(t_jac_build);
            sol.timelog.jacbld = sol.timelog.jacbld + t_jac_build;
            [Jx] = buildJacobian(Jf,Nx);
            J = Jx + Jt;
            
            
            if train_mode
                % collect residual snapshots
                sol.Rx = [sol.Rx,W];
                sol.Jnz = [sol.Jnz,Jf];
            end
            
            if and(is_encoded,~is_step_fom)
                % POD derivative
                t_projection = tic;
                J = J * GP;
                t_projection = toc(t_projection);
                sol.timelog.proj = sol.timelog.proj + t_projection;
            end
            
            % convergence
            delta = norm(J' * R);
            if k == 1
                R0 = delta;
            end
            is_converged = delta <= crit_conv*norm(R0);
            
            if is_converged
                sol.timelog.kn(nstep) = k;
                nstep = nstep + 1;
                t = t+dt;
                if and(is_encoded,~is_step_fom)
                    Y(:,1) = Y(:,2);
                    X(:,1) = POD.X.decode(Y(:,1));
                else
                    X(:,1) = X(:,2);
                end
                % timestep increase
                if k <= kmin
                    dt =  dt * 2;
                    dt = min(dt,dt_max);
                end
                if and(is_encoded,is_step_fom)
                    Y = backward(X,POD.X);
                end 
                waitbar(t/1,wb);
                break;
            elseif k == kmax
                dt = dt / 4;
                break;
            end
            
            % linsolve
%             [Qj,Rj] = qr(J);
            t_lin_solve = tic;
            dX = (J'*J) \ (-J' * R);
%             dX = Rj \ (-Qj' * R);
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
        if and(~is_step_fom,dt < 1/Nt)
            is_step_fom = true;
            disp('FOM steps: on');
        end
        if and(is_step_fom,dt >= 1/Nt)
            is_step_fom = false;
            disp('FOM steps: off');
        end
        if dt + 1e-10 < dt_min
            if(is_step_fom)
                disp(['not converged, relative residual = ',num2str(delta/norm(R0))]);
                break;
            end
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
    sol.timelog.iters.lin_avg = sol.timelog.iters.lin / Nt;
    sol.timelog.iters.bld_avg = sol.timelog.iters.bld / Nt;
    delete(wb);
end

function [is_encoded,is_gappy] = setRomStage(rom_stage,is_step_fom)
    is_encoded = false;
    is_gappy = false;
    if ~is_step_fom
        if rom_stage >= 1
            is_encoded = true;
        end
        if rom_stage >= 2
            is_gappy = true;
        end
    end
end
