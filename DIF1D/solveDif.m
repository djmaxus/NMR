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
    
    %% rom stage
    if length(varargin) >= 1
        lstage = varargin{1};
    else
        lstage = 0;
    end
    %% train state
    if length(varargin) >= 2
        train = varargin{2};
    else
        train = 0;
    end
    
    %% POD data
    if length(varargin) >= 3
        POD = varargin{3};
    else
        POD = [];
    end
    
    %% Data-driven approximation
    if length(varargin) >= 4
        Ydd = varargin{4};
    end
    
    %% stage name
    switch lstage
        case 0
            sol.name = 'FOM';
        case 1
            sol.name = 'ROM_I';
        case 2
            sol.name = 'ROM_{II}';
        case 3
            sol.name = 'ROM_{III}';
        otherwise
            sol.name = 'ROM';
    end
    
    
    %% variable POD
    v_pod = lstage >= 1;
    if v_pod
        if isempty(POD)
            error('Variable subspace POD.X required');
        end
        if isempty(POD.X)
            error('Variable subspace POD.X required');
        end
    end
    %% residual POD
    r_pod = lstage >= 2;
    if r_pod
        PODr = POD.Rx;
        PODj = POD.Jnz;
    else
        PODr = [];
        PODj = [];
    end
    dd_app = lstage >= 3;
    
    %
    sol.v_pod = v_pod;
    sol.r_pod = r_pod;
    sol.dd_app = dd_app;
    
    %% problem statement
    % time settings
    t = 0;
    dt = 1/Nt;
    if train
        dt_min = min(dt,1/Nx);
        dt_max = dt;
    else
        dt_min = dt;
        dt_max = max(dt,10/Nx);
    end
    % initial condition
    if ~v_pod
        a = 1.0;
        p = 2;
        x0 = a * sin(p*2*pi*(1:Nx)'./Nx);
        X = [x0,x0];
    else
        Y = zeros(POD.X.N,2);
        X = zeros(Nx,2);
        GP = POD.X.U(:,1:POD.X.N);
    end
    
    % paramter distribution
    
    % space discretization
    D = D*Nx^2;
    
    % convergence parameters
    kmax = 12;
    kmin = 4;
    crit_conv = 1e-3;
    
    % dataset initialization
    sol.t = zeros(1,round(1/dt_min)+1);
    sol.X = zeros(Nx,round(1/dt_min)+1);
    if v_pod
        sol.Y = zeros(POD.X.N,round(1/dt_min)+1);
    end
    if train
        sol.Rx = [];
        sol.Jnz = [];
    end
    
    % solution process
    wb = waitbar(t,'Forward problem');
    R0 = 1;
    nstep = 1;
    while t < 1
        dt = min(dt, 1-t);
        
        if ~v_pod
            X(:,2) = X(:,1);
        else
            X(:,1) = forward(Y(:,1),POD.X);
            if dd_app
                Yi = Ydd(:,nstep+1);
                alpha = rand;
                Y(:,2) = alpha * Yi + (1 - alpha)*Y(:,1);
            else
                Y(:,2) = Y(:,1);
            end
        end

        sol.t(nstep) = t;
        sol.X(:,nstep) = X(:,1);
        
        if v_pod
            sol.Y(:,nstep) = Y(:,1);
        end
        
        for k = 1:kmax
            
            if v_pod
                % decoding
                t_pod = tic;
                X(:,2) = forward(Y(:,2),POD.X);
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
            
            
            if train
                % collect residual snapshots
                sol.Rx = [sol.Rx,W];
                sol.Jnz = [sol.Jnz,Jf];
            end
            
            if v_pod
                % POD derivative
                t_projection = tic;
                J = J * GP;
                t_projection = toc(t_projection);
                sol.timelog.proj = sol.timelog.proj + t_projection;
            end
            
            R = J' * R;
            
            % convergence
            delta = norm(R);
            if k == 1
                R0 = delta;
            end
            is_converged = delta <= crit_conv*max(R0,1);

            if is_converged
                sol.timelog.kn(nstep) = k;
                nstep = nstep + 1;
                t = t+dt;
                if v_pod
                    Y(:,1) = Y(:,2);
                else
                    X(:,1) = X(:,2);
                end
                % timestep increase
                if k <= kmin
                    dt =  dt * 2;
                    dt = min(dt,dt_max);
                end
                waitbar(t/1,wb);
                break;
            elseif k == kmax
                dt = dt / 4;
                break;
            end
            
            % linsolve
            t_lin_solve = tic;
            dX = (J'*J) \ (-R);
            t_lin_solve = toc(t_lin_solve);
            sol.timelog.linsol = sol.timelog.linsol + t_lin_solve;
            
            is_weak = sum(isnan(dX)) + sum(isinf(dX));
            if is_weak > 0
                dX(:) = 0;
            end
            
            % update
            if v_pod
                Y(:,2) = Y(:,2) + dX;
            else
                X(:,2) = X(:,2) + dX;
            end
        end
        if dt + 1e-10 < dt_min
            disp(['not converged, relative residual = ',num2str(delta/max(R0,1))]);
            break;
        end
    end
    
    sol.t(nstep) = t;
    sol.t(nstep+1:end) = [];
    
    sol.X(:,nstep) = X(:,1);
    sol.X(:,nstep+1:end) = [];
    if v_pod
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
