function [sol] = solveDif(Nx,Nt,bc,D,varargin)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    %% timelog
    sol.timelog.linsol = 0;
    sol.timelog.proj = 0;
    sol.timelog.pod = 0;
    sol.timelog.resbld = 0;
    sol.timelog.jacbld = 0;
    sol.timelog.kn = zeros(1,Nt);
    %% argument unpacking
    if length(varargin) >= 1
        lstage = varargin{1};
    else
        lstage = 0;
    end
    %% stage name
    switch lstage
        case 0
            sol.name = 'FOM';
        case 2
            sol.name = 'ROM_I';
        case 3
            sol.name = 'ROM_{II}';
            if length(varargin) >= 3
                sol.name = 'ROM_{III}';
            else
                sol.name = 'ROM_{II}';
            end
        otherwise
            sol.name = 'ROM';
    end
    %% Parameter save
    sol.D = D;
    sol.bc = bc;
    %%
    if length(varargin) >= 2
        POD = varargin{2};
    else
        POD = [];
    end
    if length(varargin) >= 3
        Ydd = varargin{3};
    end
    %% stage requirements check
    v_pod = lstage > 1;
    if lstage == 2
        if isempty(POD)
            error('Variable subspace POD.X required');
        end
        if isempty(POD.X)
            error('Variable subspace POD.X required');
        end
    end
    r_pod = lstage == 3;
    if r_pod
        PODr = POD.Rx;
        PODj = POD.Jnz;
    else
        PODr = [];
        PODj = [];
    end
    %% problem statement
    % time settings
    t = 0;
    dt = 1/Nt;
    dt_min = dt;
    dt_max = dt;
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
    
    %% convergence parameters
    kmax = 8;
    kmin = 3;
    crit_conv = 1e-3;
    %% dataset initialization
    sol.t = zeros(1,Nt+1);
    sol.X = zeros(Nx,Nt+1);
    if lstage == 2
        sol.Y = [];
        sol.Rx = [];
        sol.Jnz = [];
    end
    %% solution process
    wb = waitbar(t,'Forward problem');
    R0 = 1;
    nstep = 1;
    while t < 1
        dt = min(dt, 1-t);
        
        if ~v_pod
            X(:,2) = X(:,1);
        else
            X(:,1) = forward(Y(:,1),POD.X);
            Y(:,2) = Y(:,1);
            if length(varargin) >= 3
                Yi = Ydd(:,nstep+1);
                alpha = rand;
                Y(:,2) = alpha * Yi + (1 - alpha)*Y(:,1);
            end
        end
        %         if r_pod
        %             if t==0
        %                 PODr = [];
        %                 PODj = [];
        %             end
        %         end
        sol.t(nstep) = t;
        sol.X(:,nstep) = X(:,1);
        if lstage == 2
            sol.Y(:,nstep) = Y(:,1);
        end
        
        for k = 1:kmax
            %% residual
            if v_pod
                t_pod = tic;
                X(:,2) = forward(Y(:,2),POD.X);
                t_pod = toc(t_pod);
                sol.timelog.pod = sol.timelog.pod + t_pod;
            end
            t_res_build = tic;
            [Rt,W] = buildResidual(X,dt,D,bc,PODr);
            t_res_build = toc(t_res_build);
            sol.timelog.resbld = sol.timelog.resbld + t_res_build;
            Rx = W(2:Nx+1)-W(1:Nx);
            R = Rt + Rx;
            %% Jacobian
            t_jac_build = tic;
            [Jt,Jf] = buildFlow(X,dt,D,bc,PODj);
            t_jac_build = toc(t_jac_build);
            sol.timelog.jacbld = sol.timelog.jacbld + t_jac_build;
            [Jx] = buildJacobian(Jf,Nx);
            J = Jx + Jt;
            % collect residual snapshots
            if lstage == 2
                sol.Rx = [sol.Rx,W];
                sol.Jnz = [sol.Jnz,Jf];
            end
            % LSPG
            if v_pod
                t_projection = tic;
                J = J * GP;
                t_projection = toc(t_projection);
                sol.timelog.proj = sol.timelog.proj + t_projection;
            end
            R = J' * R;
            
            %% convergence
            delta = norm(R);
            if k == 1
                R0 = delta;
            end
            is_converged = delta <= crit_conv*max(R0,1);
            %             is_converged = or(is_converged,delta <= crit_conv);
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
            
            %% linsolve
            t_lin_solve = tic;
            dX = (J'*J) \ (-R);
            t_lin_solve = toc(t_lin_solve);
            sol.timelog.linsol = sol.timelog.linsol + t_lin_solve;
            
            is_weak = sum(isnan(dX)) + sum(isinf(dX));
            if is_weak > 0
                dX(:) = 0;
            end
            
            %% update
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
    sol.X(:,nstep) = X(:,1);
    if lstage == 2
        sol.Y(:,nstep) = Y(:,1);
    end
    
    sol.timelog.bld = sol.timelog.resbld + sol.timelog.jacbld + sol.timelog.pod;
    sol.timelog.lin = sol.timelog.linsol + sol.timelog.proj;
    sol.timelog.overall = sol.timelog.lin + ...
        sol.timelog.bld;
    sol.timelog.iters.bld = sum(sol.timelog.kn,'all');
    sol.timelog.iters.lin = sum(sol.timelog.kn-1,'all');
    sol.timelog.lin_avg = sol.timelog.lin/sol.timelog.iters.lin;
    sol.timelog.bld_avg = sol.timelog.bld/sol.timelog.iters.bld;
    sol.timelog.iters.lin_avg = sol.timelog.iters.lin / Nt;
    sol.timelog.iters.bld_avg = sol.timelog.iters.bld / Nt;
    delete(wb);
end
