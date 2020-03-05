function [sol] = solveDif(Nx,Nt,bc,D,varargin)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    %% argument unpacking
    if length(varargin) >= 1
        lstage = varargin{1};
    else
        lstage = 0;
    end
    if length(varargin) >= 2
        POD = varargin{2};
    else
        POD = [];
    end
    if length(varargin) >= 5
        vals = varargin{3};
        grid = varargin{4};
        val = varargin{5};
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
    
    %     if v_pod
    %         f = figure;
    %         f.WindowState = 'maximized';
    %         subplot(2,1,1);
    %         plot(X(:,2));
    %         title('Initial field');
    %         xlabel('x');
    %         ylabel('U0');
    %         subplot(2,1,2);
    %         plot(D);
    %         title('Diffusion coefficient');
    %         xlabel('x');
    %         ylabel('D');
    %         xlim([1,Nx+1]);
    %         ylim([0, max(D)*1.1]);
    %     end
    
    %% convergence parameters
    kmax = 8;
    kmin = 3;
    crit_conv = 1e-4;
    %% dataset initialization
    sol.t = [];
    sol.X = [];
    if lstage == 2
        sol.Rx = [];
        sol.Jnz = [];
    end
    wb = waitbar(t,'Forward problem');
    timer = 0;
    R0 = 1;
    k_all = 0;
    while t < 1
        dt = min(dt, 1-t);
        
        if ~v_pod
            X(:,2) = X(:,1);
        else
            X(:,1) = forward(Y(:,1),POD.X);
            Y(:,2) = Y(:,1);
            if length(varargin) >= 5
                t_ind = find(POD.X.data{1}.t == t);
                Xp = zeros(Nx,length(POD.X.data));
                for j = 1:length(POD.X.data)
                    Xp(:,j) = POD.X.data{j}.X(:,t_ind+1);
                end
                Yp = backward(Xp - X(:,1),POD.X);
                Y(:,2) = Y(:,2) + interpTriang(val,vals,grid,Yp')';
            end
        end
%         if r_pod
%             if t==0
%                 PODr = [];
%                 PODj = [];
%             end
%         end
        sol.t = [sol.t,t];
        sol.X = [sol.X,X(:,1)];
        
        for k = 1:kmax
            %% residual
            if v_pod
                X(:,2) = forward(Y(:,2),POD.X);
            end
            t_mat_build = tic;
            [Rt,Rx,W] = buildResidual(X,dt,D,bc,PODr);
            R = Rt + Rx;
            %% Jacobian
            [Jt,Jx,Jnz] = buildJacobian(X,dt,D,bc,PODj);
            J = Jx + Jt;
            t_mat_build = toc(t_mat_build);
            timer = timer + t_mat_build;
            % collect residual snapshots
            if lstage == 2
                sol.Rx = [sol.Rx,W];
                sol.Jnz = [sol.Jnz,Jnz];
            end
            % LSPG
            if v_pod
                t_projection = tic;
                J = J * GP;
                t_projection = toc(t_projection);
                timer = timer + t_projection;
            end
            R = J' * R;
            
            %% convergence
            delta = norm(R)*dt;
            if k == 1
                R0 = delta;
            end
            is_converged = delta <= crit_conv*max(R0,1);
            %             is_converged = or(is_converged,delta <= crit_conv);
            if is_converged
                k_all = k_all + k;
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
            [Qj,Rj] = qr(J'*J);
            dX = Rj \ (-Qj'* R);
            t_lin_solve = toc(t_lin_solve);
            timer = timer + t_lin_solve;
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
            error(['not converged, residual = ',num2str(delta)]);
        end
    end
    sol.t = [sol.t,t];
    sol.X = [sol.X,X(:,1)];
    sol.k_avg = k_all/(length(POD.X.data{1}.t)-1);
    sol.timer = timer;
    delete(wb);
end
