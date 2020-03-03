function [sol] = solveDif(Nx,bc,varargin)
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
    %% problem statement
    % time settings
    t = 0;
    T = 1.0;
    Nt = 100;
    dt = T/Nt;
    dt_min = dt;
    dt_max = dt;
    % initial condition
    if ~v_pod
        a = 0.1;
        p = 2;
        x0 = a * sin(p*2*pi*(1:Nx)'./Nx);
        X = [x0,x0];
    else
        Y = zeros(POD.X.N,2);
        GP = POD.X.U(:,1:POD.X.N);
    end

    % paramter distribution
    aD = 0.5 * 1e-2;
    pD = 2;
    dD = 1;
    D = aD * sin(pD*2*pi*(1:Nx+1)'./(Nx+1)) + dD;
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
    crit_conv = 1e-3;
    %% dataset initialization
    sol.t = [];
    sol.X = [];
    if lstage == 2
        sol.Rx = [];
        sol.Jnz = [];
    end
    wb = waitbar(t/T,'Forward problem');
    timer = 0;
    while t < T
        dt = min(dt, T-t);
        
        if ~v_pod
            X(:,2) = X(:,1);
        else
            Y(:,2) = Y(:,1);
            X(:,1) = forward(Y(:,1),POD.X);
        end
        sol.t = [sol.t,t];
        sol.X = [sol.X,X(:,1)];
        for k = 1:kmax
            %% residual
            t_mat_build = tic;
            if v_pod
                X(:,2) = forward(Y(:,2),POD.X);
            end
            [Rt,Rx] = buildResidual(X,dt,D,bc,[]);
            R = Rt + Rx;
            %% Jacobian
            [Jt,Jx,Jnz] = buildJacobian(X,dt,D,bc,[]);
            J = Jx + Jt;
            t_mat_build = toc(t_mat_build);
            timer = timer + t_mat_build;
            % collect residual snapshots
            if lstage == 2
                sol.Rx = [sol.Rx,Rx];
                sol.Jnz = [sol.Jnz,Jnz];
            end
            % LSPG
            if v_pod
                t_projection = tic;
                J = J * GP;
                t_projection = toc(t_projection);
                timer = timer + t_projection;
            end
            
            %% convergence
            delta = norm(R)*dt;
            if k == 1
                R0 = delta;
            end
            
            if delta <= crit_conv*max(R0,1)
                t = t+dt;
                X(:,1) = X(:,2);
                sol.t = [sol.t,t];
                sol.X = [sol.X,X(:,2)];
%                 sol.Rx = [sol.Rx, Rx];
%                 sol.Jnz = [sol.Jnz, Jnz];
                % timestep increase
                if k <= kmin
                    dt =  dt * 2;
                    dt = min(dt,dt_max);
                end
                waitbar(t/T,wb);
                break;
            elseif k == kmax
                dt = dt / 4;
                break;
            end
            
            %% linsolve
            t_lin_solve = tic;
            [Qj,Rj] = qr(J'*J);
            dX = Rj \ (-Qj'*J'*R);
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
        if dt < dt_min
            error(['not converged, residual = ',num2str(delta)]);
            break;
        end
    end
    sol.timer = timer;
    delete(wb);
end
