function [sol] = solveDif(Nx,POD,bc,learn)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    nmr = ~isempty(POD);
    t = 0;
    T = 1.0;
    dt = T/Nx;
    dt_min = T*1e-4;
    dt_max = T * 1e-1;
    if ~learn
        dt = dt_max;
        dt_min = dt;
    end
    %     x0 = zeros(Nx,1);
    x0 = 0.1 * sin(2*2*pi*(1:Nx)'./Nx);
    X = [x0,x0];
    if nmr
        % Galerkin projection
        GP = POD.X.Ut;
        % X decomposition
        Y = backward(X,POD.X);
        X = forward(Y,POD.X);
    end
    
    D = 0.5 * sin(2*2*pi*(1:Nx+1)'./(Nx+1)) + 1;
    D = D * 1e-2;
    D = D*Nx^2;
    
    if nmr
        f = figure;
        f.WindowState = 'maximized';
        subplot(2,1,1);
        plot(X(:,2));
        title('Initial field');
        xlabel('x');
        ylabel('U0');
        subplot(2,1,2);
        plot(D);
        title('Diffusion coefficient');
        xlabel('x');
        ylabel('D');
        xlim([1,Nx+1]);
        ylim([0, max(D)*1.1]);
    end
    
    
    kmax = 8;
    crit_conv = 1e-4;
    kmin = 3;
    sol.t = t;
    sol.X = X(:,2);
    [~,sol.Rx] = buildResidual(X,dt,D,bc,POD);
    [~,~,sol.Jnz] = buildJacobian(X,dt,D,bc,[]);
    timer = tic;
    while t < T - dt_min
        dt = min(dt, T-t);
        X(:,2) = X(:,1);
        if nmr
            Y = backward(X(:,2),POD.X);
        end
        for k = 1:kmax
            %% residual
            [Rt,Rx] = buildResidual(X,dt,D,bc,POD);
            R = Rt + Rx;
            %% Jacobian
            [Jt,Jx,Jnz] = buildJacobian(X,dt,D,bc,POD);
            J = Jx + diag(Jt);
            if nmr
                % Galerkin
                P = GP;
                % Least-squares Petrov -- Galerkin
                P = J*P;
                % projection
                R = P' * R;
                J = P' * J;
                % PCA derivative
                J = J * POD.X.Ut;
            end
            %% convergence
            delta = norm(R);
            if k == 1
                R0 = delta;
            end
            
            if delta <= crit_conv*max(R0,1)
                t = t+dt;
                X(:,1) = X(:,2);
                sol.t = [sol.t,t];
                sol.X = [sol.X,X(:,2)];
                sol.Rx = [sol.Rx, Rx];
                sol.Jnz = [sol.Jnz, Jnz];
                if ~learn
                    if k <= kmin
                        dt =  dt * 2;
                        dt = min(dt,dt_max);
                    end
                end
                break;
            elseif k == kmax
                dt = dt / 4;
                break;
            end
            
            %% linsolve
            [q,r] = qr(J);
            dX = -r\(q'*R);
            
            is_weak = sum(isnan(dX)) + sum(isinf(dX));
            if is_weak > 0
                dX(:) = 0;
            end
            
            %% update
            if nmr
                Y = Y + dX;
                X(:,2) = forward(Y,POD.X);
            else
                X(:,2) = X(:,2) + dX;
            end
        end
        if dt < dt_min
            error(['not converged, residual = ',num2str(delta)]);
            break;
        end
    end
    timer = toc(timer);
    sol.timer = timer;
end
