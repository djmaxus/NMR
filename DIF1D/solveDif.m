function [sol] = solveDif(Nx,idx_mask,POD)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    nmr = length(idx_mask) < Nx;
    t = 0;
    T = 2;
    Nt = Nx;
    dt = T/Nt;
    dt_min = 1e-5*dt;
    x0 = ones(Nx,1);
    X = [x0,x0];
    if nmr
        % Galerkin projection
        GP = POD.X.Ut;
        P = GP;
        % X decomposition
        Y = backward(X,POD.X);
        X = forward(Y,POD.X);
    end
    kmax = 12;
    D = (sin(2*2*pi*(1:Nx+1)'./(Nx+1)) + 2);
    if nmr
        f = figure;
        f.WindowState = 'maximized';
        plot(D);
        xlim([1,Nx+1]);
        ylim([0, max(D)]);
    end
    D = D * 1e-2;
    D = D*Nx^2;
    bc = [-2,2];
    crit_conv = 1e-8;
    kmin = 4;
    sol.t = t;
    sol.X = X(:,2);
    [~,sol.Rx] = buildResidual(X,dt,D,bc,idx_mask,POD);
    while t < T - dt_min
        dt = min(dt, T-t);
        if nmr
            Y(:,2) = Y(:,1);     
        else
            X(:,2) = X(:,1);
        end
        for k = 1:kmax
            [Rt,Rx] = buildResidual(X,dt,D,bc,idx_mask,POD);
            R = Rt + Rx;
            if nmr
                R_conv = P'* R; 
            else
                R_conv = R;
            end
            if max(abs(R_conv)) <= crit_conv
                t = t+dt;
                if nmr
                    Y(:,1) = Y(:,2);
                else
                    X(:,1) = X(:,2);    
                end
                sol.t = [sol.t,t];
                sol.X = [sol.X,X(:,2)];
                sol.Rx = [sol.Rx, Rx];
                if k < kmin
                    dt = dt * 2;
                end
                break;
            elseif k == kmax
                dt = dt / 4;
                break;
            end
            [Jt,Jx] = buildJacobian(X,dt,D,idx_mask,POD);
            J = diag(Jt) + Jx;
            if nmr
                P = GP;
                % LSPG
                P = J*GP;
                % Y derivative
                J = P' * J;
                R = P' * R;
                % Y derivative
                J = J * POD.X.Ut;
            end
            dX = J \ (-R);
            is_weak = sum(isnan(dX)) + sum(isinf(dX));
            if is_weak > 0
                dX = 0;
            end
            if nmr
                Y(:,2) = Y(:,2) + dX;
                X = forward(Y,POD.X);
            else
                X(:,2) = X(:,2) + dX;
            end
        end
        if dt < dt_min
            error('not converged');
            break;
        end
    end
end

