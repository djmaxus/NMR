function [sol] = solveDif(Nx,POD)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    nmr = ~isempty(POD);
    t = 0;
    T = 1.0;
    Nt = Nx;
    dt = T/Nt;
    dt_min = dt;
    x0 = zeros(Nx,1);
    X = [x0,x0];
    if nmr
        % Galerkin projection
        GP = POD.X.Ut;
        P = GP;
        % X decomposition
        x0 = POD.X.init;
        X = [x0,x0];
        Y = zeros(POD.X.N,2);
        f = figure;
        f.WindowState = 'maximized';
        plot(X(:,1)-x0);
    end
    kmax = 12;
    D = 0.5 * sin(2*2*pi*(1:Nx+1)'./(Nx+1)) + 1;
    D = D * 1e-2;
    if nmr
        f = figure;
        f.WindowState = 'maximized';
        plot(D);
        xlim([1,Nx+1]);
        ylim([0, max(D)*1.1]);
    end
    D = D*Nx^2;
    bc.first = false;
    bc.val = [1,1];
    crit_conv = 1e-3;
    kmin = 3;
    sol.t = t;
    sol.X = X(:,2);
    [~,sol.Rxl,sol.Rxr] = buildResidual(X,dt,D,bc,POD);
    sol.Rx = sol.Rxr - sol.Rxl;
    tic;
    while t < T - dt_min
        dt = min(dt, T-t);
        if nmr
            Y(:,2) = Y(:,1);     
        else
            X(:,2) = X(:,1);
        end
        for k = 1:kmax
            [Rt,Rxl,Rxr] = buildResidual(X,dt,D,bc,POD);
            Rx = Rxr - Rxl;
            R = Rt + Rx;
            if nmr
                R_conv = P'* R; 
            else
                R_conv = R;
            end
            if k == 1
                R0 = R_conv;
            end
            delta = norm(R_conv);
            if delta <= crit_conv*max(norm(R0),1)
                t = t+dt;
                if nmr
                    Y(:,1) = Y(:,2);
                else
                    X(:,1) = X(:,2);  
                    sol.Rx = [sol.Rx, Rx];
                    sol.Rxl = [sol.Rxl, Rxl];
                    sol.Rxr = [sol.Rxr, Rxr];
                end
                sol.t = [sol.t,t];
                sol.X = [sol.X,X(:,2)];
                
                if k < kmin
%                     dt =  dt * 2;
                end
                break;
            elseif k == kmax
                dt = dt / 4;
                break;
            end
            [Jt,Jx] = buildJacobian(X,dt,D,bc,POD);
            J = diag(Jt) + Jx;
            if nmr
                % LSPG
                P = J*GP;
                % projection
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
            error(['not converged, residual = ',num2str(delta)]);
            break;
        end
    end
    toc;
end

