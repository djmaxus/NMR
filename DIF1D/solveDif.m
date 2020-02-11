function [sol] = solveDif(Nx,POD)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    nmr = ~isempty(POD);
    t = 0;
    T = 1.0;
    dt = T/Nx;
    dt_min = dt*1e-2;
    x0 = zeros(Nx,1);
    X = [x0,x0];
    if nmr
        % Galerkin projection
        GP = POD.X.Ut;
        P = GP;
        % X decomposition
        X = forward(backward(X,POD.X),POD.X);
        f = figure;
        f.WindowState = 'maximized';
        plot(X(:,1)-x0);
    end
    kmax = 12;
    D = 0.5 * sin(2*2*pi*(1:Nx+1)'./(Nx+1)) + 1;
    D = D * 1e-1;
    if nmr
        f = figure;
        f.WindowState = 'maximized';
        plot(D);
        xlim([1,Nx+1]);
        ylim([0, max(D)*1.1]);
    end
    D = D*Nx^2;
    bc.first = true;
    bc.val = [1,-1];
    crit_conv = 1e-8;
    kmin = 4;
    sol.t = t;
    sol.X = X(:,2);
    [~,sol.W] = buildResidual(X,dt,D,bc,POD);
    tic;
    while t < T - dt_min
        dt = min(dt, T-t);
        X(:,2) = X(:,1);
        if nmr; Y = backward(X(:,2),POD.X); end
        for k = 1:kmax
            %% residual
            [Rt,W] = buildResidual(X,dt,D,bc,POD);
            R = Rt + W(2:Nx+1) - W(1:Nx);
            
            %% Jacobian        
            [Jt,Jx] = buildJacobian(X,dt,D,bc,POD);   
            J = Jx + + diag(Jt);
            if nmr
                % LSPG
%                 P = J*GP;
                % projection
                R = P' * R;
                J = P' * J; 
                % PCA derivative
                J = J * POD.X.Ut;
            end
            
            %% convergence
            R_conv = R;
            delta = norm(R_conv);
            
            if delta*dt <= crit_conv
                t = t+dt;
                X(:,1) = X(:,2);  
                sol.t = [sol.t,t];
                sol.X = [sol.X,X(:,2)];
                sol.W = [sol.W, W];
                if nmr
                    if k < kmin
                        dt =  dt * 2;
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
    toc;
end

