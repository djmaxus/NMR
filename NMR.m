addpath('./DIF1D');
addpath('./POD');
clear all;
close all;
clc;
Nx = 100;
Nt = Nx;
T = 2;
D = 1e-2 * ones(Nx+1,1).*Nx^2;
bc = [2,1];
kmax = 12;
%% perform FOM
dt = 1/Nt;
t = 0;
X0 = ones(Nx,1);
% X0 = sin(2*pi*(1:Nx)'./Nx);
X = [X0,X0];
sol.t(1) = t;
sol.X = X(:,1);
[~,sol.Rx] = buildResidual(X,dt,D,bc,1:Nx);
while t < T
    last_step = false;
    if dt > T-t
        last_step = true;
    end
    dt = min(dt, T-t);
    X(:,2) = X(:,1);
    for k = 1:kmax
        [Rt,Rx] = buildResidual(X,dt,D,bc,1:Nx);
        R = Rt + Rx;
        if max(abs(R)) <= 1e-8
            t = t+dt;
            X(:,1) = X(:,2);
            sol.t = [sol.t,t];
            sol.X = [sol.X,X(:,2)];
            sol.Rx = [sol.Rx, Rx];
            if k < 4
                dt = dt * 2;
            end
            break;
        elseif k == kmax
            dt = dt / 4;
            break;
        end
        [Jt,Jx] = buildJacobian(X,dt,D,1:Nx);
        J = diag(Jt) + Jx;
        dX = J \ (-R);
        is_weak = sum(isnan(dX)) + sum(isinf(dX));
        if is_weak
            dX = 0;
        end
        X(:,2) = X(:,2) + dX;
    end
    if and(dt < 1e-6/Nt,~last_step)
        error('not converged');
        break;
    end
end
FOM.sol = sol;
f = figure;
f.WindowState = 'maximized';
plot(FOM.sol.X(:,end));
%% preform model reduction
% X POD
POD.X = makePOD(sol.X,0.99);
% MSE.X = calcMSE(POD.X.U,POD.X.S,sol.X);
% semilogy(MSE.X);
% Residual POD
POD.Rx = makePOD(sol.Rx,0.99);
% MSE.Rx = calcMSE(POD.Rx.U,POD.Rx.S,sol.Rx);
% semilogy(MSE.X);
% greedy mask
mask = greedy(POD.Rx.U);
frac = 0.1;
Nm = round(frac*Nx);
Nm = max(Nm,POD.Rx.N);
mask = mask(:,1:Nm);
%oblique projector
O = oblique(mask,POD.Rx.Ut);
out = [];
for i = 1:size(mask,2)
    for j=1:size(mask,1)
        if mask(j,i) ==1
            out = [out;j];
            break;
        end
    end
end
idx_mask = out;
%% perform ROM
clear sol;
dt = 1/Nt;
t = 0;
X = [X0,X0];
Y = backward(X0,POD.X);
sol.t(1) = t;
sol.X = forward(backward(X0,POD.X),POD.X);
while t < T
    dt = min(dt, T-t);
    Y(:,2) = Y(:,1);
    for k = 1:kmax
        [Rt,Rx] = buildResidual(X,dt,D,bc,idx_mask);
        % residual oblique projection
        Rx1 = POD.Rx.mu + O * (Rx - POD.Rx.mu);
        delta = sum(abs(mask'*(Rx - Rx1)),'all');
        R = Rt + Rx1;
        Rr = POD.X.Ut'*R;
        if max(abs(Rr)) <= 1e-8
            t = t+dt;
            Y(:,1) = Y(:,2);
            sol.t = [sol.t,t];
            sol.X = [sol.X,X(:,2)];
            if k < 3
                dt = dt * 2;
            end
            break;
        elseif k == kmax
            dt = dt / 4;
            break;
        end
        [Jt,Jx] = buildJacobian(X,dt,D,idx_mask);
        % jacobian oblique projector
        delta = sum(mask'*abs(O*Jx-Jx),'all');
        J = diag(Jt) + O*Jx;
        % Galerkin
        P = POD.X.Ut;
        % Least squares Petrov – Galerkin projector
%         P = J * P;
        % projection
        J = P' * J;
        R = P' * R;
        % Y derivative
        J = J * POD.X.Ut;
        dY = J \ (-R);
        Y(:,2) = Y(:,2) + dY;
    end
    if dt < 1e-6/Nt
        error('not converged');
        break;
    end
end
ROM.sol = sol;
