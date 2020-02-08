addpath('./DIF1D');
startMRSTmodules;
Nx = 10;
Nt = 10;
D = ones(Nx-1,1);
bc = [1,1];
%% perform FOM
dt = 1/Nt;
t = 0;
X = ones(Nx,2);
sol.t(1) = t;
sol.X = X(:,2);
[~,sol.Rx] = buildResidual(X,dt,D,bc);
while t < 1
    X(:,2) = X(:,1);
    for k = 1:kmax
        [Rt,Rx] = buildResidual(X,dt,D,bc);
        R = Rt + Rx;
        if max(abs(R)) <= epsilon
            t = t+dt;
            X(:,1) = X(:,2);
            sol.t = [sol.t,t];
            sol.X = [sol.X,X(:,2)];
            sol.R = [sol.Rx, Rx];
            break;
        elseif k == kmax
            dt = dt / 4;
            break;
        end
        [Jt,Jx] = buildJacobian(X,dt,D);
        J = diag(Jt) + Jx;
        dX = J \ (-R);
        X(:,2) = X(:,2) + dX;
    end
    if dt < 1/Nt
        break;
    end
end
FOM.sol = sol;

%% preform model reduction
% X POD
POD.X = makePOD(sol.X,0.9);
MSE.X = calcMSE(POD.X.U,POD.X.S,sol.X);
semilogy(MSE.X);
% Residual POD
POD.Rx = makePOD(sol.Rx,0.9);
MSE.Rx = calcMSE(POD.Rx.U,POD.Rx.S,sol.Rx);
semilogy(MSE.X);
% greedy mask
mask = greedy(POD.Rx.U);
%oblique projector
O = oblique(maskMatrix(mask),POD.Rx.U);

%% perform ROM
delete sol;
dt = 1/Nt;
t = 0;
X = ones(Nx,2);
Y = POD.X.U' * X;
sol.t(1) = t;
sol.X = X(:,2);
while t < 1
    Y(:,2) = Y(:,1);
    for k = 1:kmax
        X = POD.X.U*Y;
        [Rt,Rx] = buildResidual(X,dt,D,bc,mask);
        % residual oblique projection
        R = Rt + O * Rx;
        if max(abs(R)) <= epsilon
            t = t+dt;
            Y(:,1) = Y(:,2);
            sol.t = [sol.t,t];
            sol.X = [sol.X,X(:,2)];
            break;
        elseif k == kmax
            dt = dt / 4;
            break;
        end
        [Jt,Jx] = buildJacobian(X,dt,D,mask);
        % jacobian oblique projector
        J = diag(Jt) + O * Jx;
        % Galerkin
        P = POD.X.U;
        % Least squares Petrov – Galerkin projector
        P = J * P;
        % projection
        J = P' * J;
        R = P' * R;
        % Y derivative
        J = J * POD.X.U;
        dY = J \ (-R);
        Y(:,2) = Y(:,2) + dY;
    end
    if dt < 1/Nt
        break;
    end
end
ROM.sol = sol;
