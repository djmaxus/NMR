%% collect data for DMD
X = FOM.X;
W = FOM.Rx{1}(:,1);
for i = 1:length(FOM.Rx)
    W = [W,FOM.Rx{i}(:,end)];
end
dt = 1/Nt;
% state variable
X1 = X(:,1:end-1);
X_dot = diff(X,1,2)./dt;
% nonlinear flows
Y1 = W(:,1:end-1);
Y_dot = diff(W,1,2)./dt;
%% DMD
[A, U, L, b] = DMD(X1,X_dot,10);
%% DMD test
X2 = zeros(size(X));
X2(:,1) = X(:,1);
for i = 2:size(X,2)
    X2(:,i) = X2(:,i-1) + dt * A*X2(:,i-1);
end
norm(X2 - X)
%% eDMD
[Ae,Ue,Le,be] = DMD([X1;Y1],[X_dot;Y_dot],20);
%% eDMD test
Z2 = zeros(size([X;W]));
Z2(:,1) = [X(:,1);W(:,1)];
for i = 2:size(Z2,2)
    Z2(:,i) = Z2(:,i-1) + dt * Ae*Z2(:,i-1);
end
X2e = Z2(1:size(X,1),:);
norm(X2e - X)
