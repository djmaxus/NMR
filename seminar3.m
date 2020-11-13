%% collect data for DMD
X = FOM.X;
W = FOM.Rx{1}(:,1);
for i = 1:length(FOM.Rx)
    W = [W,FOM.Rx{i}(:,end)];
end
dt = 1/Nt;
% state variable
X1 = X(:,1:end-1);
X2 = X(:,2:end);
X_dot = diff(X,1,2)./dt;
% nonlinear flows
W1 = W(:,1:end-1);
W2 = W(:,2:end);
Y_dot = diff(W,1,2)./dt;
Ndmd = Nt;
%% DMD
[A, U, L, b] = DMD(X1,X_dot,Ndmd);
%% DMD test
Z2 = zeros(size(X));
Z2(:,1) = X(:,1);
B = eye(size(X,1))+dt*A;
for i = 2:size(X,2)
    Z2(:,i) = B * Z2(:,i-1);
end
norm(Z2 - X,'fro')
%% eDMD
[Ae,Ue,Le,be] = DMD([X1;W1],[X_dot;Y_dot],Ndmd);
%% eDMD test
Z2e = zeros(size([X;W]));
Z2e(:,1) = [X(:,1);W(:,1)];
for i = 2:size(X,2)
    Z2e(:,i) = Z2e(:,i-1) + dt * Ae*Z2e(:,i-1);
end
X2e = Z2e(1:size(X,1),:);
norm(X2e - X,'fro')
