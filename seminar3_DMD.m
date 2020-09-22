startup;
%%
clc;
disp('seminar3 started');
%% definitions
definitions;
disp('definitions set');
%% plots
close all;
f.FOM = new_figure('FOM');
f.PCA = new_figure('PCA');
f.ROM = new_figure('ROM');
f.DMD = new_figure('DMD');
f.DMDc = new_figure('DMDc');
f.test = new_figure('test');
%% problem setup
Nx = 2000;
Nt = 100;
bc.kind = 1;
D =  D_init(0.0,3,0.5e-1);
x0 = setInit(Nx,0,0);
%% FOM simulation
source_term = q_init(0.6,-1e+3);
param = setProblem(Nx,bc,source_term,D,x0);
FOM = solveDif(Nx,Nt,param);
plotSol(FOM,f.FOM); drawnow;
%% PCA
X = FOM.X;
pca.X = PCA(FOM.X,0.99,'init', zeros(size(X,1),1),'sval_deg',1);
pca.X.plotSummary(f.PCA);
%% ROM_I
ROM = solveDif(Nx,Nt,param,'rom_stage',1,'encoder',pca);
plotComp(FOM,ROM,param.D,f.ROM);
%% explore source term decompostion
dt = 1/Nt;
Q = param.q.val*dt;
clf(f.test);
plot(axes(f.test),pca.X.truncate(Q)); figure(f.test);
%% perform DMD in encoded space

X1 = FOM.X(:,1:end-1);
X2 = FOM.X(:,2:end);
Y1 = pca.X.encode(X1);
Y2 = pca.X.encode(X2);
Qy = pca.X.encode(Q);

A = X2 * pinv(X1,pca.X.npc());
B = Y2 * pinv(Y1,50);
DMD = FOM;
x = DMD.X(:,2);
y = pca.X.encode(x);
tic;
for i = 3:size(DMD.X,2)
    % advance 1 step
    x = A*x;
    y = B * y;
    % decode
    xx = pca.X.decode(y);
    % write down
    DMD.X(:,i) = xx;
end
toc;
plotComp(FOM,DMD,param.D,f.DMD);
%% DMDc with unknown actuation matrix B
approx = @(X)(svd(X,'econ'));
C = full(Q).*ones(size(X1));
Omega = [X1;C];
[Uo,So,Vo] = approx(Omega);
So = diag(So);
Uo1 = Uo(1:2000,:);
Uo2 = Uo(2001:end,:);
[U2,S2,V2] = approx(X2);
S2 = diag(S2);
A = X2 * Vo * diag(So.^(-1))*Uo1';
B = X2 * Vo * diag(So.^(-1))*Uo2';
DMDc = FOM;
% DMDc.X(:,3:end) = A*X1 + B*Q;
x = DMDc.X(:,2);
tic;
for i = 3:size(DMDc.X,2)
    % advance 1 step
    x = A*x + B * Q;
    % write down
    DMDc.X(:,i) = x;
end
toc;
plotComp(FOM,DMDc,param.D,f.DMDc);
%%
[Phi,Lambda] = eig(A);
[PhiB,LambdaB] = eig(B);
%% plots
clf(f.test);
semilogx(axes(f.test),abs(Lambda(:)),'-');
figure(f.test);
%% eDMD
%% from databook.com
% [Phi,Lambda,b] = DMD(X1,X2,10);







