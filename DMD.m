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
f.DMD = new_figure('DMD');
f.Evo = new_figure('DMD Evo');
f.Compare = new_figure('Comparison');
% f.train = trainTestPlots('train',new_figure);
% f.test = trainTestPlots('test',new_figure);
% f.remap = trainTestPlots('remap',new_figure);
% f.green.dX = new_figure('Green difference');
% f.green.slice = new_figure('Green diff slices');
%% problem setup
Nx = 2000;
Nt = 100;
bc.kind = 1;
D =  D_init(0.0,3,0.5e-1);
x0 = setInit(Nx,0,0);
%% FOM simulation
source_term = q_init(0.5,-1e+3);
param = setProblem(Nx,bc,source_term,D,x0);
FOM = solveDif(Nx,Nt,param);
plotSol(FOM,f.FOM); drawnow;
%% perform Full DMD
X1 = FOM.X(:,2:end-1);
X2 = FOM.X(:,3:end);
[U,S,V] = svds(X1,10);
X1_inv = V * diag(diag(S).^(-1)) * U';
A = X2*X1_inv;
X2_DMD = A*X1;
FOM_DMD = FOM;
FOM_DMD.X(:,3:end) = X2_DMD;
plotSol(FOM_DMD,f.DMD); drawnow;
plotComp(FOM,FOM_DMD,param.D,f.Compare);
%% Apply Evolution
FOM_Evo = FOM;
x = FOM.X(:,2);
for i = 3:size(FOM.X,2)
    x = A*x;
    FOM_Evo.X(:,i) = x;
end
plotSol(FOM_Evo,f.Evo); drawnow;
plotComp(FOM,FOM_Evo,param.D,f.Compare);
%% Leading DMD modes
A1 = U'*X2*V*diag(diag(S).^(-1));
[W,D] = eig(A1);
UA = X2*V*diag(diag(S).^(-1))*W;
[W0,D0] = eigs(A,10);
%% test source term
q_test = q_init(0.6,-1e+3);
test.param = setProblem(Nx,bc,q_test,D,x0);
[test] = trainTestSim(Nx,Nt,test.param,f.test);
%% sompare train and test
plotComp(test.FOM,train.FOM,test.param.D,f.comp);
%% compute source positions
x0 = train.param.q.pos;
x1 =  test.param.q.pos;
g0 = train.FOM.X;
dg0 = (g0(x1,:) - g0(x0,:))/(x1-x0);
dg1 = -dg0;
x = [1,x0,x1,Nx];
dy = [zeros(1,Nt+1);dg0;dg1;zeros(1,Nt+1)];
dg = interp1(x,dy,1:Nx,'linear');
g1 = g0 + dg*(x1-x0);
y = g1(x,:);
z = interp1(x,y,1:Nx,'linear');
dz = (z - g0)./(x1-x0);
alpha = 1;
g1 = g0 + (alpha * dg + (1-alpha)*dz)*(x1-x0);
dX = test.FOM.X - train.FOM.X;
%% plot differences
clf(f.green.dX);
ax.green.dX = subplot(1,3,1,axes(f.green.dX));
ax.green.mismatch = subplot(1,3,2,axes(f.green.dX));
ax.green.dg = subplot(1,3,3,axes(f.green.dX));
imagesc(ax.green.dX,dX); colorbar(ax.green.dX);
imagesc(ax.green.mismatch,dg*(x1-x0)); colorbar(ax.green.mismatch);
imagesc(ax.green.dg,g1); colorbar(ax.green.dg);
figure(f.green.dX);
%% plot slices
clf(f.green.slice);
n = 2;
ax.green.slice = subplot(1,2,1,axes(f.green.slice));
ax.green.sol = subplot(1,2,2,axes(f.green.slice));

plot(ax.green.slice,dX(:,n));
hold(ax.green.slice,'on');
plot(ax.green.slice,dg(:,n)*(x1-x0));
hold(ax.green.slice,'on');
scatter(ax.green.slice,x',dg(x,n)*(x1-x0));
hold(ax.green.slice,'off');

plot(ax.green.sol,test.FOM.X(:,n));
hold(ax.green.sol,'on');
plot(ax.green.sol,train.FOM.X(:,n));
hold(ax.green.sol,'on');
plot(ax.green.sol,g1(:,n));
hold(ax.green.slice,'on');
scatter(ax.green.sol,x',g1(x,n));
hold(ax.green.sol,'off');

figure(f.green.slice);
%% continuity
dg0_left = g0(x0,:) - g0(x0-1,:);
dg0_right = g0(x0+1,:) - g0(x0,:);
dg_right = dg(x0+1,:) - dg(x0,:);
dg_left = dg_right + (dg0_right-dg0_left)./(x1-x0);
x0_left = x0-1;
%% remap PCA
remap.FOM.X = g1;
remap.pod.X = PCA(remap.FOM.X,0.99,'init',remap.FOM.X(:,1),'sval_deg',1);
remap.pod.X.plotSummary(f.remap.POD,'data',test.FOM.X);
drawnow;
%% remap ROM-I
[remap.ROM_1] = solveDif(Nx,Nt,test.param,'rom_stage',1,'encoder',remap.pod,'train_mode',true);
plotComp(test.FOM,remap.ROM_1,test.param.D,f.remap.ROM_1);
%%
remap.FOM_trunc = test.FOM;
remap.FOM_trunc.X = remap.pod.X.truncate(test.FOM.X);
test.FOM_trunc = test.FOM;
test.FOM_trunc.X = test.pod.X.truncate(test.FOM.X);
plotComp(test.FOM_trunc,remap.FOM_trunc,test.param.D,f.remap.FOM);
%% functions
function [data] = trainTestSim(Nx,Nt,param,fig)
    data.param = param;
    %% data FOM simulation
    data.FOM = solveDif(Nx,Nt,data.param);
    plotSol(data.FOM,fig.FOM);
    drawnow;
    %% data PCA
    data.pod.X = PCA(data.FOM.X,0.99,'init',data.FOM.X(:,1),'sval_deg',1);
    data.pod.X.plotSummary(fig.POD);
    drawnow;
    %% data ROM-I
    [data.ROM_1] = solveDif(Nx,Nt,data.param,'rom_stage',1,'encoder',data.pod,'train_mode',true);
    plotComp(data.FOM,data.ROM_1,data.param.D,fig.ROM_1);
end

function [fig] = trainTestPlots(name,new_figure)
    %% plots
    fig.FOM = new_figure(['FOM ',name]);
    fig.POD = new_figure(['PCA ',name]);
    fig.ROM_1 = new_figure(['ROM-I ',name]);
end


