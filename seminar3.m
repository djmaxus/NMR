startup;
%%
clc;
disp('seminar3 started');
%% definitions
definitions;
disp('definitions set');
%% plots
close all;
f.comp = new_figure('Train/Test comparison');
f.train = trainTestPlots('train',new_figure);
f.test = trainTestPlots('test',new_figure);
f.remap = trainTestPlots('remap',new_figure);
%% problem setup
Nx = 2000;
Nt = 100;
bc.kind = 1;
D =  D_init(0.0,3,0.5e-1);
x0 = setInit(Nx,0,0);
%% train source term
q_train = q_init(0.5,-1e+3);
train.param = setProblem(Nx,bc,q_train,D,x0);
[train] = trainTestSim(Nx,Nt,train.param,f.train);
%% test source term
q_test = q_init(0.51,-1e+3);
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
y = zeros(4,Nt+1);
y(1,:) = 0;
y(2,:) = dg0;
y(3,:) = dg1;
y(4,:) = 0;
dg_spline = spline(x,y',1:Nx)';
dg_lin = zeros(size(g0));
for i = 1:x0
    dg_lin(i,:) = dg0 .* ((i - 1)/(x0-1))^1; 
end
x_mid = round(0.5*(x0+x1));
for i = x0:x1
    dg_lin(i,:) = (dg1 - 0).* ((i - x_mid)/(x1-x_mid))^(1) + 0; 
end
for i = x1:Nx
    dg_lin(i,:) = (dg1 - 0).* ((i - Nx)/(x1 - Nx))^1 + 0; 
end
alpha = 1.0;
dg = alpha * dg_lin + (1- alpha)*dg_spline;
time_mult = exp(-(1:Nt+1));
g1 = g0 + dg*(x1-x0);
dX = test.FOM.X - train.FOM.X;
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
function Y = arrayRemapper(X,remapper)
    Y = zeros(size(X));
    for i = 1:size(X,1)
        j = round(remapper(i));
        Y(i,:) = X(j,:);
    end
end
function y = polynomialRemapper(x,lim,r0,r1,n)
    if x < lim(1)
        error('Out of bounds');
    elseif x == lim(1)
        y = lim(1);
    elseif x < r1
        y = r0*((x - lim(1))/r1)^n + lim(1);
    elseif x == r1
        y = r0;
    elseif x < lim(2)
        y = (lim(2) - r0)* ((x - r1)/(lim(2) - r1))^(n) + r0;
    elseif x == lim(2)
        y = lim(2);
    else
        error('Out of bounds')
    end
end


