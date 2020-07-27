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
f.remap.demo = new_figure('Remap Demo');
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
q_test = q_init(0.6,-1e+3);
test.param = setProblem(Nx,bc,q_test,D,x0);
[test] = trainTestSim(Nx,Nt,test.param,f.test);
%% sompare train and test
plotComp(test.FOM,train.FOM,test.param.D,f.comp);
%% compute source positions
x0 = train.param.q.pos;
x1 =  test.param.q.pos;
%% check Green function symmetry
ds = train.FOM.X(x1,:)./1000 - test.FOM.X(x0,:)./1000;
%% compute Green function point
g0 = sqrt(train.FOM.X(x0,:));
g1 = train.FOM.X(x1,:)./g0;
%% swap points between x0 and x1
dX = zeros(size(train.FOM.X));
i_min = min(x0,x1);
i_max = max(x0,x1);
for i = i_min:i_max
    j = -(i-i_max) + i_min;
    dX(i,:) = train.FOM.X(j,:);
end
for i = i_max:Nx
    j = (1 - i_min)/(Nx - i_max)*(i - i_max)+i_min;
    j = round(j);
    dX(j,:) = train.FOM.X(i,:);
    dX(i,:) = train.FOM.X(j,:);
end
for i = 1:i_min
    j = (Nx - i_max)/(1 - i_min)*(i - i_min) + i_max;
    j = round(j);
    dX(j,:) = train.FOM.X(i,:);
    dX(i,:) = train.FOM.X(j,:);
end
remap.FOM.X = dX;
%% remap PCA
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


