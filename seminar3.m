%%
clc;
startup;
disp('seminar3 started');
%% definitions
definitions;
disp('definitions set');
%% plots
close all;
f.comp = new_figure('Train/Test comparison');
f.train = trainTestPlots('train',new_figure);
f.test = trainTestPlots('test',new_figure);
%% problem setup
Nx = 2000;
Nt = 100;
bc.kind = 2;
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
