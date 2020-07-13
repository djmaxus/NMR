%% startup
close all;
clc;
addpath('./DIF1D');
addpath('./POD');
addpath('./NMR');

%% set samples
% max parameter absolute value
% fid = fopen('./data/vals.mat');
% if fid ~= -1
%     load('./data/vals','vals','pMax');
%     fclose(fid);
% else
    pMax = 1;
    vals = setParSample(pMax);
% end
% save('./data/vals','vals','pMax');

fPar = figure;
fPar.WindowStyle = 'docked';
[~] = plotParSample(vals,[],fPar);
drawnow;

%% set other parmeters for dataset
NxMin = 5000;
NxMax = 5000;
NxCount = 1;
param.grids = round(logspace(log10(NxMin),log10(NxMax),NxCount));
param.Nt = 100;
param.D.A = 0.9;
param.D.P = 1;
param.D.S = 0.5e-1;
param.bc.first = true;
param.vals = vals;

%% numerical solution

% space grid
iNx = NxCount;
Nx = param.grids(iNx);

% BC kind
kind = 1;

% set boundary condition
bc = setBC(kind,pMax);

% set diffusion
D = setDiff(Nx,param.D);
%%
% solve FOM
[FOM] = solveDif(Nx,param.Nt,bc,D);

%% plot FOM
fFOM = figure;
fFOM.WindowStyle = 'docked';
fFOM.Name = 'FOM';
%%
[fFOM] = plotSol(FOM,fFOM,true);
drawnow;

%% PCA 2D-demo
% f2D = figure;
% f2D.WindowStyle = 'docked';
% [f2D] = pcaDemo(f2D);
% drawnow;

%% PCA ND-dataset
% % load data
% load 'data/misformax';
% f3D = figure;
% f3D.WindowStyle = 'docked';
% [f3D] = pcaDemoN(m,f3D);
% drawnow;

%% PCA N-D demo
% fPOD = figure;
% fPOD.WindowStyle = 'docked';
% [POD.demo] = makePOD(m(1:60,:),false,true,0.9,fPOD);
% drawnow;

%% train or load ROM-I data
% fid_1 = fopen('./data/train_1.mat');
% if fid_1 ~= -1
%     load('./data/train_1','train_1');
%     [train_1] = genRom1(param,train_1);
%     fclose(fid_1);
% else
    [train_1] = genRom1(param);
% end

% save('./data/train_1','train_1');

%% ROM-I train
fPar = figure;
fPar.WindowStyle = 'docked';
[fPar] = plotParSample(vals,bc.val,fPar);
drawnow;

%% ROM-I train
[POD.X] = truncRom1(train_1{iNx},true);

%% ROM-I test
[ROM_1] = solveDif(Nx,param.Nt,bc,D,1,0,POD);

%% ROM-I plot
fRom1 = figure;
fRom1.WindowStyle = 'docked';
[fRom1] = plotComp(FOM,ROM_1,D,fRom1);
drawnow;

%% Gappy-POD demo
% seminar1();

%% train or load ROM-II data
% fid_2 = fopen('./data/train_2.mat');
% if fid_2 ~= -1
%     load('./data/train_2','train_2');
%     [train_2] = genRom2(param,train_1,train_2);
%     fclose(fid_2);
% else
[train_2] = genRom2(param,train_1);
% end

% save('./data/train_2','train_2');

%% ROM-II train
[POD] = truncRom2(POD,train_2{iNx});

%% Gappy-POD train
[POD] = approxRom2(POD,train_2{iNx});

%% ROM-II test
[ROM_2] = solveDif(Nx,param.Nt,bc,D,'rom_stage',2,'encoder',POD);

%% ROM-II plot
fRom2 = figure;
fRom2.WindowStyle = 'docked';
%%
[fRom2] = plotComp(FOM,ROM_2,D,fRom2);
drawnow;

%% Data-driven approximation demo
fGrid = figure;
fGrid.WindowStyle = 'docked';
[parGrid,fGrid] = getParGrid(vals,bc.val,fGrid);
drawnow;
[nV,xV,in] = findSimplex(bc.val,vals);

%% ROM-III initial
time_dd = tic;
Yp = zeros(POD.X.N,length(nV));
Ydd = zeros(POD.X.N,param.Nt + 1);
for nt = 1:param.Nt+1
    for i = 1:length(nV)
        ns = nV(i);
        Yp(:,i) = train_2{iNx}.sol{ns}.Y(:,nt);
        Ydd(:,nt) = interpTriang(bc.val,vals(nV,:),Yp')';
    end
end
Xdd = forward(Ydd,POD.X);
time_dd = toc(time_dd);

ROMdd = ROM_2;
ROMdd.X = Xdd;
ROMdd.timelog.overall = time_dd;
ROMdd.name = 'ROM_{DD}';

%% DD plot
fRomDD = figure;
fRomDD.WindowStyle = 'docked';
[fRomDD] = plotComp(FOM,ROMdd,D,fRomDD);
drawnow;

%% ROM-III test
[ROM_3] = solveDif(Nx,param.Nt,bc,D,3,0,POD,Ydd);

%% ROM-II plot
fRom3 = figure;
fRom3.WindowStyle = 'docked';
[fRom3] = plotComp(FOM,ROM_3,D,fRom3);
drawnow;


