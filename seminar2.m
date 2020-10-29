%% startup
close all;
clc;

addpath('./DIF1D');
addpath('./POD');
addpath('./NMR');

addpath('../POD');
addpath('../objective-pca');
definitions;
%% plots
close all;
f.par = new_figure('Pars');
f.FOM = new_figure('FOM');
f.pca2D = new_figure('PCA 2D');
f.pcaND = new_figure('PCA ND');
%% set samples
pMax = 1;
vals = setParSample(pMax);
plotParSample(f.par,vals,[]);
figure(f.par);


%% prepare grid sequence
NxMin = 2000;
NxMax = 2000;
NxCount = 1;
grids = round(logspace(log10(NxMin),log10(NxMax),NxCount));
param.vals = vals;
%% problem setup
% space grid
iNx = NxCount;
Nx = grids(iNx);
% time
Nt = 100;
bc.kind = 1;
bc.val = rand(1,2) * 2 * pMax - pMax;
plotParSample(f.par,vals,bc.val); figure(f.par);
D =  D_init(0.9,1,0.5e-1);
x0 = setInit(Nx,0,0,0);
source_term = q_init(0.6,0);
param = setProblem(Nx,bc,source_term,D,x0);
param.grids = grids;
%%
% solve FOM
[FOM] = solveDif(Nx,Nt,param,'train_mode',1);
plotSol(FOM,f.FOM); figure(f.FOM);
%% PCA 2D-demo
pcaDemo(f.pca2D);
figure(f.pca2D);

%% PCA ND-dataset
% load data
% rsgeng2D()

%% generate dataset
pcaDemoN(f.pcaND);
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
plotParSample(vals,bc.val,f.par);
figure(f.par);

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
