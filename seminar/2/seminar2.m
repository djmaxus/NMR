%% startup
clc;
addpath(genpath('../../POD'));
addpath(genpath('../'));
definitions;
%% plots
close all;

f.par = new_figure('Pars');
f.grid = new_figure('Grid');
f.FOM = new_figure('FOM');

f.pca2D = new_figure('PCA 2D');
f.pcaND_1 = new_figure('PCA ND 1');
f.pcaND_2 = new_figure('PCA ND 2');

f.pod.x = new_figure('POD X');
f.pod.r = new_figure('POD R');
f.pod.j = new_figure('POD J');

f.rom.i = new_figure('ROM-I');
f.rom.ii = new_figure('ROM-II');
f.rom.iii = new_figure('ROM-III');
f.rom.dd = new_figure('ROM-DD');
%% set samples
pMax = 1;
vals = setParSample(pMax);
bc.val = rand(1,2) * 2 * pMax - pMax;
%%
plotParSample(f.par,vals,bc.val);
figure(f.par);

%% prepare grid sequence
NxMin = 500;
NxMax = 500;
NxCount = 1;
grids = round(logspace(log10(NxMin),log10(NxMax),NxCount));
%% problem setup
% space grid
iNx = NxCount;
Nx = grids(iNx);
% time
Nt = 100;
bc.kind = 1;
D =  D_init(0.9,1,0.5e-1);
x0 = setInit(Nx,0.5,2,1);
source_term = q_init(0.6,1e+3);
param = setProblem(Nx,bc,source_term,D,x0);
param.grids = grids;
param.vals = vals;
param.Nx = Nx;
param.Nt = Nt;
%% solve FOM
[FOM] = solveDif(Nx,Nt,param,'train_mode',0);
plotSol(FOM,f.FOM); figure(f.FOM);
%% PCA 2D-demo
pcaDemo(f.pca2D);
figure(f.pca2D);

%% N-dimensional PCA demo
pcaDemoN(f.pcaND_1,f.pcaND_2);
figure(f.pcaND_1);
figure(f.pcaND_2);

%% train or load ROM-I data
[train.i] = genRom1(param);

%% ROM-I train
[POD.X] = truncRom1(train.i);
POD.X.plotSummary(f.pod.x);
figure(f.pod.x);

%% ROM-I test
[ROM.i] = solveDif(Nx,Nt,param,'rom_stage',1,'encoder',POD);
plotComp(FOM,ROM.i,param.D,f.rom.i);

%% train or load ROM-II data
[train.ii] = genRom2(param,POD);

%% ROM-II train
[POD,Rx,Jnz] = truncRom2(POD,train.ii);
POD.Rx.plotSummary(f.pod.r); figure(f.pod.r);
POD.Jnz.plotSummary(f.pod.j); figure(f.pod.j);
%% Gappy-POD train
[POD] = approxRom2(POD,Rx,Jnz);

%% ROM-II test
[ROM.ii] = solveDif(Nx,Nt,param,'rom_stage',2,'encoder',POD);
plotComp(FOM,ROM.ii,param.D,f.rom.ii);
figure(f.rom.ii);

%% Data-driven approximation demo
[parGrid] = getParGrid(vals,bc.val,f.grid);
figure(f.grid);
[nV,xV,in] = findSimplex(bc.val,vals);

%% ROM-III initial
time_dd = tic;
Yp = zeros(POD.X.npc(),length(nV));
Ydd = zeros(POD.X.npc(),param.Nt + 1);
for nt = 1:param.Nt+1
    for i = 1:length(nV)
        ns = nV(i);
        Yp(:,i) = train.ii{ns}.Y(:,nt);
        Ydd(:,nt) = interpTriang(bc.val,vals(nV,:),Yp')';
    end
end
Xdd = POD.X.decode(Ydd);
time_dd = toc(time_dd);


%% ROM-III test
[ROM.iii] = solveDif(Nx,Nt,param,'rom_stage',2,'encoder',POD,'dd_approx',Ydd);
plotComp(FOM,ROM.iii,param.D,f.rom.iii);
figure(f.rom.iii);

%% ROM-DD
ROM.dd = ROM.ii;
ROM.dd.X = Xdd;
ROM.dd.timelog.overall = time_dd;
ROM.dd.name = 'ROM_{DD}';

plotComp(FOM,ROM.dd,param.D,f.rom.dd);
%% seminar III DMD
seminar3;
%%
seminar_encoder;
