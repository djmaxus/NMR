startup;
%% plots
close all;
% plot FOM
fFOM = figure;
fFOM.WindowStyle = 'docked';
fFOM.Name = 'FOM';
% ROM-I figure
f_rom_1 = figure;
f_rom_1.WindowStyle = 'docked';
f_rom_1.Name = 'PCA ROM-I';
% ROM-I plot
f_comp_1 = figure;
f_comp_1.WindowStyle = 'docked';
% ROM-I plot
f_comp_vac = figure;
f_comp_vac.WindowStyle = 'docked';
% ROM-II figure
f_rom_2_r = figure;
f_rom_2_r.WindowStyle = 'docked';
f_rom_2_r.Name = 'PCA ROM-II R';
f_rom_2_j = figure;
f_rom_2_j.WindowStyle = 'docked';
f_rom_2_j.Name = 'PCA ROM-II J';
% ROM-II plot
f_sol_2 = figure;
f_sol_2.WindowStyle = 'docked';
%% set samples
pMax = 1;
vals = [];
%%
% vals = setParSample(pMax);
% [~] = plotParSample(vals,[],fPar);
% drawnow;

%% set other parmeters for dataset
NxMin = 1000;
NxMax = NxMin;
NxCount = 1;
param.grids = round(logspace(log10(NxMin),log10(NxMax),NxCount));
param.Nt = 100;
param.D.A = 0.9;
param.D.P = 1;
param.D.S = 0.5e-1;
param.bc.first = false;
param.vals = vals;

%% numerical solution

% space grid
iNx = NxCount;
Nx = param.grids(iNx);

% BC kind
kind = 2;

% set boundary condition
bc = setBC(kind,pMax);
vals = bc.val;
param.vals = vals;
param.bc = bc;
% set diffusion
D = setDiff(Nx,param.D);
%%
% solve FOM
[FOM] = solveDif(Nx,param.Nt,bc,D);
plotSol(FOM,fFOM);

%% train or load ROM-I data
[train_1] = genRom1(param);

%% ROM-I train
[POD.X] = truncRom1(train_1{iNx});
POD.X.plotSummary(f_rom_1,'layout','built-in'); drawnow;
                                    
%% ROM-I test
[ROM_1] = solveDif(Nx,param.Nt,bc,D,'rom_stage',1,'encoder',POD);
plotComp(FOM,ROM_1,D,f_comp_1);

%% VAC
vac.X = VAC(data,POD.X.npc());
%% VAC test
[vac_sol] = solveDif(Nx,param.Nt,bc,D,'rom_stage',1,'encoder',vac);
plotComp(FOM,vac_sol,D,f_comp_vac);
%% train or load ROM-II data
train_2 = train_1;...genRom2(param,train_1);

%% ROM-II train
[POD] = truncRom2(POD,train_2{iNx});
%% Gappy-POD train
[POD] = approxRom2(POD,train_2{iNx});

%% ROM-I plot
POD.Rx.plotSummary(f_rom_2_r,'layout','built-in');
POD.Jnz.plotSummary(f_rom_2_j,'layout','built-in');
drawnow;
%% ROM-II test
[ROM_2] = solveDif(Nx,param.Nt,bc,D,'rom_stage',2,'encoder',POD);

%%
plotComp(FOM,ROM_2,D,f_sol_2);

%%
% %% Data-driven approximation demo
% fGrid = figure;
% fGrid.WindowStyle = 'docked';
% [parGrid,fGrid] = getParGrid(vals,bc.val,fGrid);
% drawnow;
% [nV,xV,in] = findSimplex(bc.val,vals);
% 
% %% ROM-III initial
% time_dd = tic;
% Yp = zeros(POD.X.N,length(nV));
% Ydd = zeros(POD.X.N,param.Nt + 1);
% for nt = 1:param.Nt+1
%     for i = 1:length(nV)
%         ns = nV(i);
%         Yp(:,i) = train_2{iNx}.sol{ns}.Y(:,nt);
%         Ydd(:,nt) = interpTriang(bc.val,vals(nV,:),Yp')';
%     end
% end
% Xdd = forward(Ydd,POD.X);
% time_dd = toc(time_dd);
% 
% ROMdd = ROM_2;
% ROMdd.X = Xdd;
% ROMdd.timelog.overall = time_dd;
% ROMdd.name = 'ROM_{DD}';
% 
% %% DD plot
% fRomDD = figure;
% fRomDD.WindowStyle = 'docked';
% [fRomDD] = plotComp(FOM,ROMdd,D,fRomDD);
% drawnow;
% 
% %% ROM-III test
% [ROM_3] = solveDif(Nx,param.Nt,bc,D,3,0,POD,Ydd);
% 
% %% ROM-II plot
% fRom3 = figure;
% fRom3.WindowStyle = 'docked';
% [fRom3] = plotComp(FOM,ROM_3,D,fRom3);
% drawnow;


