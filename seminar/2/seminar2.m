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
NxMin = 2000;
NxMax = 2000;
NxCount = 1;
grids = round(logspace(log10(NxMin),log10(NxMax),NxCount));
%% problem setup
% space grid
iNx = NxCount;
Nx = grids(iNx);
% time
bc.kind = 1;
D =  D_init(0.9,1,0.5e-1);
x0 = setInit(Nx,0.5,2,1);
source_term = q_init(0.6,1e+3);
problem = setProblem(Nx,bc,source_term,D,x0,'FOM');
arg_jac = @(x) speye(length(x));

%% time & timestep settings
T = 1;
Nt = 200;
dt = T/Nt;

%% nonlinear solver settings
nlin_solver.kmax = 10;
nlin_solver.kmin = 3;
nlin_solver.crit_rel = 1e-3;
nlin_solver.crit_abs = 1e-3;
nlin_solver.eval = @(x,Res,Jac) newton_common(x,Res,Jac,nlin_solver.kmax,nlin_solver.crit_abs,nlin_solver.crit_rel,arg_jac) ;

%% solve FOM
clc;
dt_lim = [dt,dt];
[FOM] = solveDif(T,Nt,dt_lim,problem,nlin_solver);
plotSol(FOM,f.FOM); figure(f.FOM);

%% PCA 2D-demo
pcaDemo(f.pca2D);
figure(f.pca2D);

%% N-dimensional PCA demo
pcaDemoN(f.pcaND_1,f.pcaND_2);
figure(f.pcaND_1);
figure(f.pcaND_2);

%% train or load ROM-I data
[train.i] = genRom1(problem);

%% ROM-I train
[POD.X] = truncRom1({FOM});
POD.X.plotSummary(f.pod.x);
figure(f.pod.x);

%% ROM-I test
% encoder
arg_1     = @(y) POD.X.decode(y);
arg_jac_1 = @(y) POD.X.JD(y) ;
% problem
problem_rom_1 = setProblem(Nx,bc,source_term,D,POD.X.encode(x0),'ROM-I');
problem_rom_1.residual = @(x,dt) buildResidual(arg_1(x),dt,problem_rom_1.param);
problem_rom_1.jacobian = @(x,dt) buildJacobian(arg_1(x),dt,problem_rom_1.param);
% nonlinear solver
nlin_rom_1 = nlin_solver;
nlin_rom_1.eval = @(x,Res,Jac)(newton_common(x,Res,Jac,nlin_rom_1.kmax,nlin_rom_1.crit_abs,nlin_rom_1.crit_rel,arg_jac_1));
% solve
[ROM.i] = solveDif(T,Nt,[dt,dt],problem_rom_1,nlin_rom_1,'rom_stage',1,'encoder',POD);
ROM.i.X = arg_1(ROM.i.X);
plotComp(FOM,ROM.i,problem.param.D,f.rom.i);

%% train or load ROM-II data
[train.ii] = genRom2(problem,POD);

%% ROM-II train
[POD,Rx,Jnz] = truncRom2(POD,{ROM.i});
POD.Rx.plotSummary(f.pod.r); figure(f.pod.r);
POD.Jnz.plotSummary(f.pod.j); figure(f.pod.j);
%% Gappy-POD train
[POD] = approxRom2(POD,Rx,Jnz);

%% ROM-II test
% encoder
arg_2    = @(y) POD.X.decode(y);
arg_jac_2 = @(y) POD.X.JD(y);


% Nonlinear Function Approximation
nfa.res.mask = POD.Rx.idx_mask;
nfa.res.eval = @(W)(applyOblique(W,POD.Rx));

nfa.jac.mask = POD.Jnz.idx_mask;
nfa.jac.eval = @(Jf)(applyOblique(Jf.*[1,-1],POD.Jnz).*[1,-1]);

% problem
problem_rom_2 = setProblem(Nx,bc,source_term,D,POD.X.encode(x0),'ROM-II');
problem_rom_2.residual = @(x,dt) buildResidual(arg_2(x),dt,problem_rom_2.param,nfa.res);
problem_rom_2.jacobian = @(x,dt) buildJacobian(arg_2(x),dt,problem_rom_2.param,nfa.jac);

% nonlinear solvrer
nlin_rom_2 = nlin_solver;
nlin_rom_2.eval = @(x,Res,Jac)(newton_common(x,Res,Jac,nlin_rom_2.kmax,nlin_rom_2.crit_abs,nlin_rom_2.crit_rel,arg_jac_2));

% solve
[ROM.ii] = solveDif(T,Nt,[dt,dt],problem_rom_2,nlin_rom_2,'rom_stage',2,'encoder',POD);
ROM.ii.X = arg_2(ROM.ii.X);

plotComp(FOM,ROM.ii,problem.param.D,f.rom.ii);
figure(f.rom.ii);

%% Adaptive FOM
% FOM
fom_res = @(x,dt) buildResidual(x,dt,problem.param);
fom_jac = @(x,dt) buildJacobian(x,dt,problem.param);
fom_dec_jac = @(x) speye(length(x));
fom_nlin = nlin_solver;
fom_nlin.eval = @(x,Res,Jac) newton_common(x,Res,Jac,nlin_solver.kmax,nlin_solver.crit_abs,nlin_solver.crit_rel,fom_dec_jac);
% ROM-I    
rom_dec     = @(y) POD.X.decode(y);
rom_dec_jac = @(y) POD.X.JD(y) ;
rom_res = @(x,dt) buildResidual(rom_dec(x),dt,problem.param);
rom_jac = @(x,dt) buildJacobian(rom_dec(x),dt,problem.param);
rom_nlin = nlin_solver;
rom_nlin.eval = @(x,Res,Jac) newton_common(x,Res,Jac,nlin_solver.kmax,nlin_solver.crit_abs,nlin_solver.crit_rel,rom_dec_jac);
% Adaptive 
ada_res = @(is_rom) adaptive_fom(is_rom,fom_res,rom_res);
ada_jac = @(is_rom) adaptive_fom(is_rom,fom_jac,rom_jac);
ada_nlin = @(is_rom) adaptive_fom(is_rom,fom_jac,rom_jac);
%% Data-driven approximation demo
[parGrid] = getParGrid(vals,bc.val,f.grid);
figure(f.grid);
[nV,xV,in] = findSimplex(bc.val,vals);

%% ROM-III initial
time_dd = tic;
Yp = zeros(POD.X.npc(),length(nV));
Ydd = zeros(POD.X.npc(),problem.Nt + 1);
for nt = 1:problem.Nt+1
    for i = 1:length(nV)
        ns = nV(i);
        Yp(:,i) = train.ii{ns}.Y(:,nt);
        Ydd(:,nt) = interpTriang(bc.val,vals(nV,:),Yp')';
    end
end
Xdd = POD.X.decode(Ydd);
time_dd = toc(time_dd);


%% ROM-III test
[ROM.iii] = solveDif(Nx,Nt,problem,'rom_stage',2,'encoder',POD,'dd_approx',Ydd);
plotComp(FOM,ROM.iii,problem.D,f.rom.iii);
figure(f.rom.iii);

%% ROM-DD
ROM.dd = ROM.ii;
ROM.dd.X = Xdd;
ROM.dd.timelog.overall = time_dd;
ROM.dd.name = 'ROM_{DD}';

plotComp(FOM,ROM.dd,problem.D,f.rom.dd);
%% seminar III DMD
seminar3;
%%
seminar_encoder;

function [fun] = adaptive_rom(is_rom,fom,rom)
    if is_rom
        fun = rom;
    else
        fun = fom;
    end
end
