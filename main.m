%% startup
clc;
close all;
addpath('./DIF1D');
addpath('./POD');
addpath('./NMR');

%% train parameter samples
% boundary condition scale
pMax = 1;
setParSample;

% boundary condition test value
val = rand(1,2) * 2 * pMax - pMax;
DG = delaunayn(vals);
F = sin(0.5 * vals(:,1) + vals(:,2)* 0.5);
Ftest= interpTriang(val,vals,DG,F);
%% grid
Nt = 100;
Ng = 1;
grids = round(logspace(log10(1000),log10(1000),Ng));
sval = zeros(Ng,1);
fs = figure;
for iG = 1:Ng  
    Nx =  grids(iG);
    %% problem
    setProblem;
    
    %% model reduction
    reduceAll;
    
    %% online stage
    bc.val = val;
    [FOM] = solveDif(Nx,Nt,bc,D);
    [ROM] = solveDif(Nx,Nt,bc,D,3,POD);
    disp(ROM.k_avg);
    [ROM] = solveDif(Nx,Nt,bc,D,3,POD,vals,DG,val);
    plotNMR;
    
    %% speed-up plot
    sval(iG) = speed_up;
    if exist('fs','var')
        close(fs);
    end
    fs = figure;
    fs.WindowState = 'maximized';
    semilogx(grids(1:iG),sval(1:iG),'--o');
    xlabel('Grid size');
    ylabel('Speed-up');
    title('Speed-up plot');
    drawnow;
end
