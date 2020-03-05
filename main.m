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

%% grid
Nt = 100;
Ng = 3;
grids = round(logspace(log10(100),log10(500),Ng));
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
    plotNMR;
    
    %% speed-up plot
    sval(iG) = speed_up;
    if exist('fs','var')
        close(fs);
    end
    fs = figure;
    fs.WindowState = 'maximized';
    semilogx(grids,sval,'--o');
    xlabel('Grid size');
    ylabel('Speed-up');
    title('Speed-up plot');
    drawnow;
end
