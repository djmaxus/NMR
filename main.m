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
Ng = 10;
grids = round(logspace(log10(100),log10(10000),Ng));
sval = [];
for Nx = grids
    %% problem
    setProblem;
    
    %% model reduction
    reduceAll;
    
    %% online stage
    bc.val = val;
    [FOM] = solveDif(Nx,bc,D);
    [ROM] = solveDif(Nx,bc,D,3,POD);
    plotNMR;
    sval = [sval,speed_up];
    if exist('fs','var')
        close(fs);
    end
    fs = figure;
    fs.WindowState = 'maximized';
    semilogx(grids(1:length(sval)),sval,'--o');
    xlabel('Grid size');
    ylabel('Speed-up');
    title('Speed-up plot');
    drawnow;
end


