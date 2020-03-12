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

%% grid
Nt = 100;
Ng = 1;
grids = round(logspace(log10(500),log10(500),Ng));
sval = zeros(Ng,1);
fs = figure;
fs.WindowStyle = 'docked';
for iG = 1:Ng  
    Nx =  grids(iG);
    %% problem
    setProblem;
    
    %% model reduction
    reduceAll;
    
    %% online stage
    % boundary condition test value
    val = rand(1,2) * 2 * pMax - pMax;
    [nV,xV,in] = findSimplex(val,vals);
    bc.val = val;
    
    %% FOM
    [FOM] = solveDif(Nx,Nt,bc,D);
    
    %% ROM with default initial approximation
    try
        [ROM] = solveDif(Nx,Nt,bc,D,3,POD);
        disp(['iters per step: ',num2str(ROM.k_avg)]);
        plotNMR;
    end
    
    %% ROM with data-driven initial approximation
    try
        [ROM] = solveDif(Nx,Nt,bc,D,3,POD,xV,nV,val);
        disp(['iters per step: ',num2str(ROM.k_avg)]);
        plotNMR;
    end
    
    %% speed-up plot
    sval(iG) = speed_up;
    if exist('fs','var')
        close(fs);
    end
    
    %%
    fs = figure;
    fs.WindowStyle = 'docked';
    semilogx(grids(1:iG),sval(1:iG),'--o');
    xlabel('Grid size');
    ylabel('Speed-up');
    title('Speed-up plot');
    drawnow;
    
end
% save('./log/last');
