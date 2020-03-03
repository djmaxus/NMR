%% startup
clc;
close all;
addpath('./DIF1D');
addpath('./POD');
%% problem
% grid
Nx = 1000;
% boundary conditions
bc.first = true;
bc.val = [1,-1];
%% snapshot collection
% I: variable snapshots and POD
[data_1] = solveDif(Nx,bc,1);
POD.X = makePOD(data_1.X,true,false,0.99);
% II: residual and Jacobian snapshots and POD
[data_2] = solveDif(Nx,bc,2,POD);
% III: basis truncation
nR = 0.01;
nR = max(nR,POD.X.N/Nx);
POD.Rx = makePOD(data_2.Rx,false,false,nR);

nJ = 0.01;
nJ = max(nJ,POD.X.N/Nx/3);
POD.Jnz = makePOD(data_2.Jnz,true,true,nJ);

% IV: greedy sampling
nR = 0.1;
nR = max(nR,POD.Rx.N/Nx);
POD.Rx = makeOblique(nR,POD.Rx);
disp(norm(POD.Rx.O * data_2.Rx - data_2.Rx));

nJ = 0.1;
nJ = max(nJ,POD.Jnz.N/Nx / 3);
POD.Jnz = makeOblique(nJ,POD.Jnz);
disp(norm(POD.Jnz.O * data_2.Jnz - data_2.Jnz));

%% offline stage

% III: offline operators calculation
%% online stage
