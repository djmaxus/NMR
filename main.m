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
Nx = 100;

%% problem
setProblem;

%% model reduction
reduceAll;

%% online stage
bc.val = val;
[FOM] = solveDif(Nx,bc,D);
[ROM] = solveDif(Nx,bc,D,3,POD);
plotNMR;
