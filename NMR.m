
% TODO: calculate X as supporting variable for ROM
% TODO: extend mask properly
addpath('./DIF1D');
addpath('./POD');
clc;
close all;
ROM = [];

Nx = 1000;

bc.first = true;
bc.val = [1,-1];
[dataset] = solveDif(Nx,[],bc,true);



% bc.first = true;
% bc.val = [1,-1];
% [dataset1] = solveDif(Nx,[],bc,true);
% 
% figure;
% imagesc('CData',dataset.Rx);
% figure;
% imagesc('CData',dataset1.Rx(:,end-100:end),'XData',dataset1.t(:,end-100:end));
% 
% dataset.X = [dataset.X,dataset1.X];
% dataset.Rx = [dataset.Rx,dataset1.Rx];
% dataset.Jnz = [dataset.Jnz,dataset1.Jnz];
%% do model reduction
[POD] = doNMR(dataset,0.99,0.99,0.2);
%% perform ROM and compare
% bc.first = false;
% bc.val = [1,1];
[FOM] = solveDif(Nx,[],bc,false);
[ROM] = solveDif(Nx,POD,bc,false);
plotNMR;
