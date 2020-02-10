addpath('./DIF1D');
addpath('./POD');
clear all;
close all;
clc;
Nx = 100;
tic;
[FOM] = solveDif(Nx,1:Nx,[]);
toc;
%% X POD
POD.X = makePOD(FOM.X,0.99);
%% Residual POD
POD.Rx = makePOD(FOM.Rx,0.99);
%% greedy mask
mask = greedy(POD.Rx.U);
frac = 0.1;
Nm = round(frac*Nx);
Nm = max(Nm,POD.Rx.N);
mask = mask(:,1:Nm);
%% oblique projector
POD.Rx.O = oblique(mask,POD.Rx.Ut);
idx_mask = [];
for i = 1:size(mask,2)
    for j=1:size(mask,1)
        if mask(j,i) ==1
            idx_mask = [idx_mask;j];
            break;
        end
    end
end
find_left = false;
find_right = false;
for i = idx_mask'
    if i == 1
        find_left = true;
    elseif i == Nx
        find_right = true;
    end
end
if ~find_left
    idx_mask = [idx_mask;1];
    dmask = zeros(Nx,1);
    dmask(1) = 1;
    mask = [mask,dmask];
end
if ~find_right
    idx_mask = [idx_mask;Nx];
    dmask = zeros(Nx,1);
    dmask(Nx) = 1;
    mask = [mask,dmask];
end
%% perform ROM
tic;
[ROM] = solveDif(Nx,idx_mask,POD);
toc;
%% plot
Xmin = min([FOM.X,ROM.X],[],'all');
Xmax = max([FOM.X,ROM.X],[],'all');
f = figure;
f.WindowState = 'maximized';
subplot(2,2,1);
imagesc(FOM.X); 
caxis([Xmin Xmax]);
colorbar;
colormap jet;

subplot(2,2,3);
imagesc(ROM.X);
caxis([Xmin Xmax]);
colorbar;
colormap jet;

subplot(2,2,2);
plot(FOM.X(:,end));
hold on;
plot(ROM.X(:,end));
subplot(2,2,4);
semilogy(abs(ROM.X(:,end) - FOM.X(:,end)));
