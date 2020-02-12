
% TODO: calculate X as supporting variable for ROM
% TODO: extend mask properly
addpath('./DIF1D');
addpath('./POD');
close all;
clc;
ROM = [];

Nx = 500;
bc.first = false;
bc.val = [2,2];
[dataset] = solveDif(Nx,[],bc,true);
%% do model reduction
tic;
[POD] = doNMR(dataset,0.01,0.1,0.2);
toc;
%% perform ROM and compare
bc.val = [1,1];
[FOM] = solveDif(Nx,[],bc,false);
[ROM] = solveDif(Nx,POD,bc,false);
%% plot
if isempty(ROM)
    ROM.X = [];
    ROM.t = [];
end
Xmin = min([FOM.X,ROM.X],[],'all');
Xmax = max([FOM.X,ROM.X],[],'all');
Tmin = min([FOM.t,ROM.t],[],'all');
Tmax = max([FOM.t,ROM.t],[],'all');

f = figure;
f.WindowState = 'maximized';

subplot(2,2,1);
imagesc('CData',FOM.X,'XData',FOM.t,'YData',(1:Nx)./Nx);
caxis([Xmin Xmax]);
xlim([Tmin Tmax]);
ylim([1/Nx 1]);
colorbar;
colormap jet;

if ~isempty(ROM)
    subplot(2,2,3);
    imagesc('CData',ROM.X,'XData',ROM.t,'YData',(1:Nx)./Nx);
    caxis([Xmin Xmax]);
    xlim([Tmin Tmax]);
    ylim([1/Nx 1]);
    colorbar;
    colormap jet;
end

subplot(2,2,2);
plot(FOM.X(:,end));

if ~isempty(ROM.X)
    hold on;
    plot(ROM.X(:,end));
    subplot(2,2,4);
    enum = abs(ROM.X(:,end) - FOM.X(:,end));
    denom = Xmax - Xmin;
    denom = denom + eq(denom,0);
    res = enum./denom;
    plot(res);
end

speed_up = FOM.timer / ROM.timer;
disp(['Speed up: ',num2str(speed_up)]);
