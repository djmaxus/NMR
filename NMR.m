
% TODO: calculate X as supporting variable for ROM
% TODO: try various projections
% TODO: extend mask properly
% TODO: build Jacobian gappy
addpath('./DIF1D');
addpath('./POD');
close all;
clc;
ROM = [];

Nx = 500;
[FOM] = solveDif(Nx,[]);
%% do model reduction
tic;
[POD] = doNMR(FOM,0.1,0.1,0.1);
toc;
%% perform ROM
[ROM] = solveDif(Nx,POD);
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
disp(['Speed up: ',num2str(FOM.timer / ROM.timer)]);
