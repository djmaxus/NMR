
% TODO: calculate X as supporting variable for ROM
% TODO: extend mask properly
addpath('./DIF1D');
addpath('./POD');
clc;
close all;
ROM = [];

Nx = 100;

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
tic;
[POD] = doNMR(dataset,0.99,0.2,0.2);
toc;
%% perform ROM and compare
bc.first = true;
bc.val = [1,-1];
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
title('FOM');
xlabel('t');
ylabel('x');
if ~isempty(ROM)
    subplot(2,2,3);
    imagesc('CData',ROM.X,'XData',ROM.t,'YData',(1:Nx)./Nx);
    caxis([Xmin Xmax]);
    xlim([Tmin Tmax]);
    ylim([1/Nx 1]);
    colorbar;
    colormap jet;
    title('ROM');
    xlabel('t');
    ylabel('x');
end

subplot(2,2,2);
plot(FOM.X(:,end));
xlabel('x');
ylabel('U');
title('t = T');

if ~isempty(ROM.X)
    hold on;
    plot(ROM.X(:,end));
    legend('FOM','ROM');
    
    subplot(2,2,4);
    enum = abs(ROM.X(:,end) - FOM.X(:,end));
    denom = Xmax - Xmin;
    denom = denom + eq(denom,0);
    res = enum./denom;
    plot(res);
    xlabel('x');
    ylabel('\Delta');
    title('t = T');
end
legend('relative error');

speed_up = FOM.timer / ROM.timer;
sms = ['Speed up: ',num2str(speed_up)];
disp(sms);
f.Name = sms;

