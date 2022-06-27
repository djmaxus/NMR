function [] = plotComp(fig,ref,sol,skip_performance)
%PLOTCOMP Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    skip_performance = false;
end
num_rows = 1 + ~skip_performance;

clf(fig);

tiledlayout(num_rows,3,'TileSpacing','tight');
%% Reference
f1 = nexttile;
plotSol(f1,ref);
title(f1,ref.name,'interpreter','none');

%% Solution
f2 = nexttile;
plotSol(f2,sol);
title(f2,sol.name,'interpreter','none');

%% Difference
x_rom = interp1(sol.t,sol.X',ref.t)';
enum = (x_rom - ref.X);
res = abs(enum);

sp3 = nexttile;
Nx = size(ref.X,1);
plotField(sp3, res, ref.t, linspace(0,1,Nx),['Absolute', newline(), 'difference'],{});
set(sp3,'ColorScale','log');
colormap(sp3,'parula');
%% Diffusion
%     sp4 = subplot(2,3,4,axes(fig));
%     plot(sp4,D.val(1:Nx),(1:Nx)./Nx);
%     ylim(sp4,[0,1]);
%     xlim(sp4,[0,max(D.val)]);
%     ylabel(sp4,'x');
%     xlabel(sp4,'$D_0$','interpreter','latex');
%     title(sp4,'Diffusion coefficient $D_0(x)$','interpreter','latex');

%% Dynamic difference
%     Nt = length(ref.t)-1;
%     sp5 = subplot(2,3,4:6,axes(fig));
%     res = sqrt(mean(enum.^2,1));
%     semilogy(sp5,(0:Nt)./(Nt),res);
%     xlim(sp5,[0,1]);
%     xlabel(sp5,'t');
%     ylabel(sp5,'\Delta');
%     title(sp5,'$$MSE(t)$$','interpreter','latex');

%     spd.all = 0;%FOM.time.all / ROM.time.all;
%     spd.lin =  0;%FOM.time.solve / ROM.time.solve;
%     spd.bld =  0;%FOM.time.build / ROM.time.build;
%     spd.iter =  0;%ROM.iters_avg;
%     sms = ['ALL = ',num2str(spd.all),newline(),...
%            'LIN = ',num2str(spd.lin),newline(),...
%            'BLD = ',num2str(spd.bld),newline(),...
%            'ITR = ',num2str(spd.iter)];
%     disp(['----',fig.Name,newline(),sms,newline(),'----']);


%% comparative performace
if ~skip_performance
sp2 = nexttile;
% sp2 = axes(t);
% sp2.Layout.Tile = 4;
sp2.Layout.TileSpan = [1 3];

categories = {'linear solve','LSPG','nonlinearity','total'};
cat_array = categorical(categories);
cat_array = reordercats(cat_array,categories);

spd = calcSpeedUp(sol,ref);

perf_data = [spd.ref.solve, spd.sol.solve;...
    spd.ref.proj,  spd.sol.proj;...
    spd.ref.build, spd.sol.build;...
    spd.ref.all,   spd.sol.all  ];

b = bar(sp2, cat_array, perf_data);

set(sp2, 'YScale', 'log');
ylabel(sp2, 'Time [seconds]');
title(sp2, 'Relative performace');

for i = 1:2
xtips_i = b(i).XEndPoints;

ytips_i = b(i).YEndPoints;

labels_i = split(num2str(b(i).YData,'% .2g'));

text(xtips_i,ytips_i,labels_i,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
end

legend(sp2, {ref.name, sol.name},'Location',"best");
end

%% Focus
figure(fig);
end

function [] = plotCompOld(fig,ref,sol,fontsize)
%PLOTCOMP Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3
    fontsize = 12;
end

clf(fig);

%% Reference
f1 = subplot(2,3,1,axes(fig));
plotSol(f1,ref, fontsize);

%% Solution
f2 = subplot(2,3,2,axes(fig));
plotSol(f2,sol, fontsize);

%% Difference
x_rom = interp1(sol.t,sol.X',ref.t)';
enum = (x_rom - ref.X);
res = abs(enum);

sp3 = subplot(2,3,3,axes(fig));
Nx = size(ref.X,1);
plotField(sp3, res, ref.t, linspace(0,1,Nx),['Absolute', newline(), 'difference'],{});
set(sp3,'ColorScale','log');
colormap(sp3,'parula');
sp3.FontSize = fontsize;
%% Diffusion
%     sp4 = subplot(2,3,4,axes(fig));
%     plot(sp4,D.val(1:Nx),(1:Nx)./Nx);
%     ylim(sp4,[0,1]);
%     xlim(sp4,[0,max(D.val)]);
%     ylabel(sp4,'x');
%     xlabel(sp4,'$D_0$','interpreter','latex');
%     title(sp4,'Diffusion coefficient $D_0(x)$','interpreter','latex');

%% Dynamic difference
%     Nt = length(ref.t)-1;
%     sp5 = subplot(2,3,4:6,axes(fig));
%     res = sqrt(mean(enum.^2,1));
%     semilogy(sp5,(0:Nt)./(Nt),res);
%     xlim(sp5,[0,1]);
%     xlabel(sp5,'t');
%     ylabel(sp5,'\Delta');
%     title(sp5,'$$MSE(t)$$','interpreter','latex');

%     spd.all = 0;%FOM.time.all / ROM.time.all;
%     spd.lin =  0;%FOM.time.solve / ROM.time.solve;
%     spd.bld =  0;%FOM.time.build / ROM.time.build;
%     spd.iter =  0;%ROM.iters_avg;
%     sms = ['ALL = ',num2str(spd.all),newline(),...
%            'LIN = ',num2str(spd.lin),newline(),...
%            'BLD = ',num2str(spd.bld),newline(),...
%            'ITR = ',num2str(spd.iter)];
%     disp(['----',fig.Name,newline(),sms,newline(),'----']);


%% comparative performace
sp2 = subplot(2,3,4:6,axes(fig));

categories = {'linear problem','LSPG projection','nonlinear functions','overall'};
cat_array = categorical(categories);
cat_array = reordercats(cat_array,categories);

spd = calcSpeedUp(sol,ref);

perf_data = [spd.ref.solve, spd.sol.solve;...
    spd.ref.proj,  spd.sol.proj;...
    spd.ref.build, spd.sol.build;...
    spd.ref.all,   spd.sol.all  ];

b = bar(sp2, cat_array, perf_data);

set(sp2, 'YScale', 'log');
ylabel(sp2, 'Time [seconds]');
title(sp2, 'Relative performace');

for i = 1:2
xtips_i = b(i).XEndPoints;

ytips_i = b(i).YEndPoints;

labels_i = split(num2str(b(i).YData,'% .2g'));

text(xtips_i,ytips_i,labels_i,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
end

sp2.FontSize = fontsize;

legend(sp2, {ref.name, sol.name},'Location',"best");

%% Focus
figure(fig);
end
