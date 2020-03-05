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

subplot(2,3,1);
imagesc('CData',FOM.X,'XData',FOM.t,'YData',(1:Nx)./Nx);
caxis([Xmin Xmax]);
xlim([0 1]);
ylim([0 1]);
colorbar;
colormap jet;
title('FOM');
xlabel('t');
ylabel('x');
axis square;
if ~isempty(ROM)
    subplot(2,3,2);
    imagesc('CData',ROM.X,'XData',ROM.t,'YData',(1:Nx)./Nx);
    caxis([Xmin Xmax]);
    xlim([0 1]);
    ylim([0 1]);
    colorbar;
    title('ROM');
    xlabel('t');
    ylabel('x');
    axis square;
    
    subplot(2,3,3);
    imagesc('CData',ROM.X - FOM.X,'XData',FOM.t,'YData',(1:Nx)./Nx);
    xlim([0 1]);
    ylim([0 1]);
    colorbar;
    title('$\Delta ROM$','interpreter','latex');
    xlabel('t');
    ylabel('x');
    axis square;
end

%% Diffusion
subplot(2,3,4);
plot(D(1:Nx),(1:Nx)./Nx);
ylim([0,1]);
xlim([0,max(D)]);
ylabel('x');
xlabel('$D_0$','interpreter','latex');
title('Diffusion coefficient $D_0(x)$','interpreter','latex');

if ~isempty(ROM.X)
    
    subplot(2,3,5:6);
    enum = abs(ROM.X - FOM.X);
    denom = 3*sqrt(var(FOM.X,[],'all'));
    denom = denom + eq(denom,0);
    res = enum./denom;
    res = sum(res.^2,1)./size(enum,1);
    plot((1:Nt+1)./Nt,100*sqrt(res));
    xlim([0,1]);
    xlabel('t');
    ylabel('\Delta, %');
    title('Relatve $$MSE(t)$$','interpreter','latex');
end
legend('Relative MSE');

speed_up = FOM.timer / ROM.timer;
sms = ['Speed up: ',num2str(speed_up)];
disp(sms);
f.Name = sms;
