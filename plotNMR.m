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

%% Diffusion
subplot(2,2,2);
plot(D);
xlim([0,length(D)]);
ylim([0,max(D)]);
xlabel('x');
ylabel('$D_0$','interpreter','latex');
title('Diffusion coefficient $D_0$','interpreter','latex');

% subplot(2,2,2);
% plot(FOM.X(:,end));
% xlabel('x');
% ylabel('U');
% title('t = T');

if ~isempty(ROM.X)
%     hold on;
%     plot(ROM.X(:,end));
%     legend('FOM','ROM');
    
    subplot(2,2,4);
    enum = abs(ROM.X - FOM.X);
    denom = 1;...Xmax - Xmin;
    denom = denom + eq(denom,0);
    res = enum./denom./size(enum,1);
    res = sum(res.^2,1);
    plot(sqrt(res));
    xlabel('t');
    ylabel('\Delta');
    title('$$MSE(t)$$','interpreter','latex');
end
legend('MSE');

speed_up = FOM.timer / ROM.timer;
sms = ['Speed up: ',num2str(speed_up)];
disp(sms);
f.Name = sms;
