function [fig] = plotComp(sol_1,sol_2,D,fig)
    %PLOTCOMP Summary of this function goes here
    %   Detailed explanation goes here
    %% rename
    FOM = sol_1;
    ROM = sol_2;
    %% limits
    Xmin = min([FOM.X,ROM.X],[],'all');
    Xmax = max([FOM.X,ROM.X],[],'all');
    
    %% FOM
    f1 = subplot(2,3,1);
    plotSol(FOM,f1);
    axis square;
    %% ROM
    f2 = subplot(2,3,2);
    plotSol(ROM,f2);
    axis square;
    
    %% error
    x_rom = interp1(ROM.t,ROM.X',FOM.t)';
    enum = abs(x_rom - FOM.X);
    res = enum;
    
    subplot(2,3,3);
    Nx = size(FOM.X,1);
    imagesc('CData',res,'XData',FOM.t,'YData',(1:Nx)./Nx);
    xlim([0 1]);
    ylim([0 1]);
    colorbar;
    title('$|\Delta ROM|$','interpreter','latex');
    xlabel('t');
    ylabel('x');
    axis square;
    
    
    %% Diffusion
    subplot(2,3,4);
    plot(D(1:Nx),(1:Nx)./Nx);
    ylim([0,1]);
    xlim([0,max(D)]);
    ylabel('x');
    xlabel('$D_0$','interpreter','latex');
    title('Diffusion coefficient $D_0(x)$','interpreter','latex');
    
    
    %% error dynamics
    Nt = length(FOM.t)-1;
    subplot(2,3,5:6);
    res = sqrt(mean(enum.^2,1));
    plot((0:Nt)./(Nt),sqrt(res));
    xlim([0,1]);
    xlabel('t');
    ylabel('\Delta, %');
    title('$$MSE(t)$$','interpreter','latex');
    
    legend('Relative MSE');
    
    spd.all = FOM.timelog.overall / ROM.timelog.overall;
    spd.lin =  FOM.timelog.linsol / ROM.timelog.linsol;
    spd.bld =  FOM.timelog.bld / ROM.timelog.bld;
    spd.iter =  ROM.timelog.iters.bld_avg;
    sms = ['ALL = ',num2str(spd.all),newline(),...
           'LIN = ',num2str(spd.lin),newline(),...
           'BLD = ',num2str(spd.bld),newline(),...
           'ITR = ',num2str(spd.iter)];
    disp(['----',newline(),sms,newline(),'----']);
    fig.Name = sms;
    
end

