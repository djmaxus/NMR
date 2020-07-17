function [] = plotComp(sol_1,sol_2,D,fig)
    %PLOTCOMP Summary of this function goes here
    %   Detailed explanation goes here
    %% rename
    FOM = sol_1;
    ROM = sol_2;
    %% limits
    Xmin = min([FOM.X,ROM.X],[],'all');
    Xmax = max([FOM.X,ROM.X],[],'all');
    %%
    clf(fig);
    %% FOM
    f1 = subplot(2,3,1,axes(fig));
    plotSol(FOM,f1);
    %% ROM
    f2 = subplot(2,3,2,axes(fig));
    plotSol(ROM,f2);
    
    %% error
    x_rom = interp1(ROM.t,ROM.X',FOM.t)';
    enum = (x_rom - FOM.X);
    res = enum;
    
    sp3 = subplot(2,3,3,axes(fig));
    Nx = size(FOM.X,1);
    imagesc(sp3,'CData',res,'XData',FOM.t,'YData',(1:Nx)./Nx);
    xlim(sp3,[0 1]);
    ylim(sp3,[0 1]);
    colorbar(sp3);
    colormap(sp3,'gray');
    title(sp3,'$\Delta ROM$','interpreter','latex');
    xlabel(sp3,'t');
    ylabel(sp3,'x');
    
    
    %% Diffusion
    sp4 = subplot(2,3,4,axes(fig));
    plot(sp4,D(1:Nx),(1:Nx)./Nx);
    ylim(sp4,[0,1]);
    xlim(sp4,[0,max(D)]);
    ylabel(sp4,'x');
    xlabel(sp4,'$D_0$','interpreter','latex');
    title(sp4,'Diffusion coefficient $D_0(x)$','interpreter','latex');
    
    
    %% error dynamics
    Nt = length(FOM.t)-1;
    sp5 = subplot(2,3,5:6,axes(fig));
    res = sqrt(mean(enum.^2,1));
    plot(sp5,(0:Nt)./(Nt),res);
    xlim(sp5,[0,1]);
    xlabel(sp5,'t');
    ylabel(sp5,'\Delta');
    title(sp5,'$$MSE(t)$$','interpreter','latex');
    
    spd.all = FOM.timelog.overall / ROM.timelog.overall;
    spd.lin =  FOM.timelog.linsol / ROM.timelog.linsol;
    spd.bld =  FOM.timelog.bld / ROM.timelog.bld;
    spd.iter =  ROM.timelog.iters.bld_avg;
    sms = ['ALL = ',num2str(spd.all),newline(),...
           'LIN = ',num2str(spd.lin),newline(),...
           'BLD = ',num2str(spd.bld),newline(),...
           'ITR = ',num2str(spd.iter)];
    disp(['----',fig.Name,newline(),sms,newline(),'----']);
%     fig.Name = sms;
    drawnow;
    figure(fig);
end

