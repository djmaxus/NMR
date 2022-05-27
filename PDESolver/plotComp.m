function [] = plotComp(fig,ref,sol,fontsize)
    %PLOTCOMP Summary of this function goes here
    %   Detailed explanation goes here

    if nargin == 3
        fontsize = 14;
    end
    
    clf(fig);

    %% Reference
    f1 = subplot(1,3,1,axes(fig));
    plotSol(f1,ref, fontsize);

    %% Solution
    f2 = subplot(1,3,2,axes(fig));
    plotSol(f2,sol, fontsize);

    %% Difference
    x_rom = interp1(sol.t,sol.X',ref.t)';
    enum = (x_rom - ref.X);
    res = abs(enum);

    sp3 = subplot(1,3,3,axes(fig));
    Nx = size(ref.X,1);
    plotField(sp3, res, ref.t, linspace(0,1,Nx),'Abs. difference',{});
 
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

%% Focus
    figure(fig);

end
