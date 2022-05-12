function [] = plotSol(where,sol, fontsize)
    if isa(where,'matlab.ui.Figure')
        clf(where);
        fig = where;
    elseif isa(where,'matlab.graphics.axis.Axes')
        cla(where);
        ax = where;
    else
        error('Whrere to plot, mate?');
    end
    
    if(~exist('fig','var'))
        error('No target figure');
    end
    
    if nargin == 2
        fontsize = 14;
    end
    
    Nx = size(sol.X,1);
    verts = linspace(0,1,Nx+2);
    
    label_settings = {'interpreter','latex'};
    
    sp1 = subplot(1,2,1,axes(fig));
    hold(sp1, 'on');
        yyaxis(sp1, 'left');
            bound = [0,1];
        
            plot(sp1, bound, sol.param.bc.val,'s');
        
            U0 = sol.X(:,1);
            plot(sp1, verts(2:end-1),U0,'-b');
            
            U05 = sol.X(:,round(length(sol.t)/2));
            plot(sp1, verts(2:end-1),U05,'-g');
            
            U1 = sol.X(:,end);
            plot(sp1, verts(2:end-1), U1,'-r');
                  
            xQ = xline(sp1,verts(sol.param.q.pos+1),'--',...
                ['q = ',num2str(sol.param.q.val(sol.param.q.pos))]);
            xQ.LabelOrientation = 'horizontal';
            
            ylabel(sp1, '$u$',label_settings{:});

        yyaxis(sp1, 'right');
            area(sp1, linspace(0,1,Nx+1),sol.param.D.val,'FaceColor','#999','LineStyle','none','FaceAlpha',0.5);
            ylabel(sp1, '$\kappa$',label_settings{:});
            ylim(sp1,[0,max(sol.param.D.val)]);

        xlim(bound);
        xlabel('$x$','interpreter','latex');
        
        leg = legend({'$u_\mathrm{B}$', ...
                '$u(x,0)$', ...
                '$u(x,\frac{1}{2})$',...
                '$u(x,1)$',...
                '$q$','$\kappa$'},label_settings{:});
            
        leg.Location = 'south';
        leg.NumColumns = 3;
        leg.Box = 'on';
            
        title(['Summary [',sol.name,']']);
        
    hold(sp1, 'off');
    sp1.FontSize = fontsize;
    axis(sp1,'square');
    
    sp2 = subplot(1,2,2,axes(fig));
        Nt = length(sol.t);
        mergeBC = @(ind) repmat(sol.param.bc.val(ind),[1,Nt]);
        cdata = @() [mergeBC(1);sol.X;mergeBC(2)];
        
        imagesc(sp2,'CData',cdata(),'XData',sol.t,'YData', linspace(0,1,Nx+2),'Interpolation','bilinear');
        
        colorbar(sp2);
        colormap(sp2,'jet');
        
        axis(sp2,'equal');
        
        xlim(sp2,[0,1]);
        ylim(sp2,[0,1]);
        
        xlabel(sp2,'$t$',label_settings{:});
        ylabel(sp2,'$x$',label_settings{:});

        title(sp2, '$u(x,t)$', label_settings{:});
        
    sp2.FontSize = fontsize;

end

function [] = plotSol2(where,sol)
    %PLOTSOL Summary of this function goes here
    %   Detailed explanation goes here

        if isa(where,'matlab.ui.Figure')
        clf(where);
        sp1 = subplot(2,2,1,axes(where));
        sp2 = subplot(2,2,2,axes(where));
        sp3 = subplot(2,2,3,axes(where));
        ax = subplot(2,2,4,axes(where));
        is_ext = true;
    elseif isa(where,'matlab.graphics.axis.Axes')
        is_ext = false;
    else
        error('Where to plot, mate?');
    end
    Nx = size(sol.X,1);
    if is_ext
        %% initial
        plot(sp1,(1:Nx)./Nx,sol.X(:,1));
        xlabel(sp1,'$x$','interpreter','latex');
        ylabel(sp1,'$u$','interpreter','latex');
        ttl = '$u^0(x)$';
        title(sp1,ttl,'interpreter','latex');
        %% diffusion

        plot(sp2,(0:Nx)./(Nx+1),sol.param.D.val);
        xlabel(sp2,'$x$','interpreter','latex');
        ylabel(sp2,'$D_0$','interpreter','latex');
        ttl = '$D_0(x)$';
        title(sp2,ttl,'interpreter','latex');
        %% relative slices
        lgnd = cell(1,3);
        tidx = [2;round(length(sol.t)/2);length(sol.t)];
        for i = 1:length(tidx)
            lgnd{i} = ['n = ',num2str(tidx(i)-1)];
        end
        plot(sp3,(1:Nx)./Nx,sol.X(:,tidx));
        xlabel(sp3,'$x$','interpreter','latex');
        ylabel(sp3,'$u$','interpreter','latex');
        ttl = '$u(x,t)$';
        title(sp3,ttl,'interpreter','latex');
        legend(sp3,lgnd{1},lgnd{2},lgnd{3});
        plotSolFom(ax,sol,Nx);
        figure(where);
    else
        plotSolFom(where,sol,Nx);
    end

end

function plotSolFom(ax,sol,Nx)
    imagesc(ax,'CData',sol.X,'XData',sol.t,'YData',(1:Nx)./Nx);
    xlim(ax,[0,1]);
    ylim(ax,[0,1]);
    colorbar(ax);
    colormap(ax,'jet');
%     ttl = ['$',sol.name,', t_{CPU} = ',num2str(sol.time.all),'$'];
    title(ax,sol.name,'interpreter','latex');
    xlabel(ax,'t');
    ylabel(ax,'x');
end
