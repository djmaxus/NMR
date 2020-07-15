function [] = plotSol(sol,where)
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
        ax = where;
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
        
        plot(sp2,(0:Nx)./(Nx+1),sol.D);
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
        for nt = tidx
            plot(sp3,(1:Nx)./Nx,sol.X(:,nt));
            hold on;
        end
        xlabel(sp3,'$x$','interpreter','latex');
        ylabel(sp3,'$u$','interpreter','latex');
        ttl = '$u(x,t)$';
        title(sp3,ttl,'interpreter','latex');
        hold off;
        legend(sp3,lgnd{1},lgnd{2},lgnd{3});
    end
        imagesc(ax,'CData',sol.X,'XData',sol.t,'YData',(1:Nx)./Nx);
        xlim(ax,[0,1]);
        ylim(ax,[0,1]);
        colorbar(ax);
        colormap(ax,'jet');
        ttl = ['$',sol.name,' , t_{CPU} = ',num2str(sol.timelog.overall),'$'];
        title(ax,ttl,'interpreter','latex');
        xlabel(ax,'t');
        ylabel(ax,'x');
    drawnow;
end

