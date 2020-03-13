function [fig] = plotSol(sol,varargin)
    %PLOTSOL Summary of this function goes here
    %   Detailed explanation goes here
    if length(varargin) >= 1
        fig = varargin{1};
    else
        fig = figure;
    end
    if length(varargin) >= 2
        is_ext = varargin{2};
    else
        is_ext = 0;
    end
    Nx = size(sol.X,1);
    if is_ext
        % initial
        subplot(2,2,1);
        plot((1:Nx)./Nx,sol.X(:,1));
        xlabel('$x$','interpreter','latex');
        ylabel('$u$','interpreter','latex');
        ttl = '$u^0(x)$';
        title(ttl,'interpreter','latex');
        % diffusion
        subplot(2,2,2);
        plot((0:Nx)./(Nx+1),sol.D);
        xlabel('$x$','interpreter','latex');
        ylabel('$D_0$','interpreter','latex');
        ttl = '$D_0(x)$';
        title(ttl,'interpreter','latex');
        % relative slices
        subplot(2,2,3);
        lgnd = cell(1,3);
        tidx = [2;round(length(sol.t)/2);length(sol.t)];
        for i = 1:length(tidx)
            lgnd{i} = ['n = ',num2str(tidx(i)-1)];
        end
        for nt = tidx
            plot((1:Nx)./Nx,sol.X(:,nt));
            hold on;
        end
        xlabel('$x$','interpreter','latex');
        ylabel('$u$','interpreter','latex');
        ttl = '$u(x,t)$';
        title(ttl,'interpreter','latex');
        hold off;
        legend(lgnd{1},lgnd{2},lgnd{3});
        % whole solution
        subplot(2,2,4);
    end
        imagesc('CData',sol.X,'XData',sol.t,'YData',(1:Nx)./Nx);
        xlim([0,1]);
        ylim([0,1]);
        colorbar;
        colormap jet;
        ttl = ['$',sol.name,' , t_{CPU} = ',num2str(sol.timelog.overall),'$'];
        title(ttl,'interpreter','latex');
        xlabel('t');
        ylabel('x');
end

