function [parGrid,fig] = getParGrid(X,val,fig)
    %PLOTPARGRID Summary of this function goes here
    %   Detailed explanation goes here
    clf(fig);
    ax = axes(fig);
    %%
    parGrid = delaunayn(X);
    for simp = parGrid'
        Xv = X(simp',:);
        Xv = [Xv;Xv(1,:)];
        plot(ax,Xv(:,1),Xv(:,2),'-ok');
        hold on;
    end
    if ~isempty(val)
        hold on;
        fTest = scatter(ax,val(:,1),val(:,2));
        hold off;
    end
    hold off;
    
    xmin = min(X,[],'all');
    xmax = max(X,[],'all');
    xlim(ax,[1.1*xmin 1.1*xmax]);
    ylim(ax,[1.1*xmin 1.1*xmax]);
    axis(ax,'square');
    
    xlabel(ax,'$b_0$','interpreter','latex');
    ylabel(ax,'$b_1$','interpreter','latex');
    
    legend(fTest,'Test sample');
end

