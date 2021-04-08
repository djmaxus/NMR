function [] = plotParSample(fig,vals,val)
    %PLOTPARSAMPLE Summary of this function goes here
    %   Detailed explanation goes here
    clf(fig);
    ax = axes(fig);
    %%
    xmin = min(vals,[],'all');
    xmax = max(vals,[],'all');
    
    scatter(ax,vals(:,1),vals(:,2));
    
    if ~isempty(val)
        hold on;
        scatter(ax,val(:,1),val(:,2));
        hold off;
    end
    
    xlim(ax,[1.1*xmin 1.1*xmax]);
    ylim(ax,[1.1*xmin 1.1*xmax]);
    axis(ax,'square');
    
    xlabel(ax,'$b_0$','interpreter','latex');
    ylabel(ax,'$b_1$','interpreter','latex');
    
    legend(ax,'Train samples','Test sample');
    title(ax,'Parameter samples');
    
    box(ax,'on');
    grid(ax,'on');
end

