function [fig] = plotParSample(vals,val,varargin)
    %PLOTPARSAMPLE Summary of this function goes here
    %   Detailed explanation goes here
    if length(varargin) >= 1
        fig = varargin{1};
    else
        fig = figure;
    end
    fig.Name = 'Parameters';
    xmin = min(vals,[],'all');
    xmax = max(vals,[],'all');
    
    scatter(vals(:,1),vals(:,2));
    
    if ~isempty(val)
        hold on;
        scatter(val(:,1),val(:,2));
        hold off;
    end
    
    xlim([1.1*xmin 1.1*xmax]);
    ylim([1.1*xmin 1.1*xmax]);
    axis square;
    
    xlabel('$b_0$','interpreter','latex');
    ylabel('$b_1$','interpreter','latex');
    
    legend('Train samples','Test sample');
    title('Points in parameter samples');
end

