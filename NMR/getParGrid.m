function [parGrid,fig] = getParGrid(X,val,varargin)
    %PLOTPARGRID Summary of this function goes here
    %   Detailed explanation goes here
    if length(varargin) >= 1
        fig = varargin{1};
    else
        fig = figure;
    end
    fig.Name = 'Area triangulation';
    parGrid = delaunayn(X);
    
    for simp = parGrid'
        Xv = X(simp',:);
        Xv = [Xv;Xv(1,:)];
        plot(Xv(:,1),Xv(:,2),'-ok');
        hold on;
    end
    if ~isempty(val)
        hold on;
        fTest = scatter(val(:,1),val(:,2));
        hold off;
    end
    hold off;
    
    xmin = min(X,[],'all');
    xmax = max(X,[],'all');
    xlim([1.1*xmin 1.1*xmax]);
    ylim([1.1*xmin 1.1*xmax]);
    axis square;
    
    xlabel('$b_0$','interpreter','latex');
    ylabel('$b_1$','interpreter','latex');
    
    legend(fTest,'Test sample');
end

