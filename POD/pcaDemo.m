function [fig] = pcaDemo(varargin)
    if length(varargin) >= 1
        fig = varargin{1};
    else
        fig = figure;
    end
    fig.Name = 'PCA-2D';
    xC = [2;1;];
    sig = [2;0.5];
    theta = pi/3;
    R = [ cos(theta),-sin(theta);
        sin(theta), cos(theta) ];
    
    nPoints = 2000;
    X = R * diag(sig) * randn(2,nPoints) + diag(xC) * ones(2,nPoints);
    xmin = min(X,[],'all');
    xmax = max(X,[],'all');
    
    scatter(X(1,:),X(2,:),'k.');
    hold on;
    box on;
    grid on;
    axis square;
    xlim([xmin xmax]);
    ylim([xmin xmax]);
   
    %%
    Xavg = mean(X,2);
    B = X - Xavg;
    
    [U,S,~] = svd(B/sqrt(nPoints-1),'econ');
    
    theta = (0:0.01:1.0)*2*pi;
    L = 3*U*S;
    [Xstd] = U*S*[cos(theta);sin(theta)];
    
    %% principal components
    waitforbuttonpress;
    plot([Xavg(1),Xavg(1) + L(1,1) ],[Xavg(2),Xavg(2) + L(2,1)],'-+r','LineWidth',1);
    waitforbuttonpress;
    plot([Xavg(1),Xavg(1) + L(1,2) ],[Xavg(2),Xavg(2) + L(2,2)],'-+r','LineWidth',1);
    %% ellipsoids
    for j = 1:3
        waitforbuttonpress;
        plot(Xavg(1) + j*Xstd(1,:),Xavg(2) + j*Xstd(2,:),'g-','LineWidth',1);
    end
end
