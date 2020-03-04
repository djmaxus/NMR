function [POD] = makePOD(data, init, mu, param)
    %MAKEPOD Summary of this function goes here
    %   Detailed explanation goes here
    
    
    
    %% Substract initial
    if init
        POD.init = data(:,1);
        data(:,1) = [];
        data = data - POD.init;
    else
        POD.init = zeros(size(data(:,1)));
    end
    
    
    %% Substract mean
    if mu
        POD.mu = mean(data,2);
        data = data - POD.mu;
    else
        POD.mu = zeros(size(data(:,1)));
    end
    
    %% Crop Data
    %     cropData;
    %% SVD
    [POD.U,S,~] = svd(data,'econ');
    Nd = size(POD.U,2);
    if param > 1
        POD.N = min(round(param),Nd);
    elseif param > 0.5
        POD.N = eigEnNum(diag(S),param);
    elseif param > 0
        POD.N = round(param*Nd);
    else
        error('incorrect parameter');
    end
    %% visualize
    f = figure;
    f.WindowState = 'maximized';
    % dataset
    plot(POD.init,'b-','LineWidth',0.01);
    hold on;
    for i = 1:size(data,2)
        plot(POD.init + POD.mu + data(:,i),'k-','LineWidth',0.01);
        hold on;
    end
    % principal components
    U = POD.U*sqrt(S);
    for i = 1:POD.N
        plot(POD.init + POD.mu + U(:,i),'r-','LineWidth',1.0);
        hold on;
    end
    hold off;
    % reconstruction error
    X = POD.U * POD.U' * data;
    err = norm(data - X);
    title(['Reconstruction error: ',num2str(err)]);
    drawnow;
end
