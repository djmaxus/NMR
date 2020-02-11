function [POD] = makePOD(data, init,mu, param)
    %MAKEPOD Summary of this function goes here
    %   Detailed explanation goes here

    %% Substract initial
%     POD.init = 0;
%     POD.init = data(:,1);

    if init
        POD.init = data(:,1);
    else
        POD.init = zeros(size(data(:,1)));
    end
    data = data - POD.init;

    %% Substract mean
    if mu
        POD.mu = mean(data,2);
    else
        POD.mu = zeros(size(data(:,1)));
    end
    data = data - POD.mu;

    %% Divide by variation
%     POD.sigma = 1;
%     POD.sigma = sqrt(var(data,0,2));
%     POD.sigma = POD.sigma + eq(POD.sigma,0);
%     data = data./POD.sigma;
    %% SVD
    [POD.U,POD.S,POD.V] = svd(data);
    Nd = min(size(data));
    if and(isinteger(param), param > 1)
        POD.N = min(param,size(data,2));
    else
        POD.N = eigEnNum(diag(POD.S(1:Nd,1:Nd)),param);
    end
    POD.Ut = POD.U(:,1:POD.N);

    %% plot
    f = figure;
    f.WindowState = 'maximized';
    en_frac = diag(POD.S)./sum(diag(POD.S));

    semilogy(en_frac);
    xlim([1 Nd]);

    hold on;
    xline(double(POD.N),'-k',num2str(POD.N));

    title('POD');
    legend('eigen value', ['reduction, frac = ', num2str(param)]);
    xlabel('Number of POD components');
%     waitforbuttonpress;
end
