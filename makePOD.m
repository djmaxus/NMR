function [POD] = makePOD(data, frac_energy)
    %MAKEPOD Summary of this function goes here
    %   Detailed explanation goes here

    %% Substract initial
%     POD.init = 0;
    POD.init = data(:,1);
    data = data - POD.init;

    %% Substract mean
    POD.mu = mean(data,2);
    data = data - POD.mu;
%     POD.mu = 0;
%     data = data./ sqrt(size(data,2)-1);

    %% Divide by variation
    POD.sigma = sqrt(var(data,0,2));
%     data = data.\(POD.sigma+eq(POD.sigma,0));

    %% SVD
    [POD.U,POD.S,POD.V] = svd(data);
    Nd = min(size(data));
    if and(isinteger(frac_energy), frac_energy > 1)
        POD.N = min(frac_energy,size(data,2));
    else
        POD.N = eigEnNum(diag(POD.S(1:Nd,1:Nd)),frac_energy);
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
    legend('eigen value', ['reduction, frac = ', num2str(frac_energy)]);
    xlabel('Number of POD components');
%     waitforbuttonpress;
end
