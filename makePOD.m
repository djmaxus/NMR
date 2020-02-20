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
%         POD.mu = mean(data,2);
        POD.mu = median(data,2);
        data = data - POD.mu;
    else
        POD.mu = zeros(size(data(:,1)));
    end
    
    
    %% calculate variance
    POD.sigma = sqrt(var(data,1,2));
    POD.sigma = POD.sigma + eq(POD.sigma,0);
    %     data = data ./ POD.sigma;
    
    %% Crop Data
    cropData;
    %% SVD
    [POD.U,POD.S,POD.V] = svd(data,'econ');
    Nd = size(POD.U,2);
    if and(isinteger(param), param > 1)
        POD.N = min(param,Nd);
    else
        if param <= 0
            error('incorrect parameter');
        elseif param < 0.5
            POD.N = round(param*Nd);
        elseif param <= 1
            POD.N = eigEnNum(diag(POD.S),param);
        else
            error('incorrect parameter');
        end
    end
    POD.Ut = POD.U(:,1:POD.N);
end
