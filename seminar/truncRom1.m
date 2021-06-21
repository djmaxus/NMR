function [POD] = truncRom1(train)
    %TRAINROM1 Summary of this function goes here
    %   Detailed explanation goes here
    X = [];
    for i = 1:length(train)
        X = [X,train{i}.X];
    end
    x0 = X(:,1);
    POD = PCA(X,0.99,'init',x0);
end

