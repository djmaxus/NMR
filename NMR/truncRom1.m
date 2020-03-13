function [POD] = truncRom1(train,varargin)
    %TRAINROM1 Summary of this function goes here
    %   Detailed explanation goes here
    if length(varargin) >= 1
%         fig = varargin{1};
    else
        varargin = {};
    end
    fig.Name = 'ROM-1 Train';
    X = [];
    for i = 1:length(train.sol)
        X = [X,train.sol{i}.X];
    end
    POD = makePOD(X,true,false,0.99);
end

