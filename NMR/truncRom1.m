function [POD] = truncRom1(train,varargin)
    %TRAINROM1 Summary of this function goes here
    %   Detailed explanation goes here
    if length(varargin) >=1
        to_plot = varargin{1};
    else
        to_plot = false;
    end
    X = [];
    for i = 1:length(train.sol)
        X = [X,train.sol{i}.X];
    end
    if to_plot
        fTrain1 = figure;
        fTrain1.WindowStyle = 'docked';
        POD = makePOD(X,true,false,0.99,fTrain1);
    else
        POD = makePOD(X,true,false,0.99);
    end
end

