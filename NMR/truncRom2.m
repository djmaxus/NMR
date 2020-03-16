function [POD] = truncRom2(POD,train,varargin)
    %TRAINROM1 Summary of this function goes here
    %   Detailed explanation goes here
    R = [];
    J = [];
    for i = 1:length(train.sol)
        R = [R,train.sol{i}.Rx];
        J = [J,train.sol{i}.Jnz];
    end
    fTrainR = figure;
    fTrainR.WindowStyle = 'docked';
    POD.Rx = makePOD(R,false,false,0.99,fTrainR);
    fTrainJ = figure;
    fTrainJ.WindowStyle = 'docked';
    POD.Jnz = makePOD(J,false,false,0.99,fTrainJ);
    
    POD.Rx.N = max(POD.Rx.N,POD.X.N);
    POD.Jnz.N = max(POD.Jnz.N,POD.X.N);
end

