function [POD] = approxRom2(POD,train,varargin)
    %TRAINROM1 Summary of this function goes here
    %   Detailed explanation goes here
    R = [];
    J = [];
    for i = 1:length(train.sol)
        R = [R,train.sol{i}.Rx];
        J = [J,train.sol{i}.Jnz];
    end
    
    nR = 0.75;
    nR = round(POD.Rx.N * (1 + nR));
    nR = nR/size(POD.Rx.U,1);
    POD.Rx = makeOblique(nR,POD.Rx);
    disp(norm(applyOblique(R,POD.Rx) - R));

    nJ = 0.75;
    nJ = round(POD.Jnz.N * (1 + nJ));
    nJ = nJ/size(POD.Jnz.U,1);
    POD.Jnz = makeOblique(nJ,POD.Jnz);
    disp(norm(applyOblique(J,POD.Jnz) - J));
end
