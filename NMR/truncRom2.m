function [POD] = truncRom2(POD,train,varargin)
    %TRAINROM1 Summary of this function goes here
    %   Detailed explanation goes here
    if length(varargin) >= 1
        fig = varargin{1};
    else
        varargin = {};
    end
    fig.Name = 'ROM-2 Train';
    R = [];
    J = [];
    for i = 1:length(train.sol)
        R = [R,train.sol{i}.Rx];
        J = [J,train.sol{i}.Jnz];
    end
    POD.Rx = makePOD(R,false,false,0.99);
    POD.Jnz = makePOD(J,false,false,0.99);
    
    POD.Rx.N = max(POD.Rx.N,POD.X.N);
    POD.Jnz.N = max(POD.Jnz.N,POD.X.N);
    
    nR = 0.5;
    nR = round(POD.Rx.N * (1 + nR));
    nR = nR/size(POD.Rx.U,1);
    POD.Rx = makeOblique(nR,POD.Rx);
    disp(norm(applyOblique(R,POD.Rx) - R));

    nJ = 0.5;
    nJ = round(POD.Jnz.N * (1 + nJ));
    nJ = nJ/size(POD.Jnz.U,1);
    POD.Jnz = makeOblique(nJ,POD.Jnz);
    disp(norm(applyOblique(J,POD.Jnz) - J));
end

