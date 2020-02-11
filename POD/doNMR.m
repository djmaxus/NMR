function [POD] = doNMR(FOM,fX,fR,fS)
    %DONMR Summary of this function goes here
    %   Detailed explanation goes here
    %% X POD
    POD.X = makePOD(FOM.X,fX);
    %% Residual POD
    % flux
    POD.W = makePOD(FOM.W,fR);
    POD.W = makeOblique(fS,POD.W);
end

