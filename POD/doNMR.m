function [POD] = doNMR(FOM,fX,fR,fS)
    %DONMR Summary of this function goes here
    %   Detailed explanation goes here
    %% X POD
    POD.X = makePOD(FOM.X,true,true,fX);
    %% Residual POD
    % flux
    POD.W = makePOD(FOM.W,false,false,fR);
    POD.W = makeOblique(fS,POD.W);
end
