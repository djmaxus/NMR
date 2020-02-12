function [POD] = doNMR(FOM,fX,fR,fS)
    %DONMR Summary of this function goes here
    %   Detailed explanation goes here
    %% X POD
    POD.X = makePOD(FOM.X,true,false,fX);
    %% Residual POD
    POD.Rx = makePOD(FOM.Rx,false,false,fR);
    POD.Rx = makeOblique(fS,POD.Rx);
end
