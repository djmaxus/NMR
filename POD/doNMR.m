function [POD] = doNMR(FOM,fX,fR,fS)
    %DONMR Summary of this function goes here
    %   Detailed explanation goes here
    %% X POD
    POD.X = makePOD(FOM.X,true,false,fX);
    %% Residual POD
    POD.Rx = makePOD(FOM.Rx,true,false,fR);
    POD.Rx = makeOblique(fS,POD.Rx);
    %% Jacobian
    POD.Jnz = makePOD(FOM.Jnz,true,true,fR);
    POD.Jnz = makeOblique(fS,POD.Jnz);
end
