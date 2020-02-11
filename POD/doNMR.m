function [POD] = doNMR(FOM,fX,fR,fS)
    %DONMR Summary of this function goes here
    %   Detailed explanation goes here
    Nx = size(FOM.X,1);
    %% X POD
    POD.X = makePOD(FOM.X,fX);
    %% Residual POD
    % flow sum
    POD.Rx = makePOD(FOM.Rx,fR);
    POD.Rx = makeOblique(fS,Nx,POD.Rx);
    % left flow
    POD.Rxl = makePOD(FOM.Rxl,fR);
    POD.Rxl = makeOblique(fS,Nx,POD.Rxl);
    % right flow
    POD.Rxr = makePOD(FOM.Rxr,fR);
    POD.Rxr = makeOblique(fS,Nx,POD.Rxr);
end

