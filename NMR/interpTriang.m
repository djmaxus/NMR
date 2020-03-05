function [Fq] = interpTriang(Xq,xV,F)
    %INTERPTRIANG Summary of this function goes here
    %   Detailed explanation goes here
    pClose = dsearchn(xV,Xq);
    dF = F - F(pClose,:);
    dX = xV - xV(pClose,:);
    gradF = lscov(dX,dF);
    Fq = F(pClose,:) + (Xq -xV(pClose,:))*gradF;
end



