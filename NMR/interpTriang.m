function [Fq] = interpTriang(Xq,X,grid,F)
    %INTERPTRIANG Summary of this function goes here
    %   Detailed explanation goes here
    for j = 1:size(grid,1)
        nV = grid(j,:);
        xV = X(nV,:);
        in = inSimplex(Xq,xV); 
        if in
            break;
        end
    end
    pClose = dsearchn(xV,Xq);
    F = F(nV,:);
    dF = F - F(pClose,:);
    dX = X(nV,:) - xV(pClose,:);
    gradF = lscov(dX,dF);
    Fq = F(pClose,:) + (Xq -xV(pClose,:))*gradF;
end

function out = inSimplex(Xq,Xp)
    out = true;
    for i = 1:size(Xp,1)
        Xs = Xp;
        Xs(i,:) = [];
        xp = Xp(i,:);
        xm = mean(Xs,1);
        dq = (Xq - xm);
        dp = (xp - xm);
        out = out .* (dp * dq' > 0);
    end
end

