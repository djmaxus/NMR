function [nV,xV] = findSimplex(Xq,X)
    %FINDSIMPLEX Summary of this function goes here
    %   Detailed explanation goes here
    grid = delaunayn(X);
    nV = [];
    xV = [];
    for j = 1:size(grid,1)
        nV = grid(j,:);
        xV = X(nV,:);
        in = inSimplex(Xq,xV); 
        if in
            break;
        end
     end
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
