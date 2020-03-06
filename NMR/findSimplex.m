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
    out_m = true;
    out = true;
    for i = 1:size(Xp,1)
        Xs = Xp;
        Xs(i,:) = [];
        xp = Xp(i,:);
        xm = mean(Xs,1);
        %% to mean point
        dq = (Xq - xm);
        dp = (xp - xm);
        out_m = out_m .* (dp * dq' > 0);
        %% orthogonal projection
%         A = [Xs;xm];
%         A = [A,-ones(size(A,1),1)];
%         n = A \ zeros(size(A,1),1);
%         Dq = hpdist((Xq-xm)',n);
%         Dp = hpdist((xp-xm)',n);
%         out = out .* (Dp' * Dq >= 0);
    end
%     if out ~= out_m
%         warning('Simplex enclosure mismatch');
%     end
    out = out_m;
end

function out = hpdist(x,n)
    out = diag(n) * x ./ norm(n);
end
