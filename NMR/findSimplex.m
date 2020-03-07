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
        n = hpn(Xs');
        if sum(isinf(n),'all') + sum(isnan(n),'all') > 0
            warning('singular');
        elseif norm(n) <= 0
            warning('zero length');
        else
            x0 = Xs(1,:);
            Dq = hpdist(Xq',n,x0');
            Dp = hpdist(xp',n,x0');
            out = out .* (Dp * Dq >= 0);
        end
    end
    if out ~= out_m
        warning('Simplex enclosure definition mismatch');
    end
%     out = out_m;
end

function out = hpdist(x,n,x0)
    out = n' * (x-x0) ./ norm(n);
end

function [n] = hpn(X)
    x0 = X(:,1);
    A = (X - x0)';
    A(1,:) = 2;
%     n = zeros(size(x0));
%     n(1) = A(2,2);
%     n(2) = -A(2,1);
    n = rand(size(x0));
    for k = 1:100
        n = n ./ norm(n);
        R = A*n;
        R(1) = norm(n) - 1;
        delta = norm(R);
        if delta < 1e-8
            break;
        end
        J = A;
        J(1,:) = A(1,:) .* n';
        n = n + J'*J \ (-J'*R);  
    end
end
