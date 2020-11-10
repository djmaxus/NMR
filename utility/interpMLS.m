function [f] = interpMLS(x,X,F)
    %INTERPMLS Summary of this function goes here
    %   Detailed explanation goes here
    N = length(x);
    w = zeros(size(X,2),1);
    for i = 1:size(X,2)
        w(i) = norm(X(:,i)-x);
    end
    [~,idx] = sort(w,'ascend');
    add = 0;
    X = X(:,idx(1:N+1+add));
    F = F(idx(1:N+1+add),:);
    
    x0 = X(:,1);
    f0 = F(1,:);
    
    A = (X - x0)';
    B = F - f0;
    
    A(1,:) = [];
    B(1,:) = [];

    dF = lscov(A,B);
    f = f0 + x'*dF;
end

