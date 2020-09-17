function [Y] = inv_svd(X)
    %INV_SVD Summary of this function goes here
    %   Detailed explanation goes here
    [U,S,V] = svd(X,'econ');
    Y = V * diag(diag(S).^(-1))*U';
end

