function [Y] = pinv(X,varargin)
    %INV_SVD Summary of this function goes here
    %   Detailed explanation goes here
    [U,S,V] = svds(X,varargin{1:end});
    Y = V * diag(diag(S).^(-1))*U';
end

