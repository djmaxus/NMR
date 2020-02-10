function [O] = oblique(P,basis)
    %OBLIQUE Oblique vector projection
    %   Detailed explanation goes here
    M = P' * basis;
    [U,S,V] = svd(M);
    N = min(size(S));
    diagSx = diag(S(1:N,1:N));
    diagSx = eye(N)./ diagSx;
    Sx = zeros(size(S'));
    Sx(1:N,1:N)= diagSx;
    Mx = V * Sx * U';
    O = sparse(Mx * P');
end

