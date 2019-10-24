function [O] = oblique(P,basis)
    %OBLIQUE Oblique vector projection
    %   Detailed explanation goes here
    M = P' * basis;
    [U,S,V] = svd(M);
    Sx = zeros(size(S'));
    diagSx = diag(S);
    diagSx = eye(length(diagSx))./ diagSx;
    Sx(1:end,1:length(diagSx))= diagSx;
    Mx = V * Sx * U';
    O = sparse(basis * Mx * P');
end

