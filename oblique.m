function [O] = oblique(P,basis)
    %OBLIQUE Oblique vector projection
    %   Detailed explanation goes here
    %% operator to inverse
    M = P' * basis;
    % SVD
    [U,S,V] = svd(M,'econ');
    % inverse S matrix
    diagS = 1./ diag(S);
    Sx = zeros(size(S'));
    for i = 1:min(size(Sx))
        Sx(i,i) = diagS(i);
    end
    % generalized matrix inversion
    Mx = V * Sx * U';
    O = sparse(basis * Mx * P');
end
