    %OBLIQUE Oblique vector projection
    %   Detailed explanation goes here
    %% operator to inverse
    M = P' * basis;
    % SVD
    [U,S,V] = svd(M,'econ');
    % inverse S matrix
    S = diag(S);
    S = 1./S;
    S = diag(S);
    % generalized matrix inversion
    Mx = V * S * U';
    O = sparse(basis * Mx * P');
end
