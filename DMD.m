function [A, Phi, Lambda, b] = DMD(X1,X2,r)
    % partial pseudo-inverse
    [U,S,V] = svd(X1,'econ');
    Ur = U(:,1:r);
    Sr = S(1:r,1:r);
    Sr_inv = diag(diag(1./Sr));
    Vr = V(:,1:r);
    X1_inv = Vr * Sr_inv * Ur';
    
    % Linear propagator
    A = X2*X1_inv;
    
    % DMD spectrum
    Ar = Ur' * X2 * Vr / Sr;
    [W,Lambda] = eig(Ar);
    Phi = X2*(Vr/Sr)*W;
    
    % initial state
    alpha1 = Sr*Vr(1,:)';
    b = (W*Lambda)\alpha1;
end
