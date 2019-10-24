function [R,mse] = reconstruct(F,P,basis)
    O = oblique(P, basis);
    R = O * F;
    mse = sum((R-F).^2,'all')/numel(F);
end