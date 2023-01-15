function [O, decoder] = oblique_projector(indices,basis)
% OBLIQUE Oblique projection operator
%   Input:
%       P   - mask projection operator
%       basis - Principal components
%   Output:
%       O - oblique projection operator
%   Reference:
%       K. Carlberg, C. Bou-Mosleh, and C. Farhat,
%       “Efficient non-linear model reduction via a least-squares
%       Petrov-Galerkin projection and compressive tensor approximations,
%       ” Int. J. Numer. Methods Eng., vol. 86, no. 2, pp. 155–181, 2011,
%       doi: 10.1002/nme.3050.

basis_projected = basis(indices,:);
O = pinv(basis_projected);
decoder = basis * O;
end
