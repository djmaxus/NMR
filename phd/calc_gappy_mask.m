function [gappy_mask] = calc_gappy_mask(basis,n)
% GREEDY Greedy sampling algorithm
%   Input:
%       basis - Principal components
%       n     - number of sampling points
%   Output:
%       idx - range of sampling points
%       P   - mask projection operator
%   Reference:
%       K. Carlberg, C. Bou-Mosleh, and C. Farhat,
%       “Efficient non-linear model reduction via a least-squares
%       Petrov-Galerkin projection and compressive tensor approximations,
%       ” Int. J. Numer. Methods Eng., vol. 86, no. 2, pp. 155–181, 2011,
%       doi: 10.1002/nme.3050.

arguments
    basis (:,:) {double}
    n (1,1) {mustBeInteger,mustBePositive}
end

gappy_mask = zeros(n,1);

[~, dim_latent] = size(basis);

current_mask_size = 0;

modulus = mod(n,dim_latent);
ratio = floor(n/dim_latent);

% iterate over all vectors, but not less than the mask size
for i = 1:min([dim_latent,n])

    current_mask = gappy_mask(1:current_mask_size);
    basis_projected = basis(current_mask,1:i);
    basis_projected_prev = basis_projected(:,1:(end-1));
    
    M = pinv(basis_projected_prev);
    
    basis_prev = basis(:,1:(i-1));
    basis_projected_curr = basis_projected(:,end);
    
    distinct_error = basis(:,i) - basis_prev * M * basis_projected_curr;
    
    distinct_error = abs(distinct_error);
    [~,indices_error] = sort(distinct_error,'descend');
   
    points_per_vector = ratio + (i < (modulus + 1)); 
    for j = 1:points_per_vector
        % find first **new** index for measurement location
        for index_error = indices_error' 
            does_index_repeat = sum(index_error == gappy_mask) > 0;
            if ~does_index_repeat
                break;
            end
        end
        current_mask_size = current_mask_size+1;
        gappy_mask(current_mask_size) = index_error;
    end
end

end
