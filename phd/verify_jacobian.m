function [delta_j] = verify_jacobian(calc_linearization,u0)
[residual0, jacobian0] = calc_linearization(u0);
jacobian = sparse(size(jacobian0));
scale = 1e-4;
delta_u = (rand(size(u0))*scale).*abs(u0);
delta_u(abs(delta_u)==0) = scale;
for i = 1:length(delta_u)
    u = u0 + sparse(i,1,delta_u(i),length(u0),1);
    [residual, ~] = calc_linearization(u);
    jacobian(:,i) = (residual - residual0)./delta_u(i);
end
delta_j = jacobian0 - jacobian;
end
