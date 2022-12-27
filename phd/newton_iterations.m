function [sol] = newton_iterations(calc_linear,initial,tolerance,max_iters)

% set initial approximation
current = initial;
delta = Inf(size(current));

% launch loop
for i = 0:max_iters
    [residual, jacobian] = calc_linear(current);
    
    % calculate errors
    err_abs = norm(residual) / sqrt(length(residual));
    
    % memorize initials
    if i == 0
        err_initial = err_abs;
    end
    
    % calculate relative errors
    err_rel = err_abs / err_initial;
    
    err_delta = norm(delta) / sqrt(length(delta));
    
    % Newton convergence
    is_conv_abs = err_abs < tolerance;
    is_conv_rel = err_rel < tolerance^2;
    is_conv_delta = err_delta < tolerance^3;
    
    is_conv_any = is_conv_abs || is_conv_rel || is_conv_delta;
    
    % is_conv_any
    if is_conv_any
        break;
    end
    
    % solve linear system
    delta = jacobian \ residual;
    current = current - delta;
end

if is_conv_any
    sol = current;
else
    sol = [];
end

end

function [delta_j] = verify_jacobian(calc_linearization,u0)
[residual0, jacobian0] = calc_linearization(u0);
jacobian = zeros(size(jacobian0));
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
