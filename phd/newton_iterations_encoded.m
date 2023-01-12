function [sol,iter] = newton_iterations_encoded(calc_linear,initial,tolerance,max_iters, linear_encoder)
check_class(linear_encoder,'EncoderLinear');
jacobian_decoder = linear_encoder.decode_jacobian();

% set initial approximation
current = initial;
delta = Inf(size(current));

% launch loop
for iter = 0:max_iters
    [residual, jacobian] = calc_linear(current);
    
    jac_encoded = jacobian * jacobian_decoder;
    res_ls = jac_encoded' * residual;
    jac_ls = jac_encoded' * jac_encoded;
    
    % calculate errors
    err_abs = norm(residual) / sqrt(length(residual));
    err_abs_ls = norm(res_ls) / sqrt(length(res_ls));
    
    
    % memorize initials
    if iter == 0
        err_initial = err_abs;
        err_initial_ls = err_abs_ls;
    end
    
    % calculate relative errors
    err_rel = err_abs / err_initial;
    err_rel_ls = err_abs_ls / err_initial_ls;
    
    err_delta = norm(delta) / sqrt(length(delta));
    
    % Newton convergence
    is_conv_abs = err_abs < tolerance;
    is_conv_rel = err_rel < tolerance^2;
    is_conv_delta = err_delta < tolerance^3;
    
    is_conv_abs_ls = err_abs_ls < tolerance;
    is_conv_rel_ls = err_rel_ls < tolerance^2;
    
    is_conv_any = is_conv_abs || is_conv_rel || is_conv_delta || is_conv_abs_ls || is_conv_rel_ls;
    
    % is_conv_any
    if and(is_conv_any, iter == max_iters)
        break;
    end
    
    % solve linear system
    delta_encoded = jac_ls \ res_ls;
    
    delta = jacobian_decoder * delta_encoded;
    
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
