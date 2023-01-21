function [sol,iter, data] = newton_iterations_encoded(calc_linear,initial,tolerance,max_iters, linear_encoder)
arguments
    calc_linear (1,1)
    initial (:,1) {mustBeReal}
    tolerance (1,1) {mustBeReal, mustBeNonnegative}
    max_iters (1,1) {mustBeInteger, mustBePositive}
    linear_encoder {mustBeA(linear_encoder,'EncoderLinear'), mustBeScalarOrEmpty}
end

if isempty(linear_encoder)
    error("empty linear_encoder");
end

jacobian_decoder = linear_encoder.decode_jacobian();

% set initial approximation
current = initial;
delta_encoded = Inf(size(linear_encoder.dim_latent()));

data = cell(1,max_iters + 1);

% launch loop
for iter = 0:max_iters
    [residual, jacobian, data_iter] = calc_linear(current);
    
    data{iter+1} = data_iter;
    
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
    
    err_delta = norm(delta_encoded) / sqrt(length(delta_encoded));
    
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
