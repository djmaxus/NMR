function [sol] = ode_solve_encoded(calc_linearization,u0,Nt,use_wb,encoder)

if use_wb
    wb = waitbar(0, 'Solving ODE');
end

sol.U = zeros(length(u0),Nt+1);
sol.U(:,1) = u0;

t = 0;
sol.timesteps = zeros(1,Nt+1);
sol.timesteps(1) = t;

dt0 = 1/Nt;
dt_max = dt0;
dt = dt0;
dt_min = dt0;
timestep_dec = 0.25;

n = 1;
while t < 1

    dt = min([dt, 1-t, dt_max, dt0]);
    if dt < 1e-16
        break;
    end
    
    u_prev = sol.U(:,n);
    calc_linear_for_newton = @(u) calc_linearization(u_prev,dt,u);
    [u, iters] = newton_iterations_encoded(calc_linear_for_newton,u_prev,1e-3,10,encoder);

    is_converged = ~isempty(u);

    if ~is_converged
        if (dt < dt_min)
            warning('Not converged');
            break;
        end

        dt = dt * timestep_dec;
        continue;
    end

    t = t+dt;
    n = n + 1;

    sol.U(:,n) = u;
    sol.timesteps(n) = t;

    if iters < 3
        timestep_inc = 1;
        dt =  dt * timestep_inc;
    end

    if use_wb
        waitbar(t,wb);
    end
end
if use_wb
    delete(wb);
end
end
