function [U,timesteps] = ode_solve(calc_linearization,u0,Nt,use_wb)
if use_wb
    wb = waitbar(0, 'Solving ODE');
end

t = 0;
U = u0;
timesteps = t;

dt0 = 1/Nt;
while t < 1
    dt_max = dt0;
    dt = min([dt0, 1-t, dt_max]);
    if dt < 1e-16
        break;
    end
    
    [u, iters] = ode_step(calc_linearization,U(:,end),dt);
    
    is_converged = ~isempty(u);
    
    if ~is_converged
        dt_min = dt0;
        if (dt < dt_min)
            warning('Not converged');
            break;
        end
        timestep_dec = 0.25;
        dt = dt * timestep_dec;
        continue;
    end
    
    t = t+dt;
    
    U = [U,u];
    timesteps = [timesteps,t];
    
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

function [u,iters] = ode_step(calc_linearization,u0,dt)
calc_linear_for_newton = @(u) calc_linearization(u0,dt,u);
[u, iters] = newton_iterations(calc_linear_for_newton,u0,1e-3,10);
end
