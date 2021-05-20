function [problem] = setProblem(Nx,bc,q,D,x0,arg)
    % boundary conditions
    param.bc.first = bc.isKindI(bc.kind);
    param.bc.val = bc.val;
    % source terms
    param.q.val = spalloc(Nx,1,1);
    param.q.pos = round(q.pos*Nx);
    param.q.val(param.q.pos) = q.val;
    % diffusion
    param.D.val = setDiff(Nx,D);
    param.D.model = @(u)([max(u.^2,0)+1,max(2*u,0)]);
    
    problem.param = param;
    %% initial condition
    problem.x0 = x0;
    %% equations
    problem.residual = @(x,dt,POD) buildResidual(arg(x),dt,problem.param,POD);
    problem.jacobian = @(x,dt,POD) buildJacobian(arg(x),dt,problem.param,POD);
end
