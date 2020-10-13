function [param] = setProblem(Nx,bc,q,D,x0)
    %% boundary conditions
    param.bc.first = bc.isKindI(bc.kind);
    param.bc.val = bc.val;
    %% source terms
    param.q.val = spalloc(Nx,1,1);
    param.q.pos = round(q.pos*Nx);
    param.q.val(param.q.pos) = q.val;
    %% diffusion
    param.D = setDiff(Nx,D);
    %% initial condition
    param.x0 = x0;
end
