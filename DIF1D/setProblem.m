function [param] = setProblem(Nx,bc,q,D,x0)
    %% boundary conditions
    param.bc.first = bc.isKindI(bc.kind);
    param.bc.val = [0,0];
    %% source terms
    param.q = spalloc(Nx,1,1);
    q.pos = round(q.pos*Nx);
    param.q(q.pos) = q.val;
    %% diffusion
    param.D = setDiff(Nx,D);
    %% initial condition
    param.x0 = x0;
end
