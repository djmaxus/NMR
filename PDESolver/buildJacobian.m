function [argout] = buildJacobian(X,dt,param,varargin)
    Nx = size(X,1);
    
    % time derivative
    Jt = speye(Nx)./dt;
    
    % spatial part
    [Jf] = buildFlowJac(X,param.D,param.bc,varargin{:});
    Jx = spalloc(Nx,Nx,3*Nx-2);
    for i=1:Nx
        if i ~=1
            J1 = Jf(i,1);
            Jx(i,i-1) = -J1;
        end
        J1 = Jf(i,2);
        J2 = Jf(i+1,1);
        Jx(i,i) = J2 - J1;
        if i ~= Nx
            J2 = Jf(i+1,2);
            Jx(i,i+1) = J2;
        end
    end
    
    J = Jx + Jt;
    
    argout.value = J;
    argout.Jf = Jf;
end
