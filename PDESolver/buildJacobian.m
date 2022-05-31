function [argout] = buildJacobian(X,dt,param,varargin)
    Nx = size(X,1);
    
    % time derivative
    Jt = speye(Nx)./dt;
    
    % spatial part
    evalTime = tic;
        Jf = buildFlowJac(X,param.D,param.bc,varargin{:}); 
    evalTime = toc(evalTime);
    argout.time = evalTime;
    
    Jx = spalloc(Nx,Nx,3*Nx-2);
    Jxb = sparse(Nx,2);
    for i=1:Nx
        % lower diagonal
        J1 = Jf(i,1);
        if i == 1
            Jxb(i,1) = - J1;
        else
            Jx(i,i-1) = -J1;
        end
        % diagonal
        J1 = Jf(i,2);
        J2 = Jf(i+1,1);
        Jx(i,i) = J2 - J1;
        % upper diagonal
        J2 = Jf(i+1,2);
        if i == Nx
            Jxb(i,2) = J2;
        else
            Jx(i,i+1) = J2;
        end
    end
    
    J = Jx + Jt;
    
    % boundary
    Jb = param.bc.jacob(param.bc.state,X,dt,param);
    
    % output
    argout.value = J;
    argout.Jt = Jt;
    argout.Jf = Jf;
    argout.Jx = Jx;
    argout.Jxb = Jxb;
    argout.Jb = Jb;
end
