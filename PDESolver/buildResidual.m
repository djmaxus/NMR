function [argout]=buildResidual(X,dt,D,param,PODr)
    Nx = size(X,1);
    % time derivative
    Rt = (X(:,2)-X(:,1))./dt;
    
    % spatial part
    W = buildFlowRes(X,D,param.bc,PODr);
    Rx = W(2:Nx+1)-W(1:Nx);
    
    % source term
    Rs = param.q.val;
    R = Rt + Rx + Rs;

    argout.value = R;
    argout.W = W;
end
