function [argout]=buildResidual(X,dt,param,varargin)

    Nx = size(X,1);
    % time derivative
    Rt = (X(:,2)-X(:,1))./dt;
    
    % spatial part
    W = buildFlowRes(X,param.D,param.bc,varargin{:});
    Rx = W(2:Nx+1)-W(1:Nx);
    
    % source term
    Rs = param.q.val;
    R = Rt + Rx - Rs;
    
    % boundary law
    Rb = param.bc.resid(param.bc.state,X,dt,param);
    
    % output
    argout.value = R;
    argout.W = W;
end
