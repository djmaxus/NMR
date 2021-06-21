function [bc] = genBC(kind,lim)
    %BC Summary of this function goes here
    %   Detailed explanation goes here
    bc.kind = kind;
    bc.val = rand(2,1) * (lim(2)-lim(1)) + lim(1);
    bc = genBoundModel(bc,bc.val);
end

function [bc] = genBoundModel(bc,b0)
    bc.state = [b0,b0];
    bc.resid = @(b,u,dt,param) (b(:,2)-b(:,1));
    bc.jacob = @(b,u,dt,param)([ sparse(size(b,1),size(u,1)), speye(length(b)) ]);
end

