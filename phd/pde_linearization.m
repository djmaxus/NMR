function [residual, jacobian] = pde_linearization(u0,dt,u,dx,ul,ur,kappa)
ddx = [-1;1]/dx;
dudx = conv(u,flip(ddx),'valid');

Nu_all = length(u);
Nf = Nu_all-1;
sign_dudx = dudx > 0;
upstream_index = (1:Nf)' + sign_dudx;
g = u.^2 + 1;
g_upstream = g(upstream_index);

f = -kappa .* g_upstream .* dudx;
dfdx = conv(f,flip(ddx),'valid');

ddt = [-1; 1]./dt;
dudt = [u0,u] * ddt;

residual = dudt(2:end-1) + dfdx;
residual = [u(1) - ul;...
    residual;...
    u(end) - ur];

Nu = Nu_all - 2;

j_dudt = eye(Nu)./dt;

dgdu = 2.*u;
dgdu = sparse(1:Nu_all,1:Nu_all,dgdu);
dgdu_upstream = dgdu(upstream_index,:);

sp_dudu = ones(Nf,1);
sp_dududx = conv2(sp_dudu,ddx');
sp_dgdu_upstream = spdiags(dgdu_upstream);
sp_dfdu = -kappa .* (sp_dgdu_upstream.*dudx + g_upstream .* sp_dududx);
sp_dfdu = [[sp_dfdu(:,1);0],[0;sp_dfdu(:,2)]];

sp_dfdudx = conv2(sp_dfdu,ddx');
sp_dfdudx = [sp_dfdudx(1:Nu,1),sp_dfdudx((1:Nu)+1,2),sp_dfdudx((1:Nu)+2,3)];

jacobian = sparse(repmat(1:Nu,1,3),[1:Nu,(1:Nu)+1,(1:Nu)+2],...
    [sp_dfdudx(:,1),sp_dfdudx(:,2),sp_dfdudx(:,3)],...
    Nu_all,Nu_all);

jacobian(2:end-1,2:end-1) = jacobian(2:end-1,2:end-1) + j_dudt;
jacobian(1,1) = 1;
jacobian(end,end) = 1;
end

