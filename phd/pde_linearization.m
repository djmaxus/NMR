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

dudt = (u - u0) ./dt;

residual = dudt(2:end-1) + dfdx;
residual = [u(1) - ul;...
    residual;...
    u(end) - ur];

Nu = Nu_all - 2;

j_dudt = 1/dt;

dgdu = 2.*u;
dgdu_diag = sparse(1:length(dgdu),1:length(dgdu),dgdu);
dgdu_upstream = dgdu_diag(upstream_index,:);

sp_dudu = ones(Nf,1);
sp_dududx = conv2(sp_dudu,ddx');
sp_dgdu_upstream = [diag(dgdu_upstream,0),diag(dgdu_upstream,1)];
sp_dfdu = -kappa .* (sp_dgdu_upstream.*dudx + g_upstream .* sp_dududx);
sp_dfdu = [[sp_dfdu(:,1);0],[0;sp_dfdu(:,2)]];

sp_dfdudx = conv2(sp_dfdu,ddx');
% sp_dfdudx_final = [sp_dfdudx(1:Nu,1),sp_dfdudx((1:Nu)+1,2),sp_dfdudx((1:Nu)+2,3)];

flat_jacobian.i = repmat((1:Nu)+1,1,3);
flat_jacobian.j = [1:Nu,(1:Nu)+1,(1:Nu)+2];
flat_jacobian.v = [sp_dfdudx(1:Nu,1);sp_dfdudx((1:Nu)+1,2);sp_dfdudx((1:Nu)+2,3)];

indices_main_diag = Nu + (1:Nu);
flat_jacobian.v(indices_main_diag)=flat_jacobian.v(indices_main_diag)+j_dudt;

jacobian = sparse(...
    flat_jacobian.i,...
    flat_jacobian.j,...
    flat_jacobian.v,...
    Nu_all,Nu_all);

jacobian(1,1) = 1;
jacobian(end,end) = 1;
end
