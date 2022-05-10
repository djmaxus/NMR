function [res, jac, G] = gen_both(U,dt,params)
%GEN_BOTH Summary of this function goes here
%   Detailed explanation goes here

%% both
dx = params.dx;
dx2 = dx^2;

Dw  = params.Dw;

K   = params.K;

Udx = params.Dxe * U(:,2); % dU/dx

sX = round(0.5 * (params.Xup + sign(Udx))); % upstream index

G = params.G(U(:,2)); % nonlinear measurement

GG = sparse(diag(G(sX,1)));

%% resid use
q   = @() params.q;

Udt = @() U(2:end-1,:) * params.Dt; % time derivative

W = @() -K * GG * Udx;

%% inner residual
res = dx2 * Udt() + dt * Dw * W() - dt * dx * q();

%% jac use
Dx = @() params.Dx;

Jt = @() params.Jt;

dG = sparse(diag(G(:,2)));

Jg = @() dG(sX,2:end-1); % upstream approximation

Jw = @() -K * (GG * Dx() + Jg() .* Udx); % flow jacobian

%% inner Jacobian
jac = dx2 * Jt() + dt * Dw * Jw(); % inner Jacobian

end

