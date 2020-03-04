%% Diffusion
% amplitude
aD = 0.5;
% period number
pD = 2;
% distribution
D = aD * sin(pD*2*pi*(1:Nx+1)'./(Nx+1)) + 1;
% scale
D = D * 1e-2;

% boundary condition type
bc.first = true;
