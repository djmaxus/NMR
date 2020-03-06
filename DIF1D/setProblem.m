%% Diffusion
% amplitude
aD = 0.9;
% period number
pD = 1;
% distribution
D = aD * sin(pD*2*pi*(1:Nx+1)'./(Nx+1)) + 1;
% scale
D = D * 1.0e-2;

% boundary condition type
bc.first = false;
