function [x0] = setInit(Nx,a,p,b)
    %SETINIT Summary of this function goes here
    %   Detailed explanation goes here
    x0 = a * sin(p*2*pi*(1:Nx)'./Nx)+b;
end

