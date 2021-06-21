function [x0] = genInit(Nx,amp,per,dev)
    %GENINIT Summary of this function goes here
    %   Detailed explanation goes here
    x0 = amp * sin(per * 2 * pi * (1:Nx)'./ Nx) + dev;
end

