function [D] = setDiff(Nx,param)
    %SETDIFF Summary of this function goes here
    %   Detailed explanation goes here
    % distribution
    A = param.A;
    Np = param.P;
    S = param.S;
    D = A * sin(Np*2*pi*(1:Nx+1)'./(Nx+1)) + 1;
    % scale
    D = D * S;
end

